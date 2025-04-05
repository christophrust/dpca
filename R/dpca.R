#' @title Dynamic Principal Component Analysis and Dynamic Factor Model Estimation
#'
#' @description \code{dpca} is used to estimate a Generalized Dynamic Factor Model
#' in the spirit of Forni et al. (2000) and Forni & Lippi (2001) via dynamic principal
#' components analysis (DPCA) a la Brillinger (2001). The number of principal components
#' can be chosen in a data-driven way using the method suggested by Hallin & Liska (2007).
#'
#' \code{dpca} has a cousin \code{\link{spca}} for estimating static principal components as
#' they are used in the dynamic factor model literature around Stock & Watson (2001),
#' also making the Hallin & Liska method for selecting the number of factors available.
#'
#' @param x Input data supplied either as a matrix (rows correspond to cross-sectional
#' units and columnts to observations in the time domain) or a multivariate object of
#' class \code{\link[stats]{ts}} or \code{\link[zoo]{zoo}}.
#'
#' @param q Number of dynamic factors. If no value is supplied for \code{q},
#' it is chosed data-driven by using the criterion of Hallin & Liska (2007).
#'
#' @param freqs Numeric vector in \eqn{[-\pi, \pi]} giving the frequencies where the
#' spectral density is evaluated.
#'
#' @param bandwidth Single integer, giving the width of the
#' lag window estimator. If unspecified, 0.75 times the square root of the number
#' of time observations is used as default.
#'
#' @param weights Kernel used for the lag window estimation of spectrum.
#'
#' @param qsel_crit Criterion to select the number of factors using the
#' Hallin & Liska (2007, JASA) method. Either \code{"IC1"} or \code{"IC2"}.
#'
#' @param q_max Maximum numer of dynamic factors considered in the data-driven
#' selection.
#'
#' @param n_path Integer vector specifying which (nested) subsets of the
#' cross section are used in the Hallin & Liska procedure. If unspecified,
#' a regular sequence of length \code{20} from \code{n/2} to \code{n} is used.
#'
#' @param penalties Evaluated values of the penalty function at
#' each value of the n_path. In case this is missing, the penalty \eqn{p_1(n,T)}
#' suggested in Hallin & Liska (2007, p. 611) is used:
#'
#' \deqn{
#'   (h^{-2} + \sqrt{h/T} + 1/n) \cdot \log{\operatorname{min}(n, h^2, \sqrt{T/h})}
#' }{%
#'   (bandwidth^(-2) + sqrt(bandwidth / ncol(x)) + 1 / n_path)
#'    * log(pmin(n_path, bandwidth^2, sqrt(ncol(x) / bandwidth)))
#' }
#'
#' @param penalty_scales Tuning values for the penalty scaling parameter
#' \eqn{c} over which the \code{q}-path is optimized to stability.
#'
#'
#' @return A object of class "dpca" wrapping a list with the entries
#' \itemize{
#'   \item \code{xmean}: a vector holding the mean of each cross-sectional unit
#'   \item \code{spectrum}: the estimated spectral density without the
#'     normalization by a factor of 2 \code{pi}. Therefore, a white noise
#'     process would have a spectrum equal to one a.e.
#'   \item \code{eig}: eigen decomposition of the spectral density
#'   \item \code{filter}: a list holding the filter coefficients for the filter
#'     returning input and dynamic common component.
#'   \item \code{input}: the input series.
#'   \item \code{dcc}: (dynamic) common component
#'   \item \code{dic}: (dynamic) idiosyncratic component
#'   \item \code{HL_select}: results of the selection methodology of Hallin & Liska (2007),
#' }
#' also see \code{\link{select_r}}.
#'
#' @examples
#' data(fredmd)
#' fredmd <- scale(fredmd)
#'
#' freqs <- -50:50 / 50 * pi
#' res <- dpca::dpca(fredmd, freqs = freqs, q_max = 10)
#'
#' ## eigenvalues
#' matplot(x = freqs, y = t(res$eig$values), type = "l")
#'
#' ## q selection
#' cat(sprintf("Number of selected dynamic components: %s\n", res$HL_select$q))
#'
#' ## sample variability of the criterion S^2_C (Hallin & Liska 2007, equation 10)
#' plot(
#'   x = res$HL_select$penalty_scales, y = res$HL_select$sample_var, type = "l",
#'   col = "red"
#' )
#' par(new = TRUE)
#' plot(
#'   x = res$HL_select$penalty_scales, y = res$HL_select$q_path,
#'   type = "l", axes = FALSE, bty = "n", xlab = "", ylab = "",
#'   col = "blue"
#' )
#' axis(4)
#' mtext("q_path", side = 4)
#'
#' @references
#' Brillinger, D. R. (2001). Time Series: Data Analysis and Theory. SIAM.
#'
#' Forni, M., Hallin, M., Lippi, M., & Reichlin, L. (2000). The generalized
#' dynamic-factor model: Identification and estimation. Review of Economics
#' and statistics, 82(4), 540-554.
#'
#' Forni, M., & Lippi, M. (2001). The Generalized Dynamic Factor Model:
#' Representation Theory. Econometric Theory, 17(6), 1113-1141.
#'
#' Hallin, M. and Liska, R. (2007). Determining the Number of Factors in
#' the General Dynamic Factor Model. Journal of the American Statistical
#' Association, 102(478).
#'
#' @useDynLib dpca
#' @importFrom stats is.ts
#' @export
dpca <- function(
    x,
    q,
    freqs = -bandwidth:bandwidth / bandwidth * pi,
    bandwidth = floor(0.75 * sqrt(tx)),
    weights = c(
      "bartlett", "trunc", "tukey", "parzen",
      "bohman", "daniell", "parzen_cogburn_davis"
    ),
    qsel_crit = c("IC1", "IC2"),
    q_max = 15,
    n_path = floor(seq(nx / 2, nx, nx / 20)),
    penalties = (
      bandwidth^(-2) + sqrt(bandwidth / tx) + 1 / n_path
    ) * log(pmin(n_path, bandwidth^2, sqrt(tx / bandwidth))),
    penalty_scales = seq(0, 2, by = 0.01)) {
  if (length(weights) > 1) {
    weights <- "bartlett"
  }

  ## our internal object is n times t
  x <- if (is.ts(x) || "zoo" %in% class(x)) {
    t(x)
  } else if (is.matrix(x)) {
    x
  } else {
    stop("x must either a \"ts\" or \"zoo\" object or a matrix!")
  }
  nx <- nrow(x)
  tx <- ncol(x)
  cl <- match.call()

  if (!missing(q) && (length(q) > 1 || floor(abs(q)) != q)) {
    stop("\"q\" has to be a single positive integer!")
  }

  if (length(bandwidth) > 1 || floor(abs(bandwidth)) != bandwidth) {
    stop("\"bandwidth\" has to be a single positive integer!")
  }

  if (!is.numeric(freqs) || any(abs(freqs) > pi)) {
    stop("\"freqs\" must be a numeric vector with values in [-pi, pi]!")
  }

  if (!isTRUE(all.equal(freqs, rev(-freqs)))) {
    warning(
      paste0(
        "\"freqs\" has to be a symmetric sequence around zero!",
        "Making symmetric using only the positive elements!",
        collapse = "\n"
      )
    )
    pfreqs <- freqs[freqs > 0]
    freqs <- c(rev(-pfreqs), 0, pfreqs)
  }

  ## centering
  mx <- rowMeans(x)
  x <- x - mx

  wghts <- get(paste0("weights_", weights))(-bandwidth:bandwidth / bandwidth)

  select_q <-
    if (!missing(q)) {
      0L
    } else if (qsel_crit[1] == "IC1") {
      1L
    } else if (qsel_crit[1] == "IC2") {
      2L
    }

  mode(x) <- "numeric"
  res <- .Call(
    "R_dpca",
    x,
    if (missing(q)) 0L else as.integer(q),
    as.numeric(freqs),
    as.integer(bandwidth),
    .Machine$double.eps,
    as.numeric(wghts),
    as.integer(q_max),
    select_q,
    as.integer(n_path),
    as.integer(rep(tx, length(n_path))),
    as.numeric(penalties),
    as.numeric(penalty_scales),
    PACKAGE = "dpca"
  )
  ## add cross-sectional means
  res$xmean <- mx
  if (select_q) {
    res$HL_select$penalty_scales <- penalty_scales
  }

  res$call <- cl
  res$freqs <- freqs
  class(res) <- "dpca"
  res
}
