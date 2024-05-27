#' Dynamic Principal Component Analysis and Large Dynamic Factor Model Estimation
#'
#' @param x Input data supplied either as a matrix (rows correspond to cross-sectional
#' units and columnts to observations in the time domain) or a multivariate object of
#' class `ts` or `zoo`.
#'
#' @param q Number of dynamic factors. If qsel = TRUE, this is
#' the maximum number of dynamic factors.
#'
#' @param freqs Numeric vector in [-pi, pi] giving the frequencies where the
#' spectral density is evaluated.
#'
#' @param bandwidth Single integer, giving the width of the
#' lag window estimator. If unspecified, the cube root of the number of time
#' observations is used as default.
#'
#' @param weights Kernel used for the lag window estimation of spectrum.
#'
#' @param qsel Logical, if TRUE one of the Hallin & Liska (2007) criteria are used to
#' choose q from the data.
#'
#' @param qsel_crit Criterion to select the number of factors using the
#'   Hallin & Liska (2007, JASA) method. Either "IC1" or "IC2".
#'
#' @param n_path Integer vector specifying which (nested) subsets of the
#' cross section are used in the Hallin & Liska procedure. If unspecified,
#' a regular sequence of length 20 from n/2 to n is used.
#'
#' @param t_path Integer vector specifying Which (nested) subsets of the
#' time domain are used in the Hallin & Liska procedure. If unspecified,
#' a regular sequence of length 20 from T/2 to T is used.
#'
#' @param penalties Evaluated values of the penalty function at
#' each value of the n_path. In case this is missing, the penalies suggested
#' in Hallin & Liska (2007) are used:
#' (bandwidth^(-2) + sqrt(bandwidth/bandwidth) + 1/nrow(x) *
#' log(min(n, bandwidth^2, sqrt(T/bandwidth)))).
#'
#' @param penalty_scales Tuning values for the penalty scaling parameter
#' c over which the q-path is optimized to stability.
#'
#' @return A list with entries
#' \itemize{
#'   \item \code{xmean}: a vector holding the mean of each cross-sectional unit
#'   \item \code{spectrum}: the estimated spectral density
#'   \item \code{eig}: eigen decomposition of the spectral density
#'   \item \code{filter}: a list holding the filter coefficients for the filter
#'     returning input and dynamic common component.
#'   \item \code{input}: the input series.
#'   \item \code{dcc}: (dynamic) common component
#'   \item \code{dic}: (dynamic) idiosyncratic component
#'   \item \code{HL_select}: results of the selection methodology of Hallin & Liska (2007),
#' }
#' see also \code{select_r}.
#'
#' @references Hallin, M. and Liska, R. (2007). Determining the Number of
#' Factors in the General Dynamic Factor Model. Journal of the American
#' Statistical Association, 102(478).
#'
#' @useDynLib dpca
#' @importFrom stats is.ts
#' @export
dpca <- function(
  x,
  q,
  freqs = -100:100 / 100 * pi,
  bandwidth = NULL,
  weights = c(
    "Bartlett", "trunc", "Tukey", "Parzen",
    "Bohman", "Daniell", "ParzenCogburnDavis"
  ),
  qsel = FALSE,
  qsel_crit = c("IC1", "IC2"),
  n_path = NULL,
  t_path = NULL,
  penalties,
  penalty_scales = seq(0, 2, by = 0.01)
) {

  if (length(weights) > 1)
    weights <- "Bartlett"

  x <- if (is.ts(x) || "zoo" %in% class(x)) {
         t(x)
       } else if (is.matrix(x)) {
         x
       } else {
         stop("x must either a \"ts\" or \"zoo\" object or a matrix!")
       }

  ## centering
  mx <- rowMeans(x)
  x <- x - mx

  if (missing(q)) {
    warning(
      "No number of dynamic principal components supplied. Using q = 1..."
    )
    q <- 1L
  }

  if (length(q) > 1 || floor(abs(q)) != q)
    stop("\"q\" has to be a single positive integer!")

  if (is.null(bandwidth)) {
    bandwidth <- floor(ncol(x)^(1 / 3))
  }

  if (length(bandwidth) > 1 || floor(abs(bandwidth)) != bandwidth)
    stop("\"bandwidth\" has to be a single positive integer!")

  if (!is.numeric(freqs) || any(abs(freqs) > pi))
    stop("\"freqs\" must be a numeric vector with values in [-pi, pi]!")

  wghts <- get(paste0("weights.", weights))(-bandwidth:bandwidth / bandwidth)

  if (is.null(n_path)) {
    n_path <- floor(seq(nrow(x) / 2, nrow(x), nrow(x) / 20))
  }
  if (is.null(t_path)) {
    floor(seq(ncol(x) / 2, ncol(x), ncol(x) / 20))
  }

  if (qsel && missing(penalties)) {
    penalties <- (
      bandwidth^(-2) + sqrt(bandwidth / ncol(x)) + 1 / n_path
    ) * log(pmin(n_path, bandwidth^2, sqrt(ncol(x) / bandwidth)))
  } else if (isFALSE(qsel)) {
    penalties <- rep(0, length(n_path))
  }

  select_q <-
    if (isFALSE(qsel)) {
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
    as.integer(q),
    as.numeric(freqs),
    as.integer(bandwidth),
    .Machine$double.eps,
    as.numeric(wghts),
    as.integer(q),
    select_q,
    as.integer(n_path),
    as.integer(t_path),
    as.numeric(penalties),
    as.numeric(penalty_scales),
    PACKAGE = "dpca"
  )

  ## add cross-sectional means
  res$xmean <- mx

  res
}
