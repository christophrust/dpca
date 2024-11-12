#' Static Principal Component Analysis and Large Dynamic Factor Model Estimation
#'
#' @param x Input data supplied either as a matrix (rows correspond to cross-sectional
#' units and columnts to observations in the time domain) or a multivariate object of
#' class \code{\link[stats]{ts}} or \code{\link[zoo]{zoo}}.
#'
#' @param r Number of static factors. If \code{rsel} is set to \code{TRUE}, this is
#' the maximum number of static factors.
#'
#' @param rsel Logical, if \code{TRUE} one of the Hallin & Liska criteria are used to
#' choose \code{r} from the data.
#'
#' @param rsel_crit Criterion to select the number of factors using the
#'   Hallin & Liska (2007, JASA) method. Either \code{"IC1"} or \code{"IC2"}.
#'
#' @param n_path Integer vector specifying which (nested) subsets of the
#' cross section are used in the Hallin & Liska procedure. If unspecified,
#' a regular sequence of length \code{20} from \code{n/2} to \code{n} is used.
#'
#' @param penalty_scales Tuning values for the penalty scaling parameter
#' \eqn{c} over which the \code{q}-path is optimized to stability.
#'
#' @return An object of class "spca", wrapping a list with the entries
#' \itemize{
#'   \item \code{xmean}: a vector holding the mean of each cross-sectional unit
#'   \item \code{cov}: variance-covariance-matrix of \code{x}
#'   \item \code{eig}: eigen decomposition of \code{cov}
#'   \item \code{factors}: an \eqn{r} times \eqn{T} dimensional matrix with the computed factors
#'   \item \code{cc}: (static) common component
#'   \item \code{ic}: (static) idiosyncratic component
#'   \item \code{HL_select}: results of the selection methodology of Hallin & Liska (2007),
#' }
#' see also \code{\link{select_r}}.
#'
#' @importFrom stats is.ts
#' @export
spca <- function(
    x,
    r,
    rsel = FALSE,
    rsel_crit = c("IC1", "IC2", "IC3"),
    n_path = NULL,
    penalty_scales = seq(0, 2, by = 0.01)) {
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

  if (missing(r)) {
    warning(
      "No number of dynamic principal components supplied. Using r = 1..."
    )
    r <- 1L
  }

  if (length(r) > 1 || floor(abs(r)) != r) {
    stop("\"r\" has to be a single positive integer!")
  }

  if (is.null(n_path)) {
    n_path <- floor(seq(nrow(x) / 2, nrow(x), nrow(x) / 20))
  }

  if (rsel) {
    hl_select <- dpca:::select_r(x,
      crit = rsel_crit, penalty_scales = penalty_scales,
      n_path = n_path, max_r = r
    )
    r <- hl_select$r
  }

  cx <- cov(t(x))
  mode(cx) <- "complex"

  edec <- .Call(
    "R_arnoldi_eigs",
    mat = cx,
    dim = nrow(cx),
    as.integer(r),
    .Machine$double.eps,
    1L,
    0L,
    0L,
    0L
  )

  factors <- t(Re(edec$vectors)) %*% x
  cc <- Re(edec$vectors) %*% factors + mx
  ic <- x - cc

  res <- list(
    xmean = mx,
    cov = Re(cx),
    eig = lapply(edec, Re),
    factors = factors,
    cc = cc,
    ic = ic
  )
  if (rsel) {
    res$HL_select <- hl_select
  }

  res
}
