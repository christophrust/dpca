#' Static Principal Component Analysis and Large Dynamic Factor Model Estimation
#'
#' @param x Input data supplied either as a matrix (rows correspond to cross-sectional
#' units and columnts to observations in the time domain) or a multivariate object of
#' class `ts` or `zoo`.
#'
#' @param r Number of static factors. If rsel = TRUE, this is
#' the maximum number of static factors.
#'#'
#' @param rsel Logical, if TRUE one of the Hallin & Liska criteria are used to
#' choose q from the data.
#'
#' @param rsel_crit Criterion to select the number of factors using the
#'   Hallin & Liska (2007, JASA) method. Either "IC1" or "IC2".
#'
#' @param n_path Integer vector specifying which (nested) subsets of the
#' cross section are used in the Hallin & Liska procedure.
#'
#' @param penalty_scales Tuning values for the penalty scaling parameter
#' c over which the q-path is optimized to stability.
#'
#' @return An object of class "spca".
#'
#' @importFrom stats is.ts
#' @export
spca <- function(
  x,
  r,
  rsel = FALSE,
  rsel_crit = c("IC1", "IC2", "IC3"),
  n_path = floor(seq(nrow(x) / 2, nrow(x), nrow(x) / 20)),
  penalty_scales = seq(0, 2, by = 0.01)
) {

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

  if (length(r) > 1 || floor(abs(r)) != r)
    stop("\"r\" has to be a single positive integer!")

  if (rsel) {
    hl_select <- select_r(x, crit = rsel_crit, penalty_scales = penalty_scales, n_path = n_path, max_r = r)
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
    xmean <- mx,
    cov = cx,
    eig = edec,
    factors = factors,
    cc = cc,
    ic = ic,
    HL_select = hl_select
  )

  class(res) <- "dfm"
  res
}
