#' Multivariate filter
#'
#' Applies linear filter to a multivariate time series.
#'
#' @param x A multivariate time series provided either as a matrix (rows correspond
#' to cross-sectional units and columnts to observations in the time domain) or a
#' multivariate object of class `ts` or `zoo`.
#'
#' @param filter_coefficients An n by q by l array holding
#' filter coefficients.
#'
#' @param lags A vector of length dim(filter_coefficients)[3]
#' specifying the lag index.
#'
#' @return A matrix of dimension n by T - max(c(0, lags)) - max(c(0, -lags)) where
#' n and T are the input dimension, and
#'
#'   y[t] = f[,,1] * x[t - lags[1]] + f[,,2] * x[t-lags[2]] ...
#'
#' @importFrom stats is.ts
#' @export
multivariate_filter <- function(
  x,
  filter_coefficients,
  lags = seq_len(dim(filter_coefficients)[3])
){

  x <- if (is.ts(x) || "zoo" %in% class(x)) {
         t(x)
       } else if (is.matrix(x)) {
         x
       } else {
         stop("x must either a \"ts\" or \"zoo\" object or a matrix!")
       }

  if (dim(x)[1] != dim(filter_coefficients)[2])
    stop("Incompatible dimension of \"x\" and \"filter_coefficients\"!")

  if (length(lags) != dim(filter_coefficients)[3])
    stop("Length of \"lags\" must be equal to the third dimension of \"filter_coefficients\"!")

  l_idx <- max(c(0, lags)) + 1
  u_idx <- ncol(r) - max(c(0, -lags))

  r <- .Call(
    "R_filter_process",
    as.numeric(filter_coefficients),
    as.numeric(x),
    as.integer(lags),
    as.integer(nrow(filter_coefficients)),
    as.integer(nrow(x)),
    as.integer(nrow(x)),
    as.integer(ncol(x)),
    as.integer(length(lags)), 0L, 0L, 0L
  )

  r[ ,l_idx:u_idx]
}
