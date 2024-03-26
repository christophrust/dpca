#' Multivariate filter
#'
#' Applies linear filter to a multivariate time series.
#'
#' @param x A multivariate time series provided as a q by T matrix.
#' @param filter_coefficients An n by q by l array holding
#' filter coefficients.
#' @param lags A vector of length dim(filter_coefficients)[3]
#' specifying the lag index.
#'
#'
#' @return A matrix of dimension n by T where
#'   y[t] = f[,,1] * x[t - lags[1]] + f[,,2] * x[t-lags[2]] ...
#'
#' @export
multivariate_filter <- function(
  x,
  filter_coefficients,
  lags = seq_len(dim(filter_coefficients)[3])
){

  if (dim(x)[1] != dim(filter_coefficients)[2])
    stop("Incompatible dimension of \"x\" and \"filter_coefficients\"!")

  if (length(lags) != dim(filter_coefficients)[3])
    stop("Length of \"lags\" must be equal to the third dimension of \"filter_coefficients\"!")

  res1 <- .Call("R_filter_process",
                as.numeric(filter_coefficients),
                as.numeric(x),
                as.integer(lags),
                as.integer(nrow(filter_coefficients)),
                as.integer(nrow(x)),
                as.integer(nrow(x)),
                as.integer(ncol(x)),
                as.integer(length(lags)), 0L, 0L, 0L)
}
