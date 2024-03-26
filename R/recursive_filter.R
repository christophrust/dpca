#' Multivariate recursive filter.
#'
#' Apply a recursive filter to a multivariate time series.
#'
#' @param x A multivariate time series provided as a q by T matrix.
#' @param filter_coefficients An n by q by p array holding
#' filter coefficients, where p is the maximum lag order.
#'
#' @export
recursive_filter <- function(x, filter_coefficients) {

  if (dim(x)[1] != dim(filter_coefficients)[2])
    stop("Incompatible dimension of \"x\" and \"filter_coefficients\"!")

  nlags <- if (is.matrix(filter_coefficients)) 1L else dim(filter_coefficients)[3]

  .Call("R_recursive_filter", x, as.numeric(filter_coefficients), nrow(x), nlags)
}
