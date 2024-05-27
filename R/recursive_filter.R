#' Multivariate recursive filter.
#'
#' Apply a recursive filter to a multivariate time series.
#'
#' @param x A multivariate time series provided either as a matrix (rows correspond
#' to cross-sectional units and columnts to observations in the time domain) or a
#' multivariate object of class `ts` or `zoo`.
#'
#' @param filter_coefficients An n by q by p array holding
#' filter coefficients, where p is the maximum lag order.
#'
#' @return A matrix of dimension n by T - p where
#' n and T are the input dimension, and
#'
#'   y[t] = f[,,1] * y[t-1] + ... +  f[,,p] * y[t-p] + x[t]
#'
#' @importFrom stats is.ts
#' @export
recursive_filter <- function(x, filter_coefficients) {

  x <- if (is.ts(x) || "zoo" %in% class(x)) {
         t(x)
       } else if (is.matrix(x)) {
         x
       } else {
         stop("x must either a \"ts\" or \"zoo\" object or a matrix!")
       }

  if (dim(x)[1] != dim(filter_coefficients)[2])
    stop("Incompatible dimension of \"x\" and \"filter_coefficients\"!")

  nlags <- if (is.matrix(filter_coefficients)) 1L else dim(filter_coefficients)[3]

  .Call("R_recursive_filter", x, as.numeric(filter_coefficients), nrow(x), nlags)
}
