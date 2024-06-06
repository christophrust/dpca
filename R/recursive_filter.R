#' Multivariate recursive filter.
#'
#' Apply a recursive filter to a multivariate time series.
#'
#' @param x A multivariate time series provided either as a matrix (rows correspond
#' to cross-sectional units and columnts to observations in the time domain) or a
#' multivariate object of class \code{\link[stats]{ts}} or \code{\link[zoo]{zoo}}.
#'
#' @param filter_coefficients An \eqn{n} by \eqn{m} by \eqn{p} array holding
#' filter coefficients, where \eqn{p} is the maximum lag order.
#'
#' @return A matrix of dimension \eqn{n} by \eqn{T - p} where
#' \eqn{n} and \eqn{T} are the input dimension, and
#'
#' \deqn{%
#'   y_t = f_{(,,1)} * y_{t-1} + ... +  f_{(,,p)} * y_{t-p} + x_t
#' }{%
#'   y[t] = f[,,1] * y[t-1] + ... +  f[,,p] * y[t-p] + x[t]
#' }
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
