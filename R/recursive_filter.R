#' Multivariate recursive filter
#'
#' @param x A multiple time series supplied as a n by T matrix
#' @param f Filter coefficient matrices supplied as n by n by p array,
#'   where p is the maximum lag order.
#' @export
recursive_filter <- function(x, f) {

  nlags <- if (is.matrix(f)) 1L else dim(f)[3]

  .Call("R_recursive_filter", x, as.numeric(f), nrow(x), nlags)
}
