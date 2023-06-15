#' Multivariate filter
#'
#' @param x a multivariate time series provided as a q by T matrix
#' @param filter an n by q by l array
#' @param ... further parameters passed to the underlying c routine,
#'       currently without functionality
#'
#' @return a matrix of dimension n by T where
#'   y[t] = f[,,1] * x[t] + f[,,2] * x[t-1] ...
#' @export
filter_process <- function(x, filter, ...) {

  res1 <- .Call("R_filter_process",
                as.numeric(filter),
                as.numeric(x),
                as.integer(seq_len(dim(filter)[3]) - 1),
                as.integer(nrow(filter)),
                as.integer(nrow(x)),
                as.integer(nrow(x)),
                as.integer(ncol(x)),
                1L, 0L, 0L, 0L)
}
