#' Multivariate filter
#'
#' Applies linear filter to a multivariate time series.
#'
#' @param x A multivariate time series provided either as a matrix (rows correspond
#' to cross-sectional units and columnts to observations in the time domain) or a
#' multivariate object of class \code{\link[stats]{ts}} or \code{\link[zoo]{zoo}}.
#'
#' @param filter_coefficients An \eqn{n} by \eqn{m} by \eqn{l} array holding
#' filter coefficients.
#'
#' @param lags A vector of length \eqn{l} specifying the lag corresponding to the
#' \eqn{l} subviews of filter_coefficients.
#'
#' @return A matrix of dimension \eqn{n} by
#' \eqn{T - \operatorname{max}(0, \texttt{lags}_1,\dots,\texttt{lags}_l) - %
#' \operatorname{max}(0, - \texttt{lags}_1,\dots,-\texttt{lags}_l)} where
#' \eqn{T} is the number of time observations of the input, and
#'
#' \deqn{%
#'   y_t = f_{(,,1)} x_{t - \texttt{lags}_1} + f_{(,,2)} x_{t - \texttt{lags}_1} + \dots
#' }{%
#'   y[t] = f[,,1] * x[t - lags[1]] + f[,,2] * x[t-lags[2]] + ...
#' }
#'
#' @importFrom stats is.ts
#' @export
multivariate_filter <- function(
    x,
    filter_coefficients,
    lags = seq_len(dim(filter_coefficients)[3])) {
  x <- if (is.ts(x) || "zoo" %in% class(x)) {
    t(x)
  } else if (is.matrix(x)) {
    x
  } else {
    stop("x must either a \"ts\" or \"zoo\" object or a matrix!")
  }

  x <- if (is.ts(x) || "zoo" %in% class(x)) {
         t(x)
       } else if (is.matrix(x)) {
         x
       } else {
         stop("x must either a \"ts\" or \"zoo\" object or a matrix!")
       }

  if (dim(x)[1] != dim(filter_coefficients)[2])
    stop("Incompatible dimension of \"x\" and \"filter_coefficients\"!")
  }

  if (length(lags) != dim(filter_coefficients)[3]) {
    stop("Length of \"lags\" must be equal to the third dimension of \"filter_coefficients\"!")
  }

  l_idx <- max(c(0, lags)) + 1
  u_idx <- ncol(x) - max(c(0, -lags))

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

  r[, l_idx:u_idx]
}
