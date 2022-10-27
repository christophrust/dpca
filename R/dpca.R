#' Dynamic Principal Component Analysis and Large Dynamic Factor Models
#'
#' @param x Input data (n by T) matrix
#' @param q Number of dynamic eigenvalues
#' @param freqs Numeric vector in [-pi, pi] giving the frequencies where the
#' spectral density is evaluated.
#' @param bandwidth Single integer, giving the width of the lag window estimator.
#'
#' @return An object of class "dpca" with different entries.
#'
#' @useDynLib dpca
#' @export
dpca <- function(x, q, freqs = -100:100/100 * pi, bandwidth,
                 weights = c("Bartlett", "trunc", "Tukey", "Parzen", "Bohman", "Daniell", "ParzenCogburnDavis"), ...) {

  if (length(weights) > 1)
    weights = "Bartlett"

  if (!is.matrix(x))
    stop("x must be a n by T matrix!")

  if (missing(q)) {
    warning("No number of dynamic principal components supplied. Using q = 1...")
    q <- 1L
  }

  if (length(q) > 1 || floor(abs(q)) != q)
    stop("\"q\" has to be a single positive integer!")

  if (length(bandwidth) > 1 || floor(abs(bandwidth)) != bandwidth)
    stop("\"bandwidth\" has to be a single positive integer!")

  if (!is.numeric(freqs) || any(abs(freqs) > pi))
    stop("\"freqs\" must be a numeric vector with values in [-pi, pi]!")


  wghts <- get(paste0("weights.", weights))(-bandwidth:bandwidth/bandwidth)

  ##browser()
  res <- .Call("R_dpca", x,
               as.integer(q),
               as.numeric(freqs),
               as.integer(bandwidth),
               .Machine$double.eps,
               wghts,
               PACKAGE = "dpca")

  class(res) <- "dpca"
  res
}
