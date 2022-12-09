#' Dynamic Principal Component Analysis and Large Dynamic Factor Models
#'
#' @param x Input data (n by T) matrix
#' @param q Number of dynamic factors. If qsel = TRUE, this is
#' the maximum number of dynamic factors.
#' @param freqs Numeric vector in [-pi, pi] giving the frequencies where the
#' spectral density is evaluated.
#' @param bandwidth Single integer, giving the width of the lag window estimator.
#' @param weights Kernel used for the lag window estimation of spectrum.
#' @param qsel Logical, if TRUE one of the Hallin & Liska criteria are used to
#' choose q from the data.
#' @param qsel_crit Criterion to select the number of factors using the
#'   Hallin & Liska (2007, JASA) method. Either "IC1" or "IC2".
#' @param n_path Integer vector specifying Which (nested) subsets of the
#' cross section are used in the Hallin & Liska procedure.
#'
#' @return An object of class "dpca" with different entries.
#'
#' @useDynLib dpca
#' @export
dpca <- function(x,
                 q,
                 freqs = -100:100/100 * pi,
                 bandwidth = floor(nrow(x)^(1/3)),
                 weights = c("Bartlett", "trunc", "Tukey", "Parzen", "Bohman", "Daniell", "ParzenCogburnDavis"),
                 qsel = FALSE,
                 qsel_crit = c("IC1", "IC2"),
                 n_path = floor(seq(nrow(x)/2, nrow(x), nrow(x)/20))) {

  ##browser()
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
               as.integer(q),
               FALSE,
               as.integer(n_path),
               PACKAGE = "dpca")

  class(res) <- "dpca"
  res
}
