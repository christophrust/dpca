library(dpca)
##library(freqdom)

test_that("Fourier transform (C-function)", {


  ncy <- ncx <- 1000L
  nrx <- 100L
  nry <- 150L

  x <- matrix(rnorm(ncx* nrx), ncol = ncx)
  y <- matrix(rnorm(ncy* nry), ncol = ncy)

  covs <- .Call("R_lagged_covs", x, y, as.integer(-5:5), nrx, ncx, nry, ncy)

  res1 <- .Call("R_fourier_transform", as.vector(covs), dim(covs)[1], dim(covs)[2],
                -100:100/100 * pi, 201L, -5L:5L, 11L)

})
