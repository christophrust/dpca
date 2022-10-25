library(dpca)
##library(freqdom)

test_that("Fourier inverse (C-function)", {

  ## symmetric case
  ncx <- 1000L
  nrx <- 100L

  x <- matrix(rnorm(ncx* nrx), ncol = ncx)


  covs <- .Call("R_lagged_covs", x, x, as.integer(-5:5), nrx, ncx, nrx, ncy, rep(1, 11), 0L)

  spec <- .Call("R_fourier_transform", as.vector(covs), dim(covs)[1], dim(covs)[2],
                -100:100/100 * pi, 201L, -5L:5L, 11L)

  eigs <- vapply(seq_len(dim(spec)[3]), function(i) {
    edec <- .Call("R_arnoldi_eigs", mat = spec[,,i], dim = dim(spec)[1],
                  4L, .Machine$double.eps, 1L, 0L)
    edec$vectors
  }, matrix(0 + 1i, dim(spec)[1], 4))

  res1 <- .Call("R_fourier_inverse", as.vector(eigs), dim(eigs)[1], dim(eigs)[2],
                   as.integer(-5:5), 11L, -100:100/100 * pi, 201L, 0L)

  res2 <- .Call("R_fourier_inverse1", as.vector(eigs), dim(eigs)[1], dim(eigs)[2],
                   as.integer(-5:5), 11L, -100:100/100 * pi, 201L)

  sum((res2 - res1)^2)



  ## assymmetric case
  ncy <- ncx <- 1000L
  nrx <- 100L
  nry <- 150L

  x <- matrix(rnorm(ncx* nrx), ncol = ncx)
  y <- matrix(rnorm(ncy* nry), ncol = ncy)


  covs <- .Call("R_lagged_covs", x, y, as.integer(-5:5), nrx, ncx, nry, ncy, rep(1, 11), 0L)

  spec <- .Call("R_fourier_transform", as.vector(covs), dim(covs)[1], dim(covs)[2],
                -100:100/100 * pi, 201L, -5L:5L, 11L)

  eigs <- vapply(seq_len(dim(spec)[3]), function(i) {
    edec <- .Call("R_arnoldi_eigs", mat = spec[,,i], dim = dim(spec)[1],
                  4L, .Machine$double.eps, 1L, 0L)
    edec$vectors
  }, matrix(0 + 1i, dim(spec)[1], 4))

  fu1 <- function(){
    res1 <- .Call("R_fourier_inverse", as.vector(eigs), dim(eigs)[1], dim(eigs)[2],
                  as.integer(-5:5), 11L, -100:100/100 * pi, 201L, 0L)
  }
  fu2 <- function() {
    res2 <- .Call("R_fourier_inverse1", as.vector(eigs), dim(eigs)[1], dim(eigs)[2],
                  as.integer(-5:5), 11L, -100:100/100 * pi, 201L)
  }

  ## minibenchmark
  ##bm <- microbenchmark::microbenchmark(fu1(), fu2(), times = 100)

})
