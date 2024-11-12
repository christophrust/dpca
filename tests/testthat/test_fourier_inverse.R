library(dpca)
library(freqdom)

test_that("Fourier inverse (C-function)", {
  ## symmetric case
  ncx <- 1000L
  nrx <- 100L

  x <- matrix(rnorm(ncx * nrx), ncol = ncx)


  covs <- .Call("R_lagged_covs", x, x, as.integer(-5:5), nrx, ncx, nrx, ncx, rep(1, 11), 0L)

  spec <- .Call(
    "R_fourier_transform", as.vector(covs), dim(covs)[1], dim(covs)[2],
    -100:100 / 100 * pi, 201L, -5L:5L, 11L
  )

  eigs <- vapply(seq_len(dim(spec)[3]), function(i) {
    edec <- .Call("R_arnoldi_eigs",
      mat = spec[, , i], dim = dim(spec)[1],
      4L, .Machine$double.eps, 1L, 0L, 0L, 0L
    )
    edec$vectors
  }, matrix(0 + 1i, dim(spec)[1], 4))

  res1 <- suppressWarnings(.Call(
    "R_fourier_inverse", as.vector(eigs), dim(eigs)[1], dim(eigs)[2],
    as.integer(-5:5), 11L, -100:100 / 100 * pi, 201L
  ))
  res2 <- fourier.inverse(freqdom(eigs, freq = -100:100 / 100 * pi), as.integer(-5:5))$operators

  expect_lt(sum(abs(res2 - res1)^2), 1e-10)
})
