library(dpca)


## random spectrum
test_that("Lagged covariance estimation", {

  ncy <- ncx <- 1000L
  nrx <- 100L
  nry <- 150L

  x <- matrix(rnorm(ncx* nrx), ncol = ncx)
  y <- matrix(rnorm(ncy* nry), ncol = ncy)

  system.time(res1 <- .Call("R_lagged_cov", x, y, 0L, nrx, ncx, nry, ncy))
  system.time(res2 <- tcrossprod(x, y)/ncx)

  expect_lt(sum((res1 - res2)^2), 1e-10)


  system.time(res1 <- .Call("R_lagged_cov", x, y, 1L, nrx, ncx, nry, ncy))
  system.time(res2 <- tcrossprod(x[,-1], y[,-ncy])/(ncx-1))

  expect_lt(sum((res1 - res2)^2), 1e-10)

  system.time(res1 <- .Call("R_lagged_cov", x, y, -1L, nrx, ncx, nry, ncy))
  system.time(res2 <- tcrossprod(x[,-ncx], y[,-1])/(ncx-1))

  expect_lt(sum((res1 - res2)^2), 1e-10)

  system.time(res1 <- .Call("R_lagged_cov", x, y, -2L, nrx, ncx, nry, ncy))
  system.time(res2 <- tcrossprod(x[,-ncx+c(1,0)], y[,-c(1:2)])/(ncx-2))

  expect_lt(sum((res1 - res2)^2), 1e-10)

})
