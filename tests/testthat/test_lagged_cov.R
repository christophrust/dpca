library(dpca)


## random spectrum
test_that("Lagged covariance estimation", {

  ncy <- ncx <- 1000L
  nrx <- 10L
  nry <- 15L

  x <- matrix(rnorm(ncx* nrx), ncol = ncx)
  y <- matrix(rnorm(ncy* nry), ncol = ncy)
  xc <- x - rowMeans(x)
  yc <- y - rowMeans(y)

  ## contemporaneous, no center
  res1 <- .Call("R_lagged_cov", x, y, 0L, nrx, ncx, nry, ncy, 1, 0L)
  res2 <- tcrossprod(x, y)/ncx
  expect_lt(sum((res1 - res2)^2), 1e-10)

  ## contemporaneous, centered
  res1 <- .Call("R_lagged_cov", x, y, 0L, nrx, ncx, nry, ncy, 1, 1L)
  res2 <- tcrossprod(xc, yc)/ncx
  expect_lt(sum((res1 - res2)^2), 1e-10)


  ## lag1, no center
  res1 <- .Call("R_lagged_cov", x, y, 1L, nrx, ncx, nry, ncy, 1, 0L)
  res2 <- tcrossprod(x[,-1], y[,-ncy])/ncx
  expect_lt(sum((res1 - res2)^2), 1e-10)


  res1 <- .Call("R_lagged_cov", x, y, -1L, nrx, ncx, nry, ncy, 1, 0L)
  res2 <- tcrossprod(x[,-ncx], y[,-1])/ncx
  expect_lt(sum((res1 - res2)^2), 1e-10)


  res1 <- .Call("R_lagged_cov", x, y, -2L, nrx, ncx, nry, ncy, 1, 0L)
  res2 <- tcrossprod(x[,-ncx+c(1,0)], y[,-c(1:2)])/ncx
  expect_lt(sum((res1 - res2)^2), 1e-10)


  res1 <- .Call("R_lagged_cov", x, y, -2L, nrx, ncx, nry, ncy, 1, 1L)
  res2 <- .Call("R_lagged_cov", xc, yc, -2L, nrx, ncx, nry, ncy, 1, 0L)
  res3 <- tcrossprod(xc[,-ncx+c(1,0)], yc[,-c(1:2)])/ncx
  expect_lt(sum((res1 - res2)^2), 1e-10)
  expect_lt(sum((res1 - res3)^2), 1e-10)
  expect_lt(sum((res2 - res3)^2), 1e-10)

})


test_that("R_lagged_cov input checking (dimension of arrays)", {


  ncy <- ncx <- 100L
  nrx <- 10L
  nry <- 15L

  x <- matrix(rnorm(ncx* nrx), ncol = ncx)
  y <- matrix(rnorm(ncy* nry), ncol = ncy)

  expect_error(.Call("R_lagged_cov", x, y, 0L, nrx +1L, ncx, nry, ncy, 1, 0L))
  expect_error(.Call("R_lagged_cov", x, y, 0L, nrx, ncx + 1L, nry, ncy, 1, 0L))
  expect_error(.Call("R_lagged_cov", x, y, 0L, nrx, ncx, nry + 2L, ncy, 1, 0L))
  expect_error(.Call("R_lagged_cov", x, y, 0L, nrx, ncx, nry, ncy + 3L, 1, 0L))

})


test_that("R_lagged_cov input checking (equal number of columns)", {


  ncy <- ncx <- 100L
  nrx <- 10L
  nry <- 15L

  x <- matrix(rnorm(ncx* nrx), ncol = ncx)
  y <- matrix(rnorm(ncy* nry + nry), ncol = ncy + 1)

  expect_error(.Call("R_lagged_cov", x, y, 0L, nrx, ncx, nry, ncy + 1L, 1, 0L))


})


test_that("R_lagged_cov input checking (lag too high)", {


  ncy <- ncx <- 100L
  nrx <- 10L
  nry <- 15L

  x <- matrix(rnorm(ncx* nrx), ncol = ncx)
  y <- matrix(rnorm(ncy* nry), ncol = ncy)

  expect_error(.Call("R_lagged_cov", x, y, 100L, nrx, ncx, nry, ncy, 0L))

})


## random spectrum
test_that("Multiple lagged covariance estimation", {

  ncy <- ncx <- 1000L
  nrx <- 100L
  nry <- 150L

  x <- matrix(rnorm(ncx* nrx), ncol = ncx)
  y <- matrix(rnorm(ncy* nry), ncol = ncy)
  xc <- x - rowMeans(x)
  yc <- y - rowMeans(y)

  weights <- rep(1, 11)
  system.time(res1 <- .Call("R_lagged_covs", x, y, as.integer(-5:5), nrx, ncx, nry, ncy, weights, 0L))
  system.time(res2 <-
                vapply(-5:5,
                       function(l) {
                         if (l<0)
                           return(tcrossprod(x[,-ncx+(-l-1):0], y[,-c(1:abs(l))])/ncx)
                         if (l == 0)
                           return(tcrossprod(x, y)/ncx)

                         tcrossprod(x[,-(1:l)], y[,-ncy + (l-1):0])/ncx

                       }, matrix(0, nrx, nry))
              )

  expect_lt(sum((res1 - res2)^2), 1e-10)

  system.time(res1 <- .Call("R_lagged_covs", xc, yc, as.integer(-5:5), nrx, ncx, nry, ncy, weights, 0L))
  system.time(res2 <- .Call("R_lagged_covs", x, y, as.integer(-5:5), nrx, ncx, nry, ncy, weights, 1L))
  expect_lt(sum((res1 - res2)^2), 1e-10)
})
