test_that("Multivariate filter, computation", {

  nrx <- 100L
  ncx <- 1000L
  nrf <- 10L
  lags <- -2L:2L

  ## simple example (contamporaneous transformation)
  x <- matrix(rnorm(nrx * ncx), nrow = nrx)
  f <- matrix(0, nrow = nrf, ncol = nrx)
  f[1:nrf, 1:nrf] <- diag(nrf)

  res1 <- .Call("R_filter_process", f, x, 0L, nrf, nrx, nrx, ncx, 1L, 0L, 0L, 0L)
  res2 <- f %*% x
  expect_lt(sum((res1 - res2)^2), 1e-9)

  f_arr <- array(f, dim = c(dim(f), 1))
  res3 <- multivariate_filter(x = x, filter_coefficients = f_arr, lags = 1)
  expect_equal(res3, res2[,-1000])

})


test_that("multivariate filter, interface", {

  nrx <- 100L
  ncx <- 1000L
  nrf <- 10L
  lags <- -2L:2L

  ## simple example (contamporaneous transformation)
  x <- matrix(rnorm(nrx * ncx), nrow = nrx)
  f <- matrix(0, nrow = nrf, ncol = nrx)
  f[1:nrf, 1:nrf] <- diag(nrf)
  f_arr <- array(f, dim = c(dim(f), 1))
  res1 <- multivariate_filter(x, f_arr)
  res2 <- multivariate_filter(ts(t(x)), f_arr)
  expect_equal(res1, res2)

  expect_error(multivariate_filter("anc", f_arr))

  f_arr2 <- array(c(f,f), dim = c(dim(f), 2))
  f_arr3 <- array(c(f,f), dim = c(dim(f), 1) + c(0,1,0))

  expect_error(multivariate_filter(x, f_arr2, 1))
  expect_error(multivariate_filter(x, f_arr3, 1))

})
