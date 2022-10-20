library(dpca)

test_that("Filter process (C-function)", {

  nrx <- 100L
  ncx <- 1000L
  nrf <- 10L
  x <- matrix(rnorm(nrx * ncx), nrow = nrx)
  f <- matrix(0, nrow = nrf, ncol = nrx)
  f[1:nrf, 1:nrf] <- diag(nrf)

  res1 <- .Call("R_filter_process", f, x, 0L, nrf, nrx, nrx, ncx, 1L, 0L, 0L, 0L)
  res2 <- f %*% x

  expect_lt(sum((res1 - res2)^2), 1e-9)


  res1 <- .Call("R_filter_process", f, x, 1L, nrf, nrx, nrx, ncx, 1L, 0L, 0L, 0L)
  res2 <- cbind(NA, f %*% x[,-ncol(x)])

  expect_lt(sum((res1 - res2)^2, na.rm = TRUE), 1e-9)






  nrx <- 100L
  ncx <- 1000L
  nrf <- 5L
  lags <- -2L:2L
  x <- matrix(rnorm(nrx * ncx), nrow = nrx)
  f <- vapply(seq_along(lags), function(i) {
    m <- matrix(0, nrow = nrf, ncol = nrx)
    m[i,i] <- 1
    m
  }, matrix(0,ncol = nrx, nrow = nrf))

  res1 <- .Call("R_filter_process", f, x, lags, nrf, nrx, nrx, ncx, length(lags), 0L, 0L, 0L)
  tmp <- vapply(seq_along(lags), function(l) {
    if (lags[l] > 0)
    {
      tmp <- f[,,l] %*% x[, -(ncol(x) - 1:lags[l] + 1)]
      return(cbind(matrix(NA_real_, nrow = nrow(tmp), ncol = lags[l]), tmp))
    }
    else if (lags[l] < 0){
      tmp <- f[,,l] %*% x[, -(1:(-lags[l]))]
      return(cbind(tmp, matrix(NA_real_, nrow = nrow(tmp), ncol = -lags[l])))
    } else
      return(f[,,l] %*% x)
  }, matrix(0, nrow = nrf, ncol = ncx))
  res2 <- apply(tmp, 1:2, sum)


  expect_lt(sum((res1 - res2)^2, na.rm = TRUE), 1e-9)

  if (FALSE){
    rowSums((res2 - res1)^2 > 1e-9, na.rm = TRUE)
    cbind(res1[4,1:10], res2[4,1:10], x[4,1:10])
    cbind(res1[2,1:10], res2[2,1:10], x[2,1:10])
  }

  nrx <- 100L
  ncx <- 1000L
  nrf <- 5L
  lags <- -2L:2L
  x <- matrix(rnorm(nrx * ncx), nrow = nrx)
  f <- vapply(lags, function(i) matrix(rnorm(nrx * nrf), nrow = nrf), matrix(0,ncol = nrx, nrow = nrf))

  res1 <- .Call("R_filter_process", f, x, lags, nrf, nrx, nrx, ncx, length(lags), 0L, 0L, 0L)
  tmp <- vapply(seq_along(lags), function(l) {
    if (lags[l] > 0)
    {
      tmp <- f[,,l] %*% x[, -(ncol(x) - 1:lags[l] + 1)]
      return(cbind(matrix(NA_real_, nrow = nrow(tmp), ncol = lags[l]), tmp))
    }
    else if (lags[l] < 0){
      tmp <- f[,,l] %*% x[, -(1:(-lags[l]))]
      return(cbind(tmp, matrix(NA_real_, nrow = nrow(tmp), ncol = -lags[l])))
    } else
      return(f[,,l] %*% x)
  }, matrix(0, nrow = nrf, ncol = ncx))
  res2 <- apply(tmp, 1:2, sum)

  ft <- vapply(seq_along(lags), function(i) {
    t(f[,,i])
  }, t(f[,,1]))
  res3 <- .Call("R_filter_process", ft, x, lags, nrx, nrf, nrx, ncx, length(lags), 0L, 1L, 0L)

  expect_lt(sum((res1 - res2)^2, na.rm = TRUE), 1e-9)
  expect_lt(sum((res1 - res3)^2, na.rm = TRUE), 1e-9)


  res1 <- .Call("R_filter_process", f, x, lags, nrf, nrx, nrx, ncx, length(lags), 1L, 0L, 0L)
  res2[is.na(res2)] <- rowMeans(res2, na.rm = TRUE)
  expect_lt(sum((res1 - res2)^2), 1e-9)
})
