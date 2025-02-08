library(dpca)

test_that("Complex crossprod", {
  nrx <- 100
  ncx <- 10
  x <- matrix(runif(nrx * ncx) + runif(nrx * ncx) * 1i, nrow = nrx)

  res1 <- .Call("R_complex_crossprod", x, 0L)
  res2 <- Conj(t(x)) %*% x
  expect_equal(res1, res2)

  res1 <- .Call("R_complex_crossprod", x, 1L)
  res2 <- x %*% Conj(t(x))
  expect_equal(res1, res2)
})
