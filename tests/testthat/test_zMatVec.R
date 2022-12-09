library(dpca)


test_that("zMatVec", {

  m <- matrix(exp(runif(100) * 1i), ncol = 10)
  v <- exp(runif(10) * 1i)

  ## simple example
  r1 <- .Call("R_zMatVec",
              m,
              v,
              10L,
              2L,
              PACKAGE = "dpca")
  r2 <- as.vector(m %*% v)
  expect_equal(r1, r2)


  ## example where dim < ldm
  r1 <- .Call("R_zMatVec",
              m,
              v,
              8L,
              2L,
              PACKAGE = "dpca")
  r2 <- as.vector(m[1:8, 1:8] %*% v[1:8])
  expect_equal(r1, r2)

})
