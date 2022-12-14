library(dpca)


context("Hallin & Liska (2007) selection criterium")

test_that("Get q path (c function)", {

  max_q <- 12L
  ln <- 10L
  unpenalized_crit <- matrix(runif(max_q * ln, 100, 140), ncol = ln)
  penalties <- runif(ln, min = 10, max = 12)
  penalty_scale <- 0.001
  r1 <- .Call("R_hl_q_path", unpenalized_crit, max_q, penalty_scale, penalties)
  r2 <- apply(unpenalized_crit + 1:max_q %o% penalties * penalty_scale, 2, which.min)

  expect_equal(r1, r2)
})
