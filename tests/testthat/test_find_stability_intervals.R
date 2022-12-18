library(dpca)


context("Hallin & Liska (2007) selection criterium")


test_that("Find stability interval",{

  set.seed(as.integer(Sys.Date()))
  intvl_idx <- c(1,sort(sample(2:100, 3)))
  sample_vars <- rep(0, 100)
  repl_idx <- c((intvl_idx[2]+1):(intvl_idx[3]-1) , (intvl_idx[4]+1):100)
  sample_vars[repl_idx] <- abs(rnorm(length(repl_idx)))
  ## sample_vars <- cumsum(rnorm(100, mean = 1, sd = ))
  ## sample_vars <- sample_vars - seq(0, 100, len = 100)/100 * sample_vars[100]
  ## sample_vars <- dplyr::if_else(sample_vars<0, 0, sample_vars)
  ## zero_idx <- which(sample_vars == 0)
  ## zero_idx[-length(zero_idx)] == zero_idx[-1] -1
  ## plot(sample_vars, type = "l")

  r1 <- .Call("R_find_stability_intervals", sample_vars)
  expect_equal(r1, intvl_idx -1)



  intvl_idx <- c(1,sort(sample(2:100, 3)))
  sample_vars <- rep(1, 100)
  repl_idx <- c((intvl_idx[2]+1):(intvl_idx[3]-1) , (intvl_idx[4]+1):100)
  sample_vars[repl_idx] <- 1 + abs(rnorm(length(repl_idx)))
  ## sample_vars <- cumsum(rnorm(100, mean = 1, sd = ))
  ## sample_vars <- sample_vars - seq(0, 100, len = 100)/100 * sample_vars[100]
  ## sample_vars <- dplyr::if_else(sample_vars<0, 0, sample_vars)
  ## zero_idx <- which(sample_vars == 0)
  ## zero_idx[-length(zero_idx)] == zero_idx[-1] -1
  ##plot(sample_vars, type = "l")
  expect_warning(r1 <- .Call("R_find_stability_intervals", sample_vars))
  expect_equal(r1, intvl_idx -1)


  sample_vars <- (1:100)^0.2
  sample_vars[c(0:4, 25:40)] <- 1
  sample_vars[10:14] <- 14^0.2
  expect_warning(r1 <- .Call("R_find_stability_intervals", sample_vars))
  expect_equal(r1, c(0,3,24,39))



  sample_vars <- rep(0, 100)
  sample_vars[16:25] <- abs(rnorm(10))
  r1 <- .Call("R_find_stability_intervals", sample_vars)
  expect_equal(r1, c(0,14,25,99))


  sample_vars <- cumsum(c(0,rnorm(99)))
  expect_warning(r1 <- .Call("R_find_stability_intervals", sample_vars))
  expect_equal(r1[1], which.min(sample_vars) - 1)

})
