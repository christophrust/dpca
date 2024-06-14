library(dpca)

test_that("Recursive filter (C-function)", {

  x <- matrix(rnorm(300), nrow = 3)
  f <- array(matrix(rnorm(18, sd = 0.1), ncol = 3), dim = c(3, 3, 2))
  res1 <- .Call("R_recursive_filter", x, as.vector(f), 3L, 2L)
  res2 <- x
  for (i in 2:ncol(x)) {
    res2[, i] <- if (i > 2) {
      res2[, i] + f[, , 1] %*% res2[, i - 1] + f[, , 2] %*% res2[, i - 2]
    } else {
      res2[, i] + f[, , 1] %*% res2[, i - 1]
    }
  }

  expect_equal(res1, res2)

  res3 <- recursive_filter(x, f)
  expect_equal(res3, res2)

  res4 <- recursive_filter(ts(t(x)), f)
  expect_equal(res4, res2)

  expect_error(recursive_filter("test", f))

  expect_error(recursive_filter(matrix(1), f))

})
