
test_that("select_r, computation", {

  ## settings
  nrx <- 100L
  ncx <- 400L
  q <- 2L
  bw <- as.integer(floor(ncx^(1/3)))
  freqs <- -10:10/10 * pi
  weights <- weights.Bartlett(-bw:bw/bw)
  ones <- rep(1, length(weights))

  ## simulate some data
  set.seed(123456)
  epsilon <- matrix(rnorm(ncx * q), nrow = q)

  b_filter <- vapply(1:10, function(l) {
    matrix(rnorm(q * nrx, sd = 1/l), q, nrx)
  }, matrix(0, q, nrx))

  chi <- .Call("R_filter_process", b_filter, epsilon, as.integer(1:10),
               nrx, q, q, ncx, 10L, 1L, 0L, 0L)
  x <- chi + rnorm(nrx * ncx, sd = 0.1 * sd(chi))
  x <- x - rowMeans(x)

  n_path <- floor(seq(nrow(x) / 2, nrow(x), nrow(x) / 20))
  penalty_scales <- seq(0, 2, by = 0.01)

  rr <- select_r(x, n_path = n_path, max_r = 26)
  expect_named(rr)

})
