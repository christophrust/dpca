test_that("Test spca, computation.", {
  ## settings
  nrx <- 100L
  ncx <- 1000L
  q <- 4L
  bw <- as.integer(floor(ncx^(1 / 3)))
  freqs <- -10:10 / 10 * pi
  weights <- weights_bartlett(-bw:bw / bw)
  ones <- rep(1, length(weights))

  ## simulate some data
  set.seed(123456)
  epsilon <- matrix(rnorm(ncx * q), nrow = q)

  b_filter <- vapply(1:10, function(l) {
    matrix(rnorm(q * nrx, sd = 1 / l), q, nrx)
  }, matrix(0, q, nrx))

  chi <- .Call(
    "R_filter_process", b_filter, epsilon, as.integer(1:10),
    nrx, q, q, ncx, 10L, 1L, 0L, 0L
  )
  x <- chi + rnorm(nrx * ncx, sd = 0.1 * sd(chi))
  x <- x - rowMeans(x)


  res_spca <- spca(x, r = 3)

  cmp_cov <- cov(t(x))
  expect_equal(res_spca$cov, cmp_cov)

  cmp_eig <- eigen(cmp_cov)
  expect_equal(cmp_eig$values[1:3], res_spca$eig$values)
  evec_cmp <- cmp_eig$vectors[, 1:3] / res_spca$eig$vectors
  expect_equal(unique(round(abs(colMeans(evec_cmp)), digits = 9)), 1)

  cmp_cc <- cmp_eig$vectors[, 1:3] %*% t(cmp_eig$vectors[, 1:3]) %*% x
  expect_equal(cmp_cc, res_spca$cc)

  cmp_ic <- x - cmp_cc
  expect_equal(cmp_ic, res_spca$ic)
})


test_that("Test spca, interface.", {
  ## settings
  nrx <- 100L
  ncx <- 400L
  q <- 2L
  bw <- as.integer(floor(ncx^(1 / 3)))
  freqs <- -10:10 / 10 * pi
  weights <- weights_bartlett(-bw:bw / bw)
  ones <- rep(1, length(weights))

  ## simulate some data
  set.seed(123456)
  epsilon <- matrix(rnorm(ncx * q), nrow = q)

  b_filter <- vapply(1:10, function(l) {
    matrix(rnorm(q * nrx, sd = 1 / l), q, nrx)
  }, matrix(0, q, nrx))

  chi <- .Call(
    "R_filter_process", b_filter, epsilon, as.integer(1:10),
    nrx, q, q, ncx, 10L, 1L, 0L, 0L
  )
  x <- chi + rnorm(nrx * ncx, sd = 0.1 * sd(chi))
  x <- x - rowMeans(x)

  ## wrong input for x
  expect_error(spca("abc"))

  ## wrong r supplied
  expect_error(spca(x, c(1:3)))
  expect_error(spca(x, 1.3))

  ## rselection
  rr <- spca(x, max_r = 10)
  expect_named(rr$HL_select)
})
