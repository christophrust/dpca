library(dpca)
library(freqdom)


test_that("Test dpca, stepwise", {
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


  ## run dpca for both implementation (freqdom and our)
  system.time(res_dpca <- .Call("R_dpca", x, 4L, freqs, bw, .Machine$double.eps, weights, 4L, 0L, 0L, 0L, 0.1, 0.1))
  system.time(res_dpca1 <- dpca::dpca(x = x, q = 4, freqs = freqs, bandwidth = bw, weights = "bartlett"))
  system.time(res_freqdom <- freqdom::dpca(t(x), bw, freqs, 4L))
  res_dpca$xmean <- rowMeans(x)
  class(res_dpca1) <- NULL
  expect_equal(res_dpca, res_dpca1)


  ## 1. compute covariance structure without weights
  covs1 <- .Call("R_lagged_covs", x, x, -bw:bw, nrx, ncx, nrx, ncx, ones, 1L)
  covs2 <- cov.structure(t(x), t(x), -bw:bw)$operators
  expect_lt(sum(vapply(seq_len(length(weights)), function(i) sum((covs1[, , i] - covs2[, , i])^2), 0)), 1e-12)


  ## 2. compute covariance structure with weights
  covs1 <- .Call("R_lagged_covs", x, x, -bw:bw, nrx, ncx, nrx, ncx, weights, 1L)
  tmp <- cov.structure(t(x), t(x), -bw:bw)$operators
  covs2 <- vapply(
    1:dim(tmp)[3], function(i) tmp[, , i] * weights[i],
    matrix(0, dim(tmp)[1], dim(tmp)[1])
  )
  expect_lt(sum((covs1 - covs2)^2), 1e-12)


  ## 3. Compute spectrum
  Ch <- timedom(covs2, lags = -bw:bw)
  dft1 <- fourier.transform(Ch, -10:10 / 10 * pi)$operators
  dft2 <- .Call(
    "R_fourier_transform", as.vector(covs1), dim(covs1)[1], dim(covs1)[2],
    -10:10 / 10 * pi, 21L, -bw:bw, length(weights)
  )
  expect_lt(sum((Re(dft1) - Re(dft2))^2), 1e-12)
  expect_lt(sum((Im(dft1) - Im(dft2))^2), 1e-12)

  expect_lt(sum((Im(dft1) - Im(dft2))^2), 1e-12)
  expect_lt(sum((Im(res_dpca$spectrum) - Im(dft2))^2), 1e-12)
  expect_lt(sum((Re(res_dpca$spectrum) - Re(dft2))^2), 1e-12)
  expect_lt(sum((Re(res_freqdom$spec.density$operators) - Re(dft2))^2), 1e-12)
  expect_lt(sum((Im(res_freqdom$spec.density$operators) - Im(dft2))^2), 1e-12)


  ## 4. Eigendecomposition of spectrum (uses different numerical algorithms)
  spec1 <- res_dpca$spectrum
  spec2 <- res_freqdom$spec.density$operators
  expect_lt(abs(sum((spec1 - spec2)^2)), 1e-10)

  tmp1 <- lapply(
    seq_len(dim(spec1)[3]),
    function(i) {
      .Call("R_arnoldi_eigs",
        mat = spec1[, , i], dim = nrx, 4L,
        .Machine$double.eps, 1L, 0L, 1L, 1L
      )
    }
  )
  tmp2 <- freqdom:::freqdom.eigen(res_freqdom$spec.density)

  ## eigenvalues
  res1 <- vapply(tmp1, function(x) Re(x$values), numeric(4))
  res2 <- tmp2$values[1:4, ]
  expect_lt(sum((res1 - res2)^2), 1e-12)

  ## eigenvectors
  ev1 <- vapply(tmp1, function(x) x$vectors, matrix(1 + 0i, dim(tmp1[[1]]$vectors)[1], dim(tmp1[[1]]$vectors)[2]))
  ev2 <- if (dim(tmp1[[1]]$vectors)[1] == q) {
    tmp2$vectors[1:4, , ]
  } else {
    vapply(
      seq_len(dim(tmp2$vectors)[3]),
      function(i) t(tmp2$vectors[1:4, , i]),
      matrix(1 + 0i, nrx, q)
    )
  }
  diff <- vapply(seq_len(dim(ev1)[3]), function(i) sum((abs(ev1[, , i] / ev2[, , i]) - 1)^2), 0 + 0i)
  expect_lt(sum(abs(diff)), 1e-12)



  ## resulting filters
  f1 <- fourier.inverse(freqdom(tmp2$vectors[1:4, , ], freq = -10:10 / 10 * pi), -bw:bw)$operators
  f2 <- .Call(
    "R_fourier_inverse", tmp2$vectors[1:4, , ], 4L, 100L,
    as.integer(-bw:bw), length(-bw:bw), -10:10 / 10 * pi, 21L
  )
  expect_lt(sum((f1 - f2)^2), 1e-12)


  ## approach 2: compute p~ * p and apply fourier transformation thereafter.
  k_theta1 <- vapply(seq_len(dim(ev1)[3]), function(i) t(Conj(ev1[, , i])) %*% ev1[, , i], matrix(0 + 0i, nrx, nrx))
  k_theta2 <- vapply(seq_len(dim(ev2)[3]), function(i) t(Conj(ev2[, , i])) %*% ev2[, , i], matrix(0 + 0i, nrx, nrx))
  expect_lt(sum(abs(k_theta1 - k_theta2)^2), 1e-10)

  ff1 <- .Call(
    "R_fourier_inverse", k_theta1, 100L, 100L,
    as.integer(-bw:bw), length(-bw:bw), -10:10 / 10 * pi, 21L
  )

  expect_equal(ff1, res_dpca$filter$dcc)

  ## SEXP R_filter_process(SEXP r_f, SEXP r_x, SEXP r_lags,
  ##                     SEXP r_nrf, SEXP r_ncf, SEXP r_nrx, SEXP r_ncx,
  ##                     SEXP r_nlags, SEXP r_inx, SEXP r_transf, SEXP r_rev)
  dcc1 <- .Call(
    "R_filter_process", ff1, x, as.integer(-bw:bw), nrx, nrx, nrx, ncx,
    length(-bw:bw), 1L, 0L, 0L
  )
  dcc2 <- res_freqdom$Xhat
  dcc3 <- res_dpca$dcc
  expect_equal(dcc1, dcc3)


  ## check superiority of dpca compared to freqdom
  r2_dpca <- summary(lm(as.vector(chi) ~ as.vector(dcc3)))$r.squared
  r2_freqdom <- summary(lm(as.vector(chi) ~ as.vector(t(dcc2))))$r.squared
  expect_lt(r2_freqdom, r2_dpca)
})


test_that("dpca: HL switch", {
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


  res_dpca1 <- dpca::dpca(
    x = x, q_max = 10, freqs = freqs, bandwidth = bw, weights = "bartlett",
  )


  res_dpca2 <- dpca::dpca(
    x = x, q = res_dpca1$HL_select$q, freqs = freqs, bandwidth = bw, weights = "bartlett"
  )

  expect_lt(sum(abs(res_dpca1$spectrum[, , 1] %*% Conj(res_dpca1$eig$vectors[1, , 1]) -
    Conj(res_dpca1$eig$vectors[1, , 1]) * res_dpca1$eig$values[1, 1])), 1e-10)

  expect_equal(res_dpca1$eig, res_dpca2$eig)
  expect_equal(res_dpca1$dcc, res_dpca2$dcc)
  expect_equal(res_dpca1$input, res_dpca2$input)
})


test_that("dpca: inferface", {
  ## settings
  nrx <- 100L
  ncx <- 1000L
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

  res_dpca1 <- dpca::dpca(
    x = x, q = 10, freqs = freqs, bandwidth = bw
  )
  res_dpca2 <- dpca::dpca(
    x = ts(t(x)), q = 10, freqs = freqs, bandwidth = bw
  )

  expect_equal(res_dpca1$eig, res_dpca2$eig)

  expect_error(
    dpca::dpca(
      x = "test", q = 10, freqs = freqs, bandwidth = bw
    )
  )

  expect_error(
    dpca::dpca(
      x = x, q = 1.23, freqs = freqs, bandwidth = bw
    )
  )

  expect_warning(
    dpca::dpca(
      x = x, q = 4, freqs = c(-pi, -2, seq(0, pi, length.out = 5)), bandwidth = bw
    )
  )
})
