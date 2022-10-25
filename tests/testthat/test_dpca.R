library(dpca)
library(freqdom)

test_that("Test dpca (C-function)", {

  nrx <- 100L
  ncx <- 1000L
  q <- 4L

  epsilon <- matrix(rnorm(ncx * q), nrow = q)

  b_filter <- vapply(1:10, function(l) {
    matrix(rnorm(q * nrx, sd = 1/l), q, nrx)
  }, matrix(0, q, nrx))

  chi <- .Call("R_filter_process", b_filter, epsilon, as.integer(1:10),
               nrx, q, q, ncx, 10L, 1L, 0L, 0L)

  x <- chi + rnorm(nrx * ncx, sd = 0.1 * sd(chi))

  weights <- weights.Bartlett(-15:15/15)
  system.time(res_dpca <- .Call("R_dpca", x, 4L, -10:10/10 * pi, 15L, .Machine$double.eps, weights))

  system.time(res_freqdom <- freqdom::dpca(X = t(x), q = 15, freq = -10:10/10 * pi, Ndpc = 4))
  str(res_freqdom)
  str(res_dpca)

  plot(x = as.vector(x), y = as.vector(t(res_freqdom$Xhat)))
  points(x = as.vector(x), y = as.vector(res_dpca$dcc), col = "red")
  res1 <- Re(res_dpca$spectrum)
  res2 <- Re(res_freqdom$spec.density$operators)
  apply((res1 - res2)^2 > 1e-9, 3, sum)
  rowSums((res1[,,1] - res2[,,1])^2 > 1e-9)
  cbind(res1[1,,1], res2[1,,1])
  plot(x = as.vector(Re(res_dpca$spectrum)), y = as.vector(Re(res_freqdom$spec.density$operators)))
  lm(as.vector(Re(res_dpca$spectrum)) ~as.vector(Re(res_freqdom$spec.density$operators)))
  plot(x = as.vector(Im(res_dpca$spectrum)), y = as.vector(Im(res_freqdom$spec.density$operators)))
  lm(as.vector(Im(res_dpca$spectrum)) ~ as.vector(Im(res_freqdom$spec.density$operators)))
})



test_that("Test dpca, stepwise", {

  set.seed(123456)
  ## simulate some data
  nrx <- 100L
  ncx <- 1000L
  q <- 4L

  epsilon <- matrix(rnorm(ncx * q), nrow = q)

  b_filter <- vapply(1:10, function(l) {
    matrix(rnorm(q * nrx, sd = 1/l), q, nrx)
  }, matrix(0, q, nrx))

  chi <- .Call("R_filter_process", b_filter, epsilon, as.integer(1:10),
               nrx, q, q, ncx, 10L, 1L, 0L, 0L)

  x <- chi + rnorm(nrx * ncx, sd = 0.1 * sd(chi))
  bw <- as.integer(floor(ncol(x)^(1/3)))
  weights <- weights.Bartlett(-bw:bw/bw)
  system.time(res_dpca <- .Call("R_dpca", x, 4L, -10:10/10 * pi, bw, .Machine$double.eps, weights))
  system.time(res_freqdom <- freqdom::dpca(t(x), bw, -10:10/10 * pi, 4L))

  ones <- rep(1, length(weights))

  ## 1. compute covariance structure without weights
  covs1 <- .Call("R_lagged_covs", x, x, -bw:bw, nrx, ncx, nrx, ncx, ones, 1L)
  covs2 <- cov.structure(t(x), t(x), -bw:bw)$operators
  ## covs2[,,2]/covs1[,,2]
  ## covs2[,,16]/covs1[,,16]
  ## covs2[,,17]/covs1[,,17]
  expect_lt(sum(vapply(seq_len(length(weights)), function(i) sum((covs1[,,i] - covs2[,,i])^2), 0)), 1e-12)

  ## 2. compute covariance structure with weights
  covs1 <- .Call("R_lagged_covs", x, x, -bw:bw, nrx, ncx, nrx, ncx, weights, 1L)
  tmp <- cov.structure(t(x), t(x), -bw:bw)$operators
  covs2 <- vapply(1:dim(tmp)[3], function(i) tmp[,,i] * weights[i],
                  matrix(0, dim(tmp)[1], dim(tmp)[1]))
  expect_lt(sum((covs1 - covs2)^2), 1e-12)


  ## 3. Compute spectrum
  Ch <- timedom(covs2, lags = -bw:bw)
  dft1 <- fourier.transform(Ch, -10:10/10 * pi)$operators
  dft2 <- .Call("R_fourier_transform", as.vector(covs1), dim(covs1)[1], dim(covs1)[2],
                -10:10/10 * pi, 21L, -bw:bw, length(weights))
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
  abs(sum((spec1 - spec2)^2))
  tmp1 <- lapply(seq_len(dim(spec1)[3]),
                 function(i) .Call("R_arnoldi_eigs", mat = spec1[,,i], dim = nrx, 4L,
                                   .Machine$double.eps, 1L, 0L, 1L, 1L))
  tmp2 <- freqdom:::freqdom.eigen(res_freqdom$spec.density)

  ## eigenvalues
  res1 <- vapply(tmp1, function(x) Re(x$values), numeric(4))
  res2 <- tmp2$values[1:4,]
  expect_lt(sum((res1 - res2)^2), 1e-12)

  ## eigenvectors
  ev1 <- vapply(tmp1, function(x) x$vectors, matrix(1+0i, dim(tmp1[[1]]$vectors)[1], dim(tmp1[[1]]$vectors)[2]))
  ev2 <- if (dim(tmp1[[1]]$vectors)[1] == q) {
            tmp2$vectors[1:4,,]
          } else  {
            vapply(seq_len(dim(tmp2$vectors)[3]),
                   function(i) t(tmp2$vectors[1:4,,i]),
                   matrix(1+0i, nrx, q))
          }
  diff <- vapply(seq_len(dim(ev1)[3]), function(i) sum((abs(ev1[,,i]/ev2[,,i]) - 1)^2), 0 + 0i)
  expect_lt(sum(abs(diff)), 1e-12)



  ## resulting filters
  f1 <- fourier.inverse(freqdom(tmp2$vectors[1:4,,], freq = -10:10/10*pi), -bw:bw)$operators
  f2 <- .Call("R_fourier_inverse", tmp2$vectors[1:4,,], 4L, 100L,
              as.integer(-bw:bw), length(-bw:bw), -10:10/10 * pi, 21L)
  f3 <- .Call("R_fourier_inverse", ev1, 4L, 100L,
              as.integer(-bw:bw), length(-bw:bw), -10:10/10 * pi, 21L)
  expect_lt( sum((f1 - f2)^2) , 1e-12)
  expect_lt( sum((f1 - f3)^2) , 1e-12)


  ## approach 2: compute p~ * p and apply fourier transformation thereafter.
  k_theta1 <- vapply( seq_len(dim(ev1)[3]), function(i) t(Conj(ev1[,,i])) %*% ev1[,,i], matrix(0 + 0i, nrx, nrx))
  k_theta2 <- vapply( seq_len(dim(ev2)[3]), function(i) t(Conj(ev2[,,i])) %*% ev2[,,i], matrix(0 + 0i, nrx, nrx))
  sum((k_theta1 - k_theta2)^2)

  ff1 <- .Call("R_fourier_inverse", k_theta1, 100L, 100L,
              as.integer(-bw:bw), length(-bw:bw), -10:10/10 * pi, 21L)

  expect_equal(ff1, res_dpca$filter$filter_dcc)

  ## SEXP R_filter_process(SEXP r_f, SEXP r_x, SEXP r_lags,
  ##                     SEXP r_nrf, SEXP r_ncf, SEXP r_nrx, SEXP r_ncx,
  ##                     SEXP r_nlags, SEXP r_inx, SEXP r_transf, SEXP r_rev)
  dcc1 <- .Call("R_filter_process", ff1, x, as.integer(-bw:bw), nrx, nrx, nrx, ncx,
               length(-bw:bw), 1L, 0L, 0L)
  dcc2 <- res_freqdom$Xhat
  dcc3 <- res_dpca$dcc
  expect_equal(dcc1, dcc3)

  plot(x = as.vector(x), y = as.vector(t(dcc2)))
  points(x = as.vector(x), y = as.vector(dcc3), col = "red")

  summary(lm(as.vector(x) ~as.vector(dcc3)))
  summary(lm(as.vector(x) ~as.vector(t(dcc2))))

  plot(x = as.vector(chi), y = as.vector(t(dcc2)))
  points(x = as.vector(chi), y = as.vector(dcc1), col = "red")

  idx <- abs(x)> 3 & abs(dcc1) < 0.2


  res1 <- res_dpca$filters
  res2 <- .Call("R_fourier_inverse", )
  res2 <- vapply(seq_len(bw*2+1),
                 function(i) t(res_freqdom$filters$operators[,,i]),
                matrix(0, nrow = nrx, ncol = q))

  expect_lt( sum((res1-res2)^2) , 1e-12)


  m1 <- matrix(runif(10) + runif(10) * 1i, ncol = 2)
  m1 %*% Conj(t(m1))




})
