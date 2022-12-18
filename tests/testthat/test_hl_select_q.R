library(dpca)


context("Hallin & Liska (2007) selection criterium")

test_that("C function hl_select_q", {

  dim <- 150L
  max_q <- 12L
  nfreqs <- 10L
  spec_a <- array(dim = c(dim,dim, nfreqs))
  for (i in seq_len(nfreqs)) {
      m <- matrix(cos(runif(dim^2, -pi, pi)) + 1i * sin(runif(dim^2, -pi, pi)), ncol = dim)
      spec <- m %*% t(Conj(m))
      spec_a[,,i] <- spec
  }

  n_path <- as.integer(seq(90, 150, 10))
  penalty_scales <- seq(0,2, by = 0.01)
  penalties <- runif(length(n_path), min = 10, max = 12)
  system.time(r1 <- .Call("R_hl_select_q",
                          spec_a, n_path, max_q, dim, nfreqs, 1L, .Machine$double.eps,
                          penalties, penalty_scales))
  str(r1)

  system.time(r2 <- vapply(n_path, function(n) {
    ic_freqs <- vapply(seq_len(nfreqs), function(i) {
      e_dec <- .Call("R_arnoldi_eigs", mat = spec_a[seq_len(n), seq_len(n), i],
                     dim = as.integer(n),
                     max_q,
                     .Machine$double.eps, 1L, 0L, 0L, 0L)
      (sum(Re(diag(spec_a[,,i])[seq_len(n)])) - cumsum(Re(e_dec$values[seq_len(max_q)])))/(n)
    }, numeric(max_q))
    rowMeans(ic_freqs)
  }, numeric(max_q)))

  expect_equal(r1, r2)
})
