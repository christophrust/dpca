library(dpca)



test_that("HL: C function hl_ic_n_path", {
  dim <- 150L
  max_q <- 12L
  nfreqs <- 11L
  spec_a <- array(dim = c(dim, dim, nfreqs))
  for (i in seq_len(nfreqs / 2 + 1)) {
    m <- matrix(cos(runif(dim^2, -pi, pi)) + 1i * sin(runif(dim^2, -pi, pi)), ncol = dim)
    spec <- m %*% t(Conj(m))
    spec_a[, , i] <- spec
    if (i < 6) {
      spec_a[, , nfreqs - i + 1] <- t(spec)
    }
  }

  n_path <- as.integer(seq(90, 150, 10))

  system.time(r1 <- .Call("R_hl_ic_n_path", spec_a, n_path, max_q, dim, nfreqs, 1L, .Machine$double.eps))

  system.time(r2 <- vapply(n_path, function(n) {
    ic_freqs <- vapply(seq_len(nfreqs), function(i) {
      e_dec <- .Call("R_arnoldi_eigs",
        mat = spec_a[seq_len(n), seq_len(n), i],
        dim = as.integer(n),
        max_q,
        .Machine$double.eps, 1L, 0L, 0L, 0L
      )
      (sum(Re(diag(spec_a[, , i])[seq_len(n)])) - c(0, cumsum(Re(e_dec$values[seq_len(max_q)])))) / (n)
    }, numeric(max_q + 1))
    rowMeans(ic_freqs)
  }, numeric(max_q + 1)))

  expect_equal(r1, r2)
})
