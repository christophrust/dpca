library(dpca)


context("Hallin & Liska (2007) selection criterium")

test_that("C function hl_select_q", {

  dim <- 150L
  m <- matrix(cos(runif(dim^2, -pi, pi)) + 1i * sin(runif(dim^2, -pi, pi)), ncol = dim)
  spec <- m %*% t(Conj(m))
  max_q <- 12L
  nfreqs <- 10L
  spec_a <- array(dim = c(dim,dim, nfreqs))
  spec_a[] <- spec
  n_path <- as.integer(seq(90, 150, 10))

  r1 <- .Call("R_hl_select_q",spec_a, n_path, max_q, dim, nfreqs, 1L, .Machine$double.eps)
  r2 <- vapply(n_path, function(n) {
    e_dec <- eigen(spec[seq_len(n), seq_len(n)])
    (sum(Re(diag(spec)[seq_len(n)])) - cumsum(Re(e_dec$values[seq_len(max_q)])))/(n)
  }, numeric(max_q))

  expect_equal(r1, r2)
})