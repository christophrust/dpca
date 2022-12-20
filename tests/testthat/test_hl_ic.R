library(dpca)


context("Hallin & Liska (2007) selection criterium")

test_that("C function hl_ic", {

  dim <- 150L
  m <- matrix(cos(runif(dim^2, -pi, pi)) + 1i * sin(runif(dim^2, -pi, pi)), ncol = dim)
  spec <- m %*% t(Conj(m))
  max_q <- 12L

  e_dec <- .Call("R_arnoldi_eigs", mat = spec, dim = dim, max_q, .Machine$double.eps, 1L, 0L, 0L, 0L)


  ## single frequency
  r1 <- .Call("R_hl_ic", spec, e_dec$values, max_q, 1L, dim,
              dim, 1L)
  r2 <- (sum(Re(diag(spec))) - c(0,cumsum(Re(e_dec$values))))/dim
  expect_equal(r1, r2)

  r1 <- .Call("R_hl_ic", spec, e_dec$values, max_q, 1L, dim,
              dim, 2L)
  r2 <- log((sum(Re(diag(spec))) - c(0, cumsum(Re(e_dec$values))))/dim)
  expect_equal(r1, r2)


  ## multi frequency
  spec_a <- array(dim = c(dim,dim, 5))
  spec_a[] <- spec
  vals_a <- matrix(NA, max_q, 5L)
  vals_a[] <- e_dec$values
  r1 <- .Call("R_hl_ic", spec_a, vals_a, max_q, 5L, dim,
              dim, 1L)
  r2 <- (sum(Re(diag(spec))) - c(0, cumsum(Re(e_dec$values))))/dim
  expect_equal(r1, r2)

  r1 <- .Call("R_hl_ic", spec, e_dec$values, max_q, 1L, dim,
              dim, 2L)
  r2 <- log((sum(Re(diag(spec))) - c(0, cumsum(Re(e_dec$values))))/dim)
  expect_equal(r1, r2)



  ## sub spectrum
  r1 <- .Call("R_hl_ic", spec_a, vals_a, max_q, 5L, as.integer(dim/2),
              dim, 1L)
  r2 <- (sum(Re(diag(spec)[seq_len(dim/2)])) - c(0, cumsum(Re(e_dec$values))))/(dim/2)
  expect_equal(r1, r2)

})
