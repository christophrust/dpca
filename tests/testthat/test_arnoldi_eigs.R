library(dpca)


dim <- 150L
m <- matrix(cos(runif(dim^2, -pi, pi)) + 1i * sin(runif(dim^2, -pi, pi)), ncol = dim)
spec <- m %*% t(Conj(m))

test_that("Eigenvalue and -vector equality", {
  ## normalize evectors to be up to sign equivalent to result of eigen
  r1 <- .Call("R_arnoldi_eigs", mat = spec, dim = dim, 4L, .Machine$double.eps, 1L, 0L, 0L, 0L)
  r2 <- eigen(spec)
  diff_cp <- abs(Re(t(Conj(r1$vectors)) %*% r2$vectors[, 1:4])) - diag(4)
  expect_lte(sum(Re(r1$values - r2$values[1:4])^2), 1e-12)
  expect_lte(sum(diff_cp^2), 1e-12)

  ## use only principal submatrix of spec
  r1 <- .Call("R_arnoldi_eigs", mat = spec, dim = dim - 10L, 4L, .Machine$double.eps, 1L, 0L, 0L, 0L)
  r2 <- eigen(spec[1:(dim - 10), 1:(dim - 10)])
  diff_cp <- abs(Re(t(Conj(r1$vectors)) %*% r2$vectors[, 1:4])) - diag(4)
  expect_lte(sum(Re(r1$values - r2$values[1:4])^2), 1e-12)
  expect_lte(sum(diff_cp^2), 1e-12)
})


test_that("eigs API", {
  ## compute column eigenvectors given that spectrum is hermitian spectrum
  res <- .Call("R_arnoldi_eigs", mat = spec, dim = dim, 4L, .Machine$double.eps, 1L, 0L, 0L, 0L)
  for (i in seq_len(4)) {
    expect_lte(sum(abs(spec %*% res$vectors[, i] - res$values[i] * res$vectors[, i])), 1e-10)
  }

  ## compute row eigenvectors given that spectrum is hermitian spectrum
  res <- .Call("R_arnoldi_eigs", mat = spec, dim = dim, 4L, .Machine$double.eps, 1L, 0L, 1L, 0L)
  for (i in seq_len(4)) {
    expect_lte(sum(abs(res$vectors[, i] %*% spec - res$values[i] * res$vectors[, i])), 1e-10)
  }

  ## output transposition

  ## compute column eigenvectors given that spectrum is hermitian spectrum
  res <- .Call("R_arnoldi_eigs", mat = spec, dim = dim, 4L, .Machine$double.eps, 1L, 0L, 0L, 1L)
  for (i in seq_len(4)) {
    expect_lte(sum(abs(spec %*% res$vectors[i, ] - res$values[i] * res$vectors[i, ])), 1e-10)
  }
})
