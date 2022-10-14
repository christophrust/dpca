#ifndef DPCA_H_
#define DPCA_H_

#include <R.h>
#include <Rinternals.h>
#include <complex.h>
#include <arpack.h>
#include <stat_c.h>
#include <debug_c.h>
#include <R_ext/Lapack.h>
#include <math.h>

SEXP R_arnoldi_eigs(SEXP r_mat, SEXP r_dim, SEXP r_q, SEXP r_tol, SEXP r_normalize_evals, SEXP r_verbose);

SEXP R_zMatVec(SEXP r_mat, SEXP r_vec, SEXP r_dim, SEXP version);

void zMatVecLa(double _Complex *x, double _Complex* y, Rcomplex* mat, int dim);

void zMatVec(double _Complex *x, double _Complex* y, Rcomplex* mat, int dim);

#endif // DPCA_H_
