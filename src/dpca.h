#ifndef DPCA_H_
#define DPCA_H_

#include <R.h>
#include <Rinternals.h>
#include <complex.h>


SEXP R_arnoldi_eigs(SEXP r_mat, SEXP r_dim, SEXP r_q, SEXP r_tol, SEXP r_normalize_evals, SEXP r_verbose);
SEXP R_zMatVec(SEXP r_mat, SEXP r_vec, SEXP r_dim, SEXP version);
SEXP R_lagged_cov(SEXP r_x, SEXP r_y, SEXP r_lag, SEXP r_nrx, SEXP r_ncx, SEXP r_nry, SEXP r_ncy);


void lagged_cov(double *x, double *y, double *res, int lag, int nrx, int ncx, int nry, int ncy);
void zMatVecLa(double _Complex *x, double _Complex* y, Rcomplex* mat, int dim);
void zMatVec(double _Complex *x, double _Complex* y, Rcomplex* mat, int dim);

#endif // DPCA_H_
