#ifndef DPCA_H_
#define DPCA_H_

#include <R.h>
#include <Rinternals.h>
#include <complex.h>


SEXP R_arnoldi_eigs(SEXP r_mat, SEXP r_dim, SEXP r_q, SEXP r_tol,
                    SEXP r_normalize_evals, SEXP r_verbose);

SEXP R_zMatVec(SEXP r_mat, SEXP r_vec, SEXP r_dim, SEXP version);

SEXP R_lagged_cov(SEXP r_x, SEXP r_y, SEXP r_lag, SEXP r_nrx,
                  SEXP r_ncx, SEXP r_nry, SEXP r_ncy);

SEXP R_lagged_covs(SEXP r_x, SEXP r_y, SEXP r_lags, SEXP r_nrx,
                   SEXP r_ncx, SEXP r_nry, SEXP r_ncy);

SEXP R_fourier_transform1(SEXP z, SEXP dim_z1, SEXP dim_z2,
                          SEXP freq, SEXP n_freq,
                          SEXP lags, SEXP n_lags);

SEXP R_fourier_transform(SEXP z, SEXP dim_z1, SEXP dim_z2,
                         SEXP freq, SEXP n_freq,
                         SEXP lags, SEXP n_lags);

SEXP R_fourier_inverse(SEXP r_f, SEXP r_nrf, SEXP r_ncf,
                       SEXP r_lags, SEXP r_nlags,
                       SEXP r_freqs, SEXP r_nfreqs);

SEXP R_fourier_inverse1(SEXP f, SEXP dim_f1, SEXP dim_f2,
                     SEXP lags, SEXP n_lags,
                     SEXP freqs, SEXP n_freqs);

void lagged_cov(double *x, double *y, double *res, int lag,
                int nrx, int ncx, int nry, int ncy);

void lagged_covs(double *x, double *y, double *res, int *lags, int nlags,
                 int nrx, int ncx, int nry, int ncy);

void zMatVecLa(double _Complex *x, double _Complex* y,
               Rcomplex* mat, int dim);

void zMatVec(double _Complex *x, double _Complex* y,
             Rcomplex* mat, int dim);

void fourier_inverse(double _Complex * f, int nrf, int ncf,
                     int *lags, int nlags, double *freqs, int nfreqs,
                     double *res, double * cmplx_accum);

void fourier_transform(double *z, int nrz, int ncz,
                       double * freqs, int nfreq,
                       int * lags, int nlags,
                       double _Complex *res);
#endif // DPCA_H_
