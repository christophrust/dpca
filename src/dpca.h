#ifndef DPCA_H_
#define DPCA_H_

#include <R.h>
#include <Rinternals.h>
#include <complex.h>


SEXP R_arnoldi_eigs(SEXP r_mat, SEXP r_dim, SEXP r_q, SEXP r_tol,
                    SEXP r_normalize_evals, SEXP r_verbose,
                    SEXP r_row_evecs, SEXP r_transpose_out);

SEXP R_zMatVec(SEXP r_mat, SEXP r_vec, SEXP r_dim, SEXP version);

SEXP R_lagged_cov(SEXP r_x, SEXP r_y, SEXP r_lag, SEXP r_nrx,
                  SEXP r_ncx, SEXP r_nry, SEXP r_ncy, SEXP r_weight,
                  SEXP r_center);

SEXP R_lagged_covs(SEXP r_x, SEXP r_y, SEXP r_lags, SEXP r_nrx,
                   SEXP r_ncx, SEXP r_nry, SEXP r_ncy, SEXP r_weights,
                   SEXP r_center);

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

SEXP R_filter_process(SEXP r_f, SEXP r_x, SEXP r_lags,
                      SEXP r_nrf, SEXP r_ncf, SEXP r_nrx,
                      SEXP r_ncx, SEXP r_nlags, SEXP r_inx,
                      SEXP r_transf);
SEXP R_dpca(SEXP r_x, SEXP r_q, SEXP r_freqs, SEXP r_bandwidth,
            SEXP r_tol, SEXP kernel);

void lagged_cov(double *x, double *y, double *res, int lag,
                int nrx, int ncx, int nry, int center, double weight);

void lagged_covs(double *x, double *y, double *res, int *lags, int nlags,
                 int nrx, int ncx, int nry, int ncy, double * weights, int center);

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

void arnoldi_eigs(Rcomplex *mat, int dim, int q,
                  Rcomplex *eval, Rcomplex *evecs,
                  double tol, int normalize_evals, int verbose,
                  int row_evecs, int transpose_out);

void filter_process(double *f, double *x, int *lags, int nrf,
                    int ncf, int nrx, int ncx, int nlags, double *y,
                    int insert_sample_end, int transf, int rev);

#endif // DPCA_H_
