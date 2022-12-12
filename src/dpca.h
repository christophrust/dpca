#ifndef DPCA_H_
#define DPCA_H_

#include <R.h>
#include <Rinternals.h>
#include <complex.h>


SEXP R_arnoldi_eigs(SEXP r_mat, SEXP r_dim, SEXP r_q, SEXP r_tol,
                    SEXP r_normalize_evecs, SEXP r_verbose,
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
                      SEXP r_transf, SEXP r_rev);

SEXP R_recursive_filter(SEXP x, SEXP f, SEXP r_dim, SEXP r_nlags );

SEXP R_dpca(SEXP r_x, SEXP r_q, SEXP r_freqs, SEXP r_bandwidth,
            SEXP r_tol, SEXP kernel, SEXP r_max_q, SEXP r_select_q,
            SEXP r_n_path);

SEXP R_complex_crossprod(SEXP r_x, SEXP r_trans_conj);

void lagged_cov(double *x, double *y, double *res, int lag,
                int nrx, int ncx, int nry, int center, double weight);

void lagged_covs(double *x, double *y, double *res, int *lags, int nlags,
                 int nrx, int ncx, int nry, int ncy, double * weights, int center);

void zMatVecLa(double _Complex *x, double _Complex* y,
               Rcomplex* mat, int dim, int ldm);

void zMatVec(double _Complex *x, double _Complex* y,
             Rcomplex* mat, int dim);

void fourier_inverse(double _Complex * f, int nrf, int ncf,
                     int *lags, int nlags, double *freqs, int nfreqs,
                     double *res, double * cmplx_accum);

void fourier_transform(double *z, int nrz, int ncz,
                       double * freqs, int nfreq,
                       int * lags, int nlags,
                       double _Complex *res);

/**
 * @brief Compute truncated eigendecomposition of a hermitian matrix.
 *
 * @param mat An array of at least length ldm \times ldm holding the
 *   matrix to be decomposed (INPUT).
 * @param dim The dimension of the matrix mat.
 * @param ldm Leading dimension of mat.
 * @param q Number of eigenpairs to be computed.
 * @param eval Array of at least lenght q which on output will contain the
 *   computed eigenvalues.
 * @param evecs Array of at least length dim times q which will hold the
 *   eigenvectors on output.
 * @param tol Tolerance passed to the ARPACK routine.
 * @param normalize_evecs If set to a value different from zero, each of the resulting
 *   eigenvectors is normalized such that the complex part of the first entry is zero
 *   (the same is done in R's base::eigen()).
 * @param verbose If different from zero, we show some messages.
 * @param row_evecs If different from zero, the eigenvalues are roweigenvectors instead
 *   of column-eigenvectors
 * @param transpose_out If different from zero, the eigenvectors leading
 *  dimension is q, otherwise ldm
 */
void arnoldi_eigs(Rcomplex *mat, int dim, int ldm, int q,
                  Rcomplex *eval, Rcomplex *evecs,
                  double tol, int normalize_evecs, int verbose,
                  int row_evecs, int transpose_out);

void filter_process(double *f, double *x, int *lags, int nrf,
                    int ncf, int nrx, int ncx, int nlags, double *y,
                    int insert_sample_end, int transf, int rev);

void complex_crossprod(double _Complex *x, int nrx, int ncx,
                       double _Complex *res, int trans_conj);

/**
 * @brief Selection of number of dynamic factors.
 *
 * @param spec The spectral density matrix (of dimension n by n by nfreqs).
 * @param evals An max_q by nfreqs array which will (on output) hold in
 * the first q by nfreqs entries the resulting eigenvalues.
 * @param evecs An n by max_q by nfreqs array which will (on output) hold
 * in the first n by q by nfreqs entries the resulting eigenvectors.
 * @param dim Dimension of the spectral density matrix.
 * @param nfreqs Number of frequencies where the spectrum is evaluated.
 * @param max_q The maximum number of dynamic factors.
 * @param select_q At this stage one of 1 or 2. 1 indicates that the IC1
 * criterion is used to select q, and 2 indicates IC2.
 * @param n_path Array holding the different values for subspectra to compute
 * the information criteria on to do fine tuning.
 * @param ln The length of n_path.
 *
 * @return The number of chosen dynamic factors (hence, q).
 * */
int hl_select_q(double _Complex * spec, double _Complex * evals, double _Complex *evecs,
                int dim, int nfreqs, int max_q, int select_q, int * n_path, int ln, double tol);

/**
 * @brief Compute the unpenalized information criterium from Hallin & Liska (2007)
 *
 * @param spec A dim by dim by nfreqs array.
 * @param evals The first max_q eigenvalues of spec.
 * @param max_q Maximum number of factors to be considered.
 * @param nfreqs Number of frequencies along the spectrum
 * @param dim Dimension of spectrum.
 * @param select_q Which information criterion (IC1 or IC2).
 * @param ic_vals (On output) the computed criteria (must be at least of length max_q).
 * */
void hl_ic(double _Complex * spec, double _Complex * evals, int max_q, int nfreqs,
           int dim, int select_q, double * ic_vals);

#endif // DPCA_H_
