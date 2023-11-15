#ifndef R_INTERFACES_H_
#define R_INTERFACES_H_

#include <R.h>
#include <Rinternals.h>

SEXP R_arnoldi_eigs(SEXP r_mat, SEXP r_dim, SEXP r_q, SEXP r_tol,
                    SEXP r_normalize_evecs, SEXP r_verbose,
                    SEXP r_row_evecs, SEXP r_transpose_out);

SEXP R_zMatVec(SEXP r_mat, SEXP r_vec, SEXP r_dim, SEXP version);


SEXP R_complex_crossprod(SEXP r_x, SEXP r_trans_conj);

SEXP R_filter_process(SEXP r_f, SEXP r_x, SEXP r_lags,
                      SEXP r_nrf, SEXP r_ncf, SEXP r_nrx, SEXP r_ncx,
                      SEXP r_nlags, SEXP r_inx, SEXP r_transf, SEXP r_rev);

SEXP R_fourier_inverse(SEXP r_f, SEXP r_nrf, SEXP r_ncf,
                     SEXP r_lags, SEXP r_nlags,
                     SEXP r_freqs, SEXP r_nfreqs);

SEXP R_fourier_transform(SEXP z, SEXP dim_z1, SEXP dim_z2,
                         SEXP freq, SEXP n_freq,
                         SEXP lags, SEXP n_lags);

SEXP R_find_stability_intervals(SEXP r_sample_var);


SEXP R_hl_ic(SEXP r_spec, SEXP r_evals, SEXP r_max_q, SEXP r_nfreqs, SEXP r_dim,
             SEXP r_ldm, SEXP r_select_q);


SEXP R_hl_ic_n_path(SEXP r_spec, SEXP r_n_path, SEXP r_max_q, SEXP r_dim,
                   SEXP r_nfreqs, SEXP r_select_q, SEXP r_tol);


SEXP R_hl_q_path
(
    SEXP r_unpenalized_ic_vals,
    SEXP r_max_q,
    SEXP r_penalty_scale,
    SEXP r_penalties
);


SEXP R_hl_select_q
(
    SEXP r_spec,
    SEXP r_n_path,
    SEXP r_max_q,
    SEXP r_dim,
    SEXP r_nfreqs,
    SEXP r_select_q,
    SEXP r_tol,
    SEXP r_penalties,
    SEXP r_penalty_scales
);


SEXP R_lagged_cov(SEXP r_x, SEXP r_y, SEXP r_lag, SEXP r_nrx,
                  SEXP r_ncx, SEXP r_nry, SEXP r_ncy,
                  SEXP r_weight, SEXP r_center);

SEXP R_lagged_covs(SEXP r_x, SEXP r_y, SEXP r_lags, SEXP r_nrx, SEXP r_ncx,
                   SEXP r_nry, SEXP r_ncy, SEXP r_weights, SEXP r_center);

SEXP R_recursive_filter(SEXP x, SEXP f, SEXP r_dim, SEXP r_nlags );

#endif // R_INTERFACES_H_
