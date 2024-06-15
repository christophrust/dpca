#include "dpca.h"
#include "R_ext/Complex.h"
#include "R_ext/RS.h"
#include "R_ext/Rdynload.h"
#include "Rinternals.h"
#include <complex.h>

#include "complex_crossprod.h"
#include "eigs.h"
#include "filter_process.h"
#include "fourier_inverse.h"
#include "fourier_transform.h"
#include "complex_mv_product.h"
#include "hl_find_stability_interval.h"
#include "hl_ic.h"
#include "hl_ic_n_path.h"
#include "hl_q_path.h"
#include "hl_select_q.h"
#include "lagged_cov.h"


SEXP R_dpca
(
    SEXP r_x,
    SEXP r_q,
    SEXP r_freqs,
    SEXP r_bandwidth,
    SEXP r_tol,
    SEXP kernel,
    SEXP r_max_q,
    SEXP r_select_q,
    SEXP r_n_path,
    SEXP r_t_path,
    SEXP r_penalties,
    SEXP r_penalty_scales
) {

    int nrx = Rf_nrows(r_x);
    int ncx = Rf_ncols(r_x);
    double *freqs = REAL(r_freqs);
    int nfreqs = length(r_freqs);
    int bw = *INTEGER(r_bandwidth);
    int nlags = 2 * bw + 1;
    int lags[nlags];
    int q;
    for (int i = 0; i < nlags; i++)
        lags[i] = i - bw;
    double tol = *REAL(r_tol);
    int select_q = *INTEGER(r_select_q);

    int max_q;
    if (select_q) {
        max_q = *INTEGER(r_max_q);
    }

    int *n_path = INTEGER(r_n_path);

    double *covs;
    covs = (double *) R_Calloc(nrx * nrx * nlags, double);
    SEXP spec = PROTECT(alloc3DArray(CPLXSXP, nrx, nrx, nfreqs));
    SEXP filter_input;
    SEXP filter_dcc = PROTECT(alloc3DArray(REALSXP, nrx, nrx, nlags));
    SEXP input;
    SEXP dcc = PROTECT(allocMatrix(REALSXP, nrx, ncx));
    SEXP dic = PROTECT(allocMatrix(REALSXP, nrx, ncx));
    SEXP evecs, evals, unpenalized_ic_vals, sample_var, info, q_path;
    double tmp_accum;
    double _Complex * evec_crossprod;
    evec_crossprod = (double _Complex *) R_Calloc(nrx * nrx * nfreqs, double _Complex);

    /* compute autocovariances */
    lagged_covs(REAL(r_x), REAL(r_x), covs, lags, nlags, nrx, ncx, nrx, ncx, REAL(kernel), 1);


    /* compute spectrum */
    fourier_transform(covs, nrx, nrx, freqs, nfreqs, lags, nlags,
                      (double _Complex *)COMPLEX(spec));


    if (select_q) { // chose q ( number of dynamic factors) using Hallin & Liska (2007) method

        int lps = length(r_penalty_scales);
        unpenalized_ic_vals = PROTECT(allocMatrix(REALSXP, max_q + 1, length(r_n_path)));
        sample_var = PROTECT(allocVector(REALSXP, lps));
        info = PROTECT(allocVector(INTSXP, 1));
        q_path = PROTECT(allocVector(INTSXP, lps));

        _Complex double * temp_evecs, * temp_evals;
        temp_evecs = (_Complex double *) R_Calloc(max_q * nrx * nfreqs, _Complex double);
        temp_evals = (_Complex double *) R_Calloc(max_q * nfreqs, _Complex double);
        hl_select_q((_Complex double *) COMPLEX(spec),
                    temp_evals,
                    temp_evecs,
                    nrx,
                    nfreqs,
                    max_q,
                    *INTEGER(r_select_q),
                    n_path,
                    length(r_n_path),
                    tol,
                    REAL(unpenalized_ic_vals),
                    REAL(r_penalties),
                    REAL(r_penalty_scales),
                    length(r_penalty_scales),
                    REAL(sample_var),
                    INTEGER(info),
                    &q,
                    INTEGER(q_path));

        evecs = PROTECT(alloc3DArray(CPLXSXP, q, nrx, nfreqs));
        evals = PROTECT(allocMatrix(CPLXSXP, q, nfreqs));
        filter_input = PROTECT(alloc3DArray(REALSXP, nrx, q, nlags));
        input = PROTECT(allocMatrix(REALSXP, q, ncx));
        // TODO: remove copy of evecs and evals arrays
        for (int i = 0; i < q; i++) {
            for (int j = 0; j < nfreqs; j++) {
                COMPLEX(evals)[i + j * q] = ((Rcomplex *) temp_evals)[i + j * max_q];

                for (int k = 0; k < nrx; k++) {
                    COMPLEX(evecs)[i + k * q + j * nrx * q] =
                        ((Rcomplex *) temp_evecs)[i + k * max_q + j * nrx * max_q];
                }
            }
        }
        R_Free(temp_evecs);
        R_Free(temp_evals);

    } else { // user-specfied q
        q = *INTEGER(r_q);
        evecs = PROTECT(alloc3DArray(CPLXSXP, q, nrx, nfreqs));
        evals = PROTECT(allocMatrix(CPLXSXP, q, nfreqs));
        filter_input = PROTECT(alloc3DArray(REALSXP, nrx, q, nlags));
        input = PROTECT(allocMatrix(REALSXP, q, ncx));

        /* eigen decomposition of spectrum with preselected q */
        if (nfreqs % 2) {
            // center
            int i = nfreqs / 2;
            arnoldi_eigs(
                (double _Complex *) COMPLEX(spec) + nrx * nrx * i, nrx, nrx, q,
                (double _Complex *) COMPLEX(evals) + q * i,
                (double _Complex *) COMPLEX(evecs) + nrx * q * i, tol, 1, 0, 1, 1);
        }
        for (int k = 0; k < nfreqs / 2; k++) {
            int i = nfreqs / 2 + k + nfreqs % 2;
            int i_neg = nfreqs / 2 - k - 1;

            arnoldi_eigs(
                (double _Complex *) COMPLEX(spec) + nrx * nrx * i, nrx, nrx, q,
                (double _Complex *) COMPLEX(evals) + q * i,
                (double _Complex *) COMPLEX(evecs) + nrx * q * i, tol, 1, 0, 1, 1);

            // copy eigenvalues and conjugated eigenvectors into i_neg subview;
            for (int j = 0; j < q; j++) {

                *(COMPLEX(evals) + q * i_neg + j) = *(COMPLEX(evals) + q * i + j);

                for (int l = 0; l < nrx ; l++) {
                    *(COMPLEX(evecs) + q * nrx * i_neg + nrx * j + l) = *(COMPLEX(evecs) + q * nrx * i + nrx * j + l);
                    (*(COMPLEX(evecs) + q * nrx * i_neg + nrx * j + l)).i *= -1;
                }
            }
        }
    }

    for (int i = 0; i < nfreqs; i++) {
        complex_crossprod((double _Complex *) COMPLEX(evecs) + nrx * q * i,
                          q, nrx, evec_crossprod + i * nrx * nrx, 0);
    }

    /* compute filter coefficients */
    fourier_inverse((double _Complex *)COMPLEX(evecs), nrx, q, lags, nlags, freqs,
                    nfreqs, REAL(filter_input), &tmp_accum);

    /* compute filter coefficients */
    fourier_inverse(evec_crossprod, nrx, nrx, lags, nlags, freqs,
                    nfreqs, REAL(filter_dcc), &tmp_accum);

    /* apply filter on output to get input */
    filter_process(REAL(filter_input), REAL(r_x), lags, q, q, nrx, ncx, nlags,
                   REAL(input), 1, 0, 0);

    /* apply filter on output to get dcc */
    filter_process(REAL(filter_dcc), REAL(r_x), lags, nrx, nrx, nrx, ncx, nlags,
                   REAL(dcc), 1, 0, 0);

    /* compute idiosyncratic component */
    for (int i = 0; i < nrx * ncx; i++)
        REAL(dic)[i] = REAL(r_x)[i] - REAL(dcc)[i];


    // create result list objects
    SEXP eig = PROTECT(allocVector(VECSXP, 2));
    SET_VECTOR_ELT(eig, 0, evals);
    SET_VECTOR_ELT(eig, 1, evecs);
    SEXP nms_eig = PROTECT(allocVector(STRSXP, 2));
    SET_STRING_ELT(nms_eig, 0, mkChar("values"));
    SET_STRING_ELT(nms_eig, 1, mkChar("vectors"));
    setAttrib(eig, R_NamesSymbol, nms_eig);


    SEXP filter = PROTECT(allocVector(VECSXP, 2));
    SET_VECTOR_ELT(filter, 0, filter_input);
    SET_VECTOR_ELT(filter, 1, filter_dcc);
    SEXP nms_filter = PROTECT(allocVector(STRSXP, 2));
    SET_STRING_ELT(nms_filter, 0, mkChar("filter_input"));
    SET_STRING_ELT(nms_filter, 1, mkChar("filter_dcc"));
    setAttrib(filter, R_NamesSymbol, nms_filter);

    SEXP hl_select;
    if (select_q) {
        hl_select = PROTECT(allocVector(VECSXP, 5));
        SEXP chosen_q = PROTECT(allocVector(INTSXP, 1));
        *INTEGER(chosen_q) = q;
        SET_VECTOR_ELT(hl_select, 0, unpenalized_ic_vals);
        SET_VECTOR_ELT(hl_select, 1, sample_var);
        SET_VECTOR_ELT(hl_select, 2, chosen_q);
        SET_VECTOR_ELT(hl_select, 3, info);
        SET_VECTOR_ELT(hl_select, 4, q_path);
        SEXP nms_hl_select = PROTECT(allocVector(STRSXP, 5));
        SET_STRING_ELT(nms_hl_select, 0, mkChar("unpenalized_ic_vals"));
        SET_STRING_ELT(nms_hl_select, 1, mkChar("sample_var"));
        SET_STRING_ELT(nms_hl_select, 2, mkChar("q"));
        SET_STRING_ELT(nms_hl_select, 3, mkChar("info"));
        SET_STRING_ELT(nms_hl_select, 4, mkChar("q_path"));
        setAttrib(hl_select, R_NamesSymbol, nms_hl_select);
    } else {
        hl_select = PROTECT(allocVector(VECSXP, 0));
    }

    SEXP res = PROTECT(allocVector(VECSXP, 7));
    SET_VECTOR_ELT(res, 0, spec);
    SET_VECTOR_ELT(res, 1, eig);
    SET_VECTOR_ELT(res, 2, filter);
    SET_VECTOR_ELT(res, 3, input);
    SET_VECTOR_ELT(res, 4, dcc);
    SET_VECTOR_ELT(res, 5, dic);
    SET_VECTOR_ELT(res, 6, hl_select);


    SEXP nms = PROTECT(allocVector(STRSXP, 7));
    SET_STRING_ELT(nms, 0, mkChar("spectrum"));
    SET_STRING_ELT(nms, 1, mkChar("eig"));
    SET_STRING_ELT(nms, 2, mkChar("filter"));
    SET_STRING_ELT(nms, 3, mkChar("input"));
    SET_STRING_ELT(nms, 4, mkChar("dcc"));
    SET_STRING_ELT(nms, 5, mkChar("dic"));
    SET_STRING_ELT(nms, 6, mkChar("HL_select"));
    setAttrib(res, R_NamesSymbol, nms);

    R_Free(covs);
    R_Free(evec_crossprod);

    if (select_q) {
        UNPROTECT(21);
    } else{
        UNPROTECT(15);
    }
    return res;
}
