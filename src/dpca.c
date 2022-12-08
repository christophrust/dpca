#include "dpca.h"
#include "R_ext/RS.h"
#include "R_ext/Rdynload.h"
#include "Rinternals.h"
#include <complex.h>


SEXP R_dpca(SEXP r_x, SEXP r_q, SEXP r_freqs, SEXP r_bandwidth, SEXP r_tol, SEXP kernel, SEXP r_max_q, SEXP r_select_q) {

    int nrx = Rf_nrows(r_x);
    int ncx = Rf_ncols(r_x);
    double *freqs = REAL(r_freqs);
    int nfreqs = length(r_freqs);
    int bw = *INTEGER(r_bandwidth);
    int nlags = 2 * bw + 1;
    int lags[nlags];
    int q = *INTEGER(r_q);
    for (int i = 0; i < nlags; i++)
        lags[i] = i - bw;
    double tol = *REAL(r_tol);
    int select_q = *INTEGER(r_select_q);

    int spec_q;
    if (select_q) {
        spec_q = *INTEGER(r_max_q);
    } else {
        spec_q = q;
    }

    double *covs;
    covs = (double *)R_Calloc(nrx * nrx * nlags, double);
    SEXP spec = PROTECT(alloc3DArray(CPLXSXP, nrx, nrx, nfreqs));
    SEXP evecs = PROTECT(alloc3DArray(CPLXSXP, spec_q, nrx, nfreqs));
    SEXP evals = PROTECT(allocMatrix(CPLXSXP, spec_q, nfreqs));
    SEXP filter_input = PROTECT(alloc3DArray(REALSXP, nrx, q, nlags));
    SEXP filter_dcc = PROTECT(alloc3DArray(REALSXP, nrx, nrx, nlags));
    SEXP input = PROTECT(allocMatrix(REALSXP, q, ncx));
    SEXP dcc = PROTECT(allocMatrix(REALSXP, nrx, ncx));
    SEXP dic = PROTECT(allocMatrix(REALSXP, nrx, ncx));
    double tmp_accum;
    double _Complex * evec_cp;
    evec_cp = (double _Complex *) R_Calloc(nrx * nrx * nfreqs, double _Complex);

    /* compute autocovariances */
    lagged_covs(REAL(r_x), REAL(r_x), covs, lags, nlags, nrx, ncx, nrx, ncx, REAL(kernel), 1);


    /* compute spectrum */
    fourier_transform(covs, nrx, nrx, freqs, nfreqs, lags, nlags,
                      (double _Complex *)COMPLEX(spec));


    // TODO: compute eigendecomposition only on 0 to pi and get
    // eigenvalues for -pi to 0 by conjugating them!!

    /* eigen decomposition of spectrum */
    for (int i = 0; i < nfreqs; i++)
        arnoldi_eigs(COMPLEX(spec) + nrx * nrx * i, nrx, spec_q, COMPLEX(evals) + spec_q * i,
                     COMPLEX(evecs) + nrx * spec_q * i, tol, 1, 0, 1, 1);

    // do selection of number of eigenvalues using hallin & liska (2007) method
    if (select_q) {

    } else {

    }

    for (int i = 0; i < nfreqs; i++)
        complex_crossprod((double _Complex *) COMPLEX(evecs) + nrx * q * i,
                          q, nrx, evec_cp + i * nrx * nrx, 0);

    /* compute filter coefficients */
    fourier_inverse((double _Complex *)COMPLEX(evecs), nrx, q, lags, nlags, freqs,
                    nfreqs, REAL(filter_input), &tmp_accum);

    /* compute filter coefficients */
    fourier_inverse(evec_cp, nrx, nrx, lags, nlags, freqs,
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


    SEXP res = PROTECT(allocVector(VECSXP, 6));
    SET_VECTOR_ELT(res, 0, spec);
    SET_VECTOR_ELT(res, 1, eig);
    SET_VECTOR_ELT(res, 2, filter);
    SET_VECTOR_ELT(res, 3, input);
    SET_VECTOR_ELT(res, 4, dcc);
    SET_VECTOR_ELT(res, 5, dic);

    SEXP nms = PROTECT(allocVector(STRSXP, 6));
    SET_STRING_ELT(nms, 0, mkChar("spectrum"));
    SET_STRING_ELT(nms, 1, mkChar("eig"));
    SET_STRING_ELT(nms, 2, mkChar("filter"));
    SET_STRING_ELT(nms, 3, mkChar("input"));
    SET_STRING_ELT(nms, 4, mkChar("dcc"));
    SET_STRING_ELT(nms, 5, mkChar("dic"));
    setAttrib(res, R_NamesSymbol, nms);

    R_Free(covs);
    R_Free(evec_cp);
    UNPROTECT(14);
    return res;
}
