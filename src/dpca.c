#include "dpca.h"
#include "R_ext/RS.h"
#include "R_ext/Rdynload.h"
#include "Rinternals.h"


SEXP R_dpca(SEXP r_x, SEXP r_q, SEXP r_freqs, SEXP r_bandwidth, SEXP r_tol) {

    int nrx = nrows(r_x);
    int ncx = ncols(r_x);
    double * freqs = REAL(r_freqs);
    int nfreqs = length(r_freqs);
    int bw = *INTEGER(r_bandwidth);
    int nlags = 2 * bw + 1;
    int lags[nlags];
    int q = *INTEGER(r_q);
    for (int i = 0; i < nlags; i++)
        lags[i] = i-bw;
    double tol = *REAL(r_tol);


    double *covs;
    covs = (double *) R_Calloc(nrx * nrx * nlags, double);
    SEXP spec = PROTECT(alloc3DArray(CPLXSXP, nrx, nrx, nfreqs));
    SEXP evecs = PROTECT(allocMatrix(CPLXSXP, nrx, q));
    SEXP evals = PROTECT(allocVector(CPLXSXP, q));
    SEXP filters = PROTECT(alloc3DArray(REALSXP, nrx, q, nlags));
    SEXP input = PROTECT(allocMatrix(REALSXP, q, ncx));
    SEXP dcc = PROTECT(allocMatrix(REALSXP, nrx, ncx));
    SEXP dic = PROTECT(allocMatrix(REALSXP, nrx, ncx));
    double tmp_accum;

    /* compute autocovariances */
    lagged_covs(REAL(r_x), REAL(r_x), covs, lags, nlags,
                nrx, ncx, nrx, ncx);

    /* compute spectrum */
    fourier_transform(covs, nrx, nrx, freqs,
                      nfreqs, lags, nlags, (double _Complex *) COMPLEX(spec));

    /* eigen decomposition of spectrum */
    for (int i = 0; i < nfreqs; i++)
        arnoldi_eigs(COMPLEX(spec) + nrx * nrx * i, nrx, q,
                     COMPLEX(evals) + q * i, COMPLEX(evecs) + nrx * q * i,
                     tol, 1, 0);

    /* compute filter coefficients */
    fourier_inverse((double _Complex *) COMPLEX(evecs), nrx, q, lags,
                    nlags, freqs, nfreqs, REAL(filters), &tmp_accum);

    /* apply filter on output to get input */
    filter_process(REAL(filters), REAL(r_x), lags, nrx, q, nrx,
                   ncx, nlags, REAL(input), 1, 1, 0);

    /* apply filter on output to get dcc */
    filter_process(REAL(filters), REAL(input), lags, nrx, q, q,
                   ncx, nlags, REAL(dcc), 1, 0, 1);

    /* compute idiosyncratic component */
    for (int i=0; i < nrx * ncx; i++)
        REAL(dic)[i] = REAL(r_x)[i] - REAL(dcc)[i];


    // create result list objects
    SEXP eig = PROTECT(allocVector(VECSXP, 2));
    SET_VECTOR_ELT(eig, 0, evals);
    SET_VECTOR_ELT(eig, 1, evecs);
    SEXP nms_eig = PROTECT(allocVector(STRSXP, 2));
    SET_STRING_ELT(nms_eig, 0, mkChar("values"));
    SET_STRING_ELT(nms_eig, 1, mkChar("vectors"));
    setAttrib(eig, R_NamesSymbol, nms_eig);


    SEXP res = PROTECT(allocVector(VECSXP, 6));
    SET_VECTOR_ELT(res, 0, spec);
    SET_VECTOR_ELT(res, 1, eig);
    SET_VECTOR_ELT(res, 2, filters);
    SET_VECTOR_ELT(res, 3, input);
    SET_VECTOR_ELT(res, 4, dcc);
    SET_VECTOR_ELT(res, 5, dic);

    SEXP nms = PROTECT(allocVector(STRSXP, 6));
    SET_STRING_ELT(nms, 0, mkChar("spectrum"));
    SET_STRING_ELT(nms, 1, mkChar("eig"));
    SET_STRING_ELT(nms, 2, mkChar("filters"));
    SET_STRING_ELT(nms, 3, mkChar("input"));
    SET_STRING_ELT(nms, 4, mkChar("dcc"));
    SET_STRING_ELT(nms, 5, mkChar("dic"));
    setAttrib(res, R_NamesSymbol, nms);

    R_Free(covs);
    UNPROTECT(11);
    return res;
}
