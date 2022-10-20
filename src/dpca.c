#include "dpca.h"
#include "R_ext/RS.h"
#include "R_ext/Rdynload.h"
#include "Rinternals.h"


SEXP R_dpca(SEXP r_x, SEXP r_q, SEXP r_freqs, SEXP r_bandwidth) {

    int nrx = nrows(r_x);
    int ncx = ncols(r_x);
    int nfreqs = length(r_freqs);
    int bw = *INTEGER(r_bandwidth);
    int nlags = 2 * bw + 1;
    int lags[nlags];
    for (int i = 0; i < nlags; i++)
        lags[i] = i-bw;

    double *covs;
    covs = (double *) R_Calloc(nrx * nrx * nlags, double);
    double _Complex * spec;
    spec = (double _Complex *) R_Calloc(nrx * nrx * nfreqs, double _Complex);

    /* compute autocovariances */
    lagged_covs(REAL(r_x), REAL(r_x), covs, lags, nlags,
                nrx, ncx, nrx, ncx);

    /* compute spectrum */
    fourier_transform(covs, nrx, nrx, REAL(r_freqs), nfreqs, lags, nlags, spec);

    /* eigen decomposition of spectrum */



}
