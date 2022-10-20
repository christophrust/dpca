#include <R.h>
#include "R_ext/RS.h"
#include "Rinternals.h"
// #include "dpca.h"
#include "R_ext/Lapack.h"

/*
** Compute the convolution of a (multiple) filters f and a (multivariate) time series x
**
** Inputs:
**  - f: array (of length nrf * ncf * nlags) holding the filter coefficients
**  - x: array (of length nrx * ncx) where each row is one univariate time series
**       and each column one observation in time
**  - lags: array of lags of length nlags
**
** Ouput:
**  - y: array of length nrf * ncx
**
** NB: it must hold that ncf = nrx
*/
void filter_process1(double *f, double *x, int *lags, int nrf, int ncf, int nrx, int ncx, int nlags, double *y) {

    int fdim12 = nrf * ncf;
    int maxlag = 0;
    int minlag = 0;

    for (int i = 0; i < nlags; i++) {
        if (lags[i] > maxlag) maxlag = lags[i];
        if (lags[i] < minlag) minlag = lags[i];
    }

    memset(y, NA_REAL, nrf * ncx * sizeof(double));
    double tempval;

    for (int t = maxlag; t < ncx + minlag; t++) {
        for (int i=0; i < nrf; i++) {
            tempval = 0;
            for (int j = 0; j < nrx; j++){
                for (int k = 0; k < nlags; k++) {
                    tempval += f[i + j * nrf + k * fdim12 ] * x[j + nrx * (t - lags[k])];
                }
            }
            y[i + nrf * t] = tempval;
        }
    }
}


void filter_process(double *f, double *x, int *lags, int nrf, int ncf,
                    int nrx, int ncx, int nlags, double *y, int insert_sample_end,
                    int transf, int rev) {

    int fdim12 = nrf * ncf;
    int maxlag = lags[0];
    int minlag = lags[0];

    for (int i = 0; i < nlags; i++) {
        if (lags[i] > maxlag) maxlag = lags[i];
        if (lags[i] < minlag) minlag = lags[i];
    }

    int maxlag0 = 0;
    int minlag0 = 0;
    if (maxlag > 0) maxlag0 = maxlag;
    if (minlag < 0) minlag0 = minlag;

    int nry = nrf;
    if (transf) nry = ncf;

    memset(y, 0, nry * ncx * sizeof(double));

    int n, offsetf, offsetx, offsety;
    double alpha = 1.0;
    double beta = 1.0;

    offsety = nry * maxlag0;
    n = ncx - maxlag0 + minlag0;

    for (int k = 0; k < nlags; k++) {


        if (rev) {
            offsetf = (nlags - k) * fdim12;
        } else {
            offsetf = k * fdim12;
        }

        offsetx = nrx * (maxlag -lags[k]);

        if (transf) {
            F77_CALL(dgemm)("T", "N", &ncf, &n, &nrx, &alpha,
                            f + offsetf, &nrf, x + offsetx, &nrx,
                            &beta, y + offsety, &nry);
        } else {
            F77_CALL(dgemm)("N", "N", &nrf, &n, &nrx, &alpha,
                            f + offsetf, &nrf, x + offsetx, &nrx,
                            &beta, y + offsety, &nry);
        }
    }

    if (insert_sample_end) {
        // simple method: replace NA colums with row means of y
        double row_means[nrf];


        for (int i = 0; i < nry; i++) {
            row_means[i] = 0;
            for (int j = 0; j < n; j++)
                row_means[i] += y[offsety + i + j * nry] / n;
        }


        for (int i = 0; i < nry; i++) {
            for (int k = 0; k < maxlag0; k++) {
                y[i + k * nry] = row_means[i];
            }
            for (int k = 0; k < abs(minlag0); k++) {
                y[offsety + i + (n + k) * nry] = row_means[i];
            }
        }
    }
}



SEXP R_filter_process(SEXP r_f, SEXP r_x, SEXP r_lags,
                      SEXP r_nrf, SEXP r_ncf, SEXP r_nrx, SEXP r_ncx,
                      SEXP r_nlags, SEXP r_inx, SEXP r_transf, SEXP r_rev) {

    int nrf = *INTEGER(r_nrf);
    int ncf = *INTEGER(r_ncf);
    int ncx = *INTEGER(r_ncx);
    int transf = *INTEGER(r_transf);
    int rev = *INTEGER(r_rev);
    int nry = nrf;
    if (transf) nry = ncf;
    SEXP res = PROTECT(allocMatrix(REALSXP, nry , ncx));

    filter_process(REAL(r_f), REAL(r_x), INTEGER(r_lags), nrf, *INTEGER(r_ncf),
                   *INTEGER(r_nrx), *INTEGER(r_ncx), *INTEGER(r_nlags), REAL(res),
                   *INTEGER(r_inx), transf, rev);
    UNPROTECT(1);
    return res;
}
