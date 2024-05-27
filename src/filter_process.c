#include "filter_process.h"
#include <stdlib.h>

#ifndef  USE_FC_LEN_T
# define USE_FC_LEN_T
#endif

#include <Rconfig.h>
#include <R_ext/Lapack.h>

#ifndef FCONE
# define FCONE
#endif


/**
* Multivariate linear filter
*
* Compute the convolution of a (multiple) filters f and a (multivariate) time series x.
*
* @param f Array (of length nrf * ncf * nlags) holding the filter coefficients
* @param x Array (of length nrx * ncx) where each row is one univariate time series
*       and each column one observation in time
* @param lags Array of lag indices
* @param nrf Number of rows in f.
* @param ncf Number of columns in f.
* @param nrx Number of rows in x.
* @param ncf Number of columns in x.
* @param nlags Length of lags.
* @param y Pointer to memory at least of size nrf * ncx * sizeof(double).
* @param insert_sample_end If 1, the NA columns at start and/or end of
* time period will be extrapolated with the mean of each univariate series. Otherwise,
* these values will remain unchanged.
* @param transf If 1, f must be passed in transposed form.
* @param rev If 1, f must be passed in reversed order wrt
* the last dimension (time periods).
*
*
* NB: it is assumed that ncf = nrx
*/
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
            offsetf = (nlags - 1 - k ) * fdim12;
        } else {
            offsetf = k * fdim12;
        }

        offsetx = nrx * (maxlag - lags[k]);

        if (transf) {
            F77_CALL(dgemm)("T", "N", &ncf, &n, &nrx, &alpha,
                            f + offsetf, &nrf, x + offsetx, &nrx,
                            &beta, y + offsety, &nry FCONE FCONE);
        } else {
            F77_CALL(dgemm)("N", "N", &nrf, &n, &nrx, &alpha,
                            f + offsetf, &nrf, x + offsetx, &nrx,
                            &beta, y + offsety, &nry FCONE FCONE);
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
