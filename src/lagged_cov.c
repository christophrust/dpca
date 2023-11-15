#ifndef  USE_FC_LEN_T
# define USE_FC_LEN_T
#endif
#include <Rconfig.h>
#include <R.h>
#include "Rinternals.h"
#include "R_ext/RS.h"
#include <R_ext/BLAS.h>
#include "R_ext/Error.h"
#include <math.h>
#include "dpca.h"
#ifndef FCONE
# define FCONE
#endif

void lagged_cov(double *x, double *y, double *res,
                int lag, int nrx, int ncx, int nry,
                int center, double weight) {

    //double alpha = 1.0/(ncx - abs(lag)) * weight;
    double alpha = 1.0/ncx * weight;
    double beta = 0;
    int m = nrx;
    int n = nry;
    int k = ncx - abs(lag);
    int ldx = nrx;
    int ldy = nry;
    int offsetx = nrx * lag;
    if (offsetx < 0) offsetx = 0;
    int offsety = -nry * lag;
    if (offsety < 0) offsety = 0;

    if (center) {

        int one = 1;
        double ones[ncx];
        for (int i = 0; i < ncx; i++) ones[i] = 1.0;
        double meanx[nrx];
        double meany[nry];
        double alpha1 = 1.0/ncx;
        double beta1 = 0.0;

        F77_CALL(dgemv)("N", &nrx, &ncx, &alpha1, x,
                        &nrx, ones, &one, &beta1, meanx, &one FCONE);

        F77_CALL(dgemv)("N", &nry, &ncx, &alpha1, y,
                        &nry, ones, &one, &beta1, meany, &one FCONE);

        // create centered versions of x and y
        double *xc, *yc;
        xc = (double *) R_Calloc(nrx * ncx, double);
        yc = (double *) R_Calloc(nry * ncx, double);

        for (int i = 0; i< nrx * ncx; i++)
            xc[i] = x[i] - meanx[i%nrx];

        for (int i = 0; i< nry * ncx; i++)
            yc[i] = y[i] - meany[i%nry];

        F77_CALL(dgemm)("N", "T", &m, &n, &k, &alpha, xc + offsetx,
                        &ldx, yc + offsety, &ldy, &beta, res, &ldx FCONE FCONE);
        R_Free(xc);
        R_Free(yc);
    } else {

        F77_CALL(dgemm)("N", "T", &m, &n, &k, &alpha, x + offsetx,
                        &ldx, y + offsety, &ldy, &beta, res, &ldx FCONE FCONE);
    }

}



void lagged_covs(double *x, double *y, double *res, int *lags, int nlags, int nrx, int ncx, int nry, int ncy, double * weights, int center) {

    int dim = nrx * nry;

    if (center) {
        int one = 1;
        double ones[ncx];
        for (int i = 0; i < ncx; i++) ones[i] = 1.0;
        double meanx[nrx];
        double meany[nry];
        double alpha1 = 1.0/ncx;
        double beta1 = 0.0;

        F77_CALL(dgemv)("N", &nrx, &ncx, &alpha1, x,
                        &nrx, ones, &one, &beta1, meanx, &one FCONE);

        F77_CALL(dgemv)("N", &nry, &ncx, &alpha1, y,
                        &nry, ones, &one, &beta1, meany, &one FCONE);

        // create centered versions of x and y
        double *xc, *yc;

        xc = R_Calloc(nrx * ncx, double);
        yc = R_Calloc(nry * ncx, double);

        for (int i = 0; i< nrx * ncx; i++)
            xc[i] = x[i] - meanx[i%nrx];

        for (int i = 0; i< nry * ncx; i++)
            yc[i] = y[i] - meany[i%nry];

        for (int i = 0; i < nlags; i++) {
            lagged_cov(xc, yc, res + i * dim, lags[i], nrx, ncx, nry, 0, weights[i]);
        }

        R_Free(xc);
        R_Free(yc);

    } else {

        for (int i = 0; i < nlags; i++) {
            lagged_cov(x, y, res + i * dim, lags[i], nrx, ncx, nry, 0, weights[i]);
        }
    }

}


// here r_x and r_y are {n,m} by T arrays
SEXP R_lagged_cov(SEXP r_x, SEXP r_y, SEXP r_lag, SEXP r_nrx,
                  SEXP r_ncx, SEXP r_nry, SEXP r_ncy,
                  SEXP r_weight, SEXP r_center) {
    int ncx = *INTEGER(r_ncx);
    int nrx = *INTEGER(r_nrx);
    int ncy = *INTEGER(r_ncy);
    int nry = *INTEGER(r_nry);
    int lag = *INTEGER(r_lag);
    int center = *INTEGER(r_center);

    if (length(r_x) != nrx * ncx)
        error("Supplied dimension of X do not fit to passed array length!");
    if (length(r_y) != nry * ncy)
        error("Supplied dimension of Y do not fit to passed array length!");
    if (ncx != ncy)
        error("X and Y must have same number of columns!");
    if (ncx - 1 <= lag)
        error("Not enough observations to compute covariance with lag %d!", lag);

    SEXP res = PROTECT(allocMatrix(REALSXP, nrx, nry));

    lagged_cov(REAL(r_x), REAL(r_y), REAL(res), lag, nrx, ncx, nry, center, *REAL(r_weight));

    UNPROTECT(1);
    return res;
}


SEXP R_lagged_covs(SEXP r_x, SEXP r_y, SEXP r_lags, SEXP r_nrx, SEXP r_ncx,
                   SEXP r_nry, SEXP r_ncy, SEXP r_weights, SEXP r_center) {

    int ncx = *INTEGER(r_ncx);
    int nrx = *INTEGER(r_nrx);
    int ncy = *INTEGER(r_ncy);
    int nry = *INTEGER(r_nry);
    int center = *INTEGER(r_center);
    int *lags = INTEGER(r_lags);
    int nlags = length(r_lags);
    int maxlag = 0;

    for (int i = 0; i < nlags; i++)
        if (abs(lags[i]) > maxlag) maxlag = abs(lags[i]);

    if (length(r_x) != nrx * ncx)
        error("Supplied dimension of X do not fit to passed array length!");
    if (length(r_y) != nry * ncy)
        error("Supplied dimension of Y do not fit to passed array length!");
    if (ncx != ncy)
        error("X and Y must have same number of columns!");
    if (ncx - 1 <= maxlag)
        error("Not enough observations to compute covariance with maximum lag %d!", maxlag);

    SEXP res = PROTECT(alloc3DArray(REALSXP, nrx, nry, nlags));

    lagged_covs(REAL(r_x), REAL(r_y), REAL(res), INTEGER(r_lags), nlags, nrx, ncx, nry, ncy, REAL(r_weights), center);

    UNPROTECT(1);
    return res;
}
