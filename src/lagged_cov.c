#include <R.h>
#include "Rinternals.h"
#include "R_ext/RS.h"
#include <R_ext/Lapack.h>
#include "R_ext/Error.h"
#include <math.h>
#include "dpca.h"

void lagged_cov(double *x, double *y, double *res, int lag, int nrx, int ncx, int nry, int ncy) {

    double alpha = 1.0/(ncx - abs(lag));
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

    F77_CALL(dgemm)("N", "T", &m, &n, &k, &alpha, x + offsetx, &ldx, y + offsety, &ldy, &beta, res, &ldx);

}



void lagged_covs(double *x, double *y, double *res, int *lags, int nlags, int nrx, int ncx, int nry, int ncy) {

    int dim = nrx * nry;
    for (int i = 0; i < nlags; i++) {
        lagged_cov(x, y, res + i * dim, lags[i], nrx, ncx, nry, ncy);
    }

}


// here r_x and r_y are {n,m} by T arrays
SEXP R_lagged_cov(SEXP r_x, SEXP r_y, SEXP r_lag, SEXP r_nrx, SEXP r_ncx, SEXP r_nry, SEXP r_ncy) {
    int ncx = *INTEGER(r_ncx);
    int nrx = *INTEGER(r_nrx);
    int ncy = *INTEGER(r_ncy);
    int nry = *INTEGER(r_nry);
    int lag = *INTEGER(r_lag);

    if (length(r_x) != nrx * ncx)
        error("Supplied dimension of X do not fit to passed array length!");
    if (length(r_y) != nry * ncy)
        error("Supplied dimension of Y do not fit to passed array length!");
    if (ncx != ncy)
        error("X and Y must have same number of columns!");
    if (ncx - 1 <= lag)
        error("Not enough observations to compute covariance with lag %d!", lag);

    SEXP res = PROTECT(allocMatrix(REALSXP, nrx, nry));

    lagged_cov(REAL(r_x), REAL(r_y), REAL(res), lag, nrx, ncx, nry, ncy);

    UNPROTECT(1);
    return res;
}


SEXP R_lagged_covs(SEXP r_x, SEXP r_y, SEXP r_lags, SEXP r_nrx, SEXP r_ncx, SEXP r_nry, SEXP r_ncy) {
    int ncx = *INTEGER(r_ncx);
    int nrx = *INTEGER(r_nrx);
    int ncy = *INTEGER(r_ncy);
    int nry = *INTEGER(r_nry);
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

    lagged_covs(REAL(r_x), REAL(r_y), REAL(res), INTEGER(r_lags), nlags, nrx, ncx, nry, ncy);

    UNPROTECT(1);
    return res;
}
