#include "R_ext/Error.h"
#include "R_ext/RS.h"
#include <R_ext/Lapack.h>
#include "Rinternals.h"
#include <math.h>
#include <stdlib.h>
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


// here r_x and r_y are n,m by T arrays
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

    SEXP res = PROTECT(allocMatrix(REALSXP, nrx, nry));

    lagged_cov(REAL(r_x), REAL(r_y), REAL(res), lag, nrx, ncx, nry, ncy);

    Rf_unprotect(1);
    return res;
}
