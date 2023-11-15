#ifndef  USE_FC_LEN_T
# define USE_FC_LEN_T
#endif

#include "lagged_cov.h"

#include "lagged_cov.h"
#include <stdlib.h>
#include <R_ext/BLAS.h>
#include <R_ext/Error.h>
#include <math.h>

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

