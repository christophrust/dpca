#ifndef LAGGED_COV_H_
#define LAGGED_COV_H_

void lagged_cov(double *x, double *y, double *res, int lag,
                int nrx, int ncx, int nry, int center, double weight);

void lagged_covs(double *x, double *y, double *res, int *lags, int nlags,
                 int nrx, int ncx, int nry, int ncy, double * weights, int center);

#endif // LAGGED_COV_H_
