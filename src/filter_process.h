#ifndef FILTER_PROCESS_H_
#define FILTER_PROCESS_H_

void filter_process(double *f, double *x, int *lags, int nrf, int ncf, int nrx,
                    int ncx, int nlags, double *y, int insert_sample_end,
                    int transf, int rev);

#endif // FILTER_PROCESS_H_
