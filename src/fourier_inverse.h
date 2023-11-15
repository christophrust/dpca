#ifndef FOURIER_INVERSE_H_
#define FOURIER_INVERSE_H_


void fourier_inverse(double _Complex * f, int nrf, int ncf,
                     int *lags, int nlags, double *freqs, int nfreqs,
                     double *res, double * cmplx_accum);

#endif // FOURIER_INVERSE_H_
