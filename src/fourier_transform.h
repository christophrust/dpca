#ifndef FOURIER_TRANSFORM_H_
#define FOURIER_TRANSFORM_H_

void fourier_transform(double *z, int nrz, int ncz,
                       double * freqs, int nfreq,
                       int * lags, int nlags,
                       double _Complex *res);


#endif // FOURIER_TRANSFORM_H_
