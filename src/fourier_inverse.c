#include "fourier_inverse.h"

#include <complex.h>
#include <R.h>


void fourier_inverse(double _Complex * f, int nrf, int ncf,
                     int *lags, int nlags, double *freqs, int nfreqs,
                     double *res, double * cmplx_accum) {

    int flen = nrf * ncf;
    memset(res, 0, nrf * ncf * nlags * sizeof(double));

    // temporary objects
    double * acc_flen;
    acc_flen = (double *) R_Calloc(flen, double);

    double _Complex z;
    double _Complex z2;

    *cmplx_accum = 0.0;

    for (int i = 0; i < nlags; i++) {
        memset(acc_flen, 0 , flen * sizeof(double));
        for (int j = 0; j < nfreqs; j++) {

            z = cexp(freqs[j] * (double) lags[i] * _Complex_I);

            for (int k = 0; k < flen; k++) {

                z2 = z * f[j * flen + k];
                res[i * flen + k] += creal(z2) / (double) nfreqs;
                acc_flen[k] += cimag(z2);
            }

        }
        for (int k = 0; k < flen; k++ ) *cmplx_accum += fabs(acc_flen[k]);
    }
    R_Free(acc_flen);
}
