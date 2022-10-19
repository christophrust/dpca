#include "dpca.h"
#include <complex.h>
#include <string.h>


/*
  Discrete fourier transform, moving heavy stuff of fourier.transform into C code

  Inputs:
    - z: an (C-) Array of (implicit) dimension dim_z1 times dim_z2 times n_lags
    - freq: a vector of real-valued frequencies (values on [-pi, pi])
    - lags: an integer vector of lags

  Returns:
    A numeric vector of length dim_z1 * dim_z2 * n_freq

*/



void fourier_transform(double *z, int nrz, int ncz,
                       double * freqs, int nfreq,
                       int * lags, int nlags,
                       double _Complex *res) {

    double _Complex w;
    int sub_dim = nrz * ncz;
    memset(res, 0, nrz * ncz * nfreq * sizeof(double _Complex));

    for (int i = 0; i< nfreq; i++) {
        for (int j = 0; j < nlags; j++) {
            w = cexp(0.0 - lags[j] * freqs[i] * I);

            for (int k = 0; k < sub_dim; k++) {
                res[i * sub_dim + k] += w * z[j * sub_dim + k];
            }
        }
    }
}

SEXP R_fourier_transform(SEXP z, SEXP dim_z1, SEXP dim_z2,
                       SEXP freq, SEXP n_freq,
                         SEXP lags, SEXP
                         n_lags)
{


    SEXP res = PROTECT(alloc3DArray(CPLXSXP, *INTEGER(dim_z1), *INTEGER(dim_z2), *INTEGER(n_freq)));

    fourier_transform(REAL(z), *INTEGER(dim_z1), *INTEGER(dim_z2),
                       REAL(freq), *INTEGER(n_freq),
                       INTEGER(lags), *INTEGER(n_lags),
                       (double _Complex *) COMPLEX(res));

    UNPROTECT(1);
    return res;
}
