#include "Rinternals.h"
#include "dpca.h"
#include <complex.h>


// TODO: ask akif, if multiple access into array on heap for accumulation.
// better use value on stack, accumulate there and assign to heap?
void fourier_inverse1(double _Complex * f, int nrf, int ncf,
                     int *lags, int nlags, double *freqs, int nfreqs,
                     double *res, double * cmplx_accum) {

    int flen = nrf * ncf;
    memset(res, 0, nrf * ncf * nlags * sizeof(double));

    // temporary object
    double _Complex z;
    double _Complex z2;
    *cmplx_accum = 0.0;

    for (int i = 0; i < nlags; i++) {
        for (int j = 0; j < nfreqs; j++) {

            z = cexp(freqs[j] * (double) lags[i] * _Complex_I);

            for (int k = 0; k < flen; k++) {
                z2 = z * f[j * flen + k];
                res[i * flen + k] += creal(z2) / (double) nfreqs;
                *cmplx_accum += cimag(z2);
            }
        }
    }
}


void fourier_inverse(double _Complex * f, int nrf, int ncf,
                     int *lags, int nlags, double *freqs, int nfreqs,
                     double *res, double * cmplx_accum) {

    int flen = nrf * ncf;
    memset(res, 0, nrf * ncf * nlags * sizeof(double));

    // temporary objects
    double * acc_flen;
    acc_flen = (double *) Calloc(flen, double);

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
        for (int k = 0; k< flen; k++ ) *cmplx_accum += fabs(acc_flen[k]);
    }
}


/* SEXP R_fourier_inverse(SEXP r_f, SEXP r_nrf, SEXP r_ncf, */
/*                      SEXP r_lags, SEXP r_nlags, */
/*                      SEXP r_freqs, SEXP r_nfreqs) { */

/*     int nrf = *INTEGER(r_nrf); */
/*     int ncf = *INTEGER(r_ncf); */
/*     int nlags = *INTEGER(r_nlags); */
/*     int nfreqs = *INTEGER(r_nfreqs); */
/*     double accum; */


/*     SEXP res = PROTECT(alloc3DArray(REALSXP, nrf, ncf, nlags)); */

/*     fourier_inverse((double _Complex *) COMPLEX(r_f), nrf, ncf, */
/*                     INTEGER(r_lags), nlags, REAL(r_freqs), nfreqs, REAL(res), &accum); */

/*     printf("Accum %f\n", accum); */

/*     if (accum > 1e-9) { */
/*         warning("The imaginary part of the coefficients was not zero, probably due to an assymmetric spectrum!"); */
/*     } */

/*     UNPROTECT(1); */
/*     return res; */
/* } */

SEXP R_fourier_inverse(SEXP r_f, SEXP r_nrf, SEXP r_ncf,
                     SEXP r_lags, SEXP r_nlags,
                     SEXP r_freqs, SEXP r_nfreqs) {

    int nrf = *INTEGER(r_nrf);
    int ncf = *INTEGER(r_ncf);
    int nlags = *INTEGER(r_nlags);
    int nfreqs = *INTEGER(r_nfreqs);
    double accum;


    SEXP res = PROTECT(alloc3DArray(REALSXP, nrf, ncf, nlags));

    fourier_inverse((double _Complex *) COMPLEX(r_f), nrf, ncf,
                    INTEGER(r_lags), nlags, REAL(r_freqs), nfreqs, REAL(res), &accum);

    Rprintf("Accum %f\n", accum);

    if (accum > 1e-9) {
        warning("The imaginary part of the coefficients was not zero, probably due to an assymmetric spectrum!");
    }

    UNPROTECT(1);
    return res;
}





SEXP R_fourier_inverse1(SEXP f, SEXP dim_f1, SEXP dim_f2,
                     SEXP lags, SEXP n_lags,
                     SEXP freqs, SEXP n_freqs) {

    int res_array_len = *INTEGER(dim_f1) * *INTEGER(dim_f2) * *INTEGER(n_lags);
    int sub_dim = *INTEGER(dim_f1) * *INTEGER(dim_f2);

    // initialize res object
    SEXP res = PROTECT(alloc3DArray(REALSXP, *INTEGER(dim_f1), *INTEGER(dim_f2), *INTEGER(n_lags)));
    double complex *tmp_cmplx_array;
    tmp_cmplx_array = (double complex *) Calloc(res_array_len, double complex);
    for (int i = 0; i < res_array_len; i++) {
        tmp_cmplx_array[i] = 0.0 + 0.0 * I;
    }


    // temporary object
    double complex z;

    for (int i = 0; i < *INTEGER(n_lags); i++){
        for (int j = 0; j < *INTEGER(n_freqs); j++) {
            z = cexp(REAL(freqs)[j] * (double)INTEGER(lags)[i] * I);
            for (int k = 0; k < sub_dim; k++) {
                tmp_cmplx_array[i * sub_dim + k] +=
                    z * (COMPLEX(f)[j * sub_dim + k].r + COMPLEX(f)[j * sub_dim + k].i * I);
            }
        }
    }

    //check that all values of tmp_cmplx_array are real and copy into res
    double accum = 0.0;
    for (int i = 0; i < res_array_len; i++){
        accum += fabs(cimag(tmp_cmplx_array[i]));
        REAL(res)[i] = creal(tmp_cmplx_array[i])/(double) *INTEGER(n_freqs);
    }

    Rprintf("Accum %f\n", accum);
    if (accum > 1e-9) {
        warning("The imaginary part of the coefficients was not zero, probably due to an assymmetric spectrum!");
        Rprintf("The accumulated value of the complex part is: %f \n", accum);
    }

    Free(tmp_cmplx_array);

    UNPROTECT(1);
    return res;
}
