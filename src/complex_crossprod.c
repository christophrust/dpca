#include "complex_crossprod.h"

#ifndef  USE_FC_LEN_T
# define USE_FC_LEN_T
#endif

#include <Rconfig.h>
#include <Rinternals.h>
#include <R_ext/BLAS.h>

#ifndef FCONE
# define FCONE
#endif

void complex_crossprod(double _Complex *x, int nrx, int ncx, double _Complex *res, int trans_conj) {

    int m, n, k, lda, ldc;
    double _Complex alpha = 1.0;
    double _Complex beta = 0;


    if (trans_conj) {
        m = nrx;
        n = nrx;
        k = ncx;
        lda = nrx;
        ldc = nrx;

        F77_CALL(zgemm)("N", "C",  &m, &n, &k, (Rcomplex *) &alpha, (Rcomplex *) x, &lda,
                        (Rcomplex *) x, &lda, (Rcomplex *) &beta, (Rcomplex *) res, &ldc FCONE FCONE);

    } else {

        m = ncx;
        n = ncx;
        k = nrx;
        lda = nrx;
        ldc = ncx;

        F77_CALL(zgemm)("C", "N",  &m, &n, &k, (Rcomplex *) &alpha, (Rcomplex *) x, &lda,
                        (Rcomplex *) x, &lda, (Rcomplex *) &beta, (Rcomplex *) res, &ldc FCONE FCONE);

    }

}

