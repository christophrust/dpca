#ifndef  USE_FC_LEN_T
# define USE_FC_LEN_T
#endif
#include <Rconfig.h>
#include <Rinternals.h>
#include <complex.h>
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


SEXP R_complex_crossprod(SEXP r_x, SEXP r_trans_conj) {

    int nrx = nrows(r_x);
    int ncx = ncols(r_x);
    int trans_conj = *INTEGER(r_trans_conj);

    int dim;
    if (trans_conj) {
        dim = nrx;
    } else {
        dim = ncx;
    }

    SEXP res = PROTECT(allocMatrix(CPLXSXP, dim, dim));

    complex_crossprod((_Complex double *) COMPLEX(r_x), nrx, ncx, (_Complex double *) COMPLEX(res), trans_conj);

    UNPROTECT(1);
    return res;
}
