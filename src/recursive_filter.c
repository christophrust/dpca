#ifndef  USE_FC_LEN_T
# define USE_FC_LEN_T
#endif
#include <Rconfig.h>
#include "dpca.h"
#include <R.h>
#include "R_ext/RS.h"
#include "Rinternals.h"
// #include "dpca.h"
#include "R_ext/Lapack.h"
#ifndef FCONE
# define FCONE
#endif


SEXP R_recursive_filter(SEXP x, SEXP f, SEXP r_dim, SEXP r_nlags ) {

    int nlags = *INTEGER(r_nlags);
    int dim = *INTEGER(r_dim);
    int ncx = ncols(x);
    SEXP res = PROTECT(allocMatrix(REALSXP, dim, ncx));
    memcpy(REAL(res), REAL(x), dim*ncx* sizeof(double));

    int offsetx, offsety, offsetf;
    int one = 1;
    double alpha = 1.0;
    double beta = 1.0;
    int i,j;
    for (i = 0; i < ncx; i++) {
        offsety = i * dim;

        for (j = 0; j < nlags && j < i; j++) {
            offsetf = j * dim * dim;
            offsetx = (i - j - 1) * dim;
            //printf("%i; %i; %i; %i\n", i, j, offsetx, offsetf);
            F77_CALL(dgemv)("N", &dim, &dim, &alpha,
                            REAL(f) + offsetf, &dim, REAL(res) + offsetx, &one, &beta,
                            REAL(res) + offsety, &one FCONE);
        }
    }

    UNPROTECT(1);
    return res;
}
