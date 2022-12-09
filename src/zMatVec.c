#include <R.h>
#include <Rinternals.h>
#include <R_ext/Lapack.h>
#include <complex.h>
#include "dpca.h"


void zMatVec(double _Complex *x, double _Complex* y, Rcomplex* mat, int dim) {

  int i, j;

  for (i = 0; i < dim; i++) {
    double _Complex accum = 0.0 + 0.0 * _Complex_I;
    for (j = 0; j < dim; j++) {
      accum += x[j] * (mat[j * dim + i].r + mat[j * dim + i].i * _Complex_I);
      /* printf("Curr mat[%d,%d]: %f%+f\n", i,j, creal(mat[i * dim + j]), cimag(mat[i * dim + j])); */
      /* printf("Curr v[%d]: %f%+f\n", j, creal(x[j]), cimag(x[j])); */
      /* printf("Curr sum: %f%+f\n", creal(x[j] * mat[i * dim + j]), cimag(x[j] * mat[i * dim + j])); */
      /* printf("Curr accum: %f%+f\n", creal(accum), cimag(accum)); */
    }
    y[i] = accum;
  }

}

void zMatVecLa(double _Complex *x, double _Complex* y, Rcomplex* mat, int dim, int ldm) {

  Rcomplex alpha;
  alpha.r = 1.0; alpha.i = 0.0;
  Rcomplex beta;
  beta.r = 0.0; beta.i = 0.0;
  int inc = 1;

  F77_CALL(zgemv)("N", &dim, &dim, &alpha, mat,
                  &ldm,
                  (Rcomplex*) x, &inc, &beta,
                  (Rcomplex*) y, &inc);
}




SEXP R_zMatVec(SEXP r_mat, SEXP r_vec, SEXP r_dim, SEXP version) {
  int dim = *INTEGER(r_dim);
  SEXP res = PROTECT(allocVector(CPLXSXP, dim));
  int ldm = nrows(r_mat);

  if (ldm < dim) {
    error("ldm must be larger than dim!");
  }

  if (*INTEGER(version) == 1) {

    double _Complex x[dim];
    double _Complex y[dim];

    // printf("dim: %d\n", dim);
    for (int i = 0; i < dim; i ++)
      x[i] = COMPLEX(r_vec)[i].r + COMPLEX(r_vec)[i].i * _Complex_I;

    zMatVec(x, y, COMPLEX(r_mat), dim);
    for (int i = 0; i < dim; i ++){
      COMPLEX(res)[i].r = creal(y[i]);
      COMPLEX(res)[i].i = cimag(y[i]);

    }
  } else {
    // printf("blube\n");
    zMatVecLa((double _Complex *) COMPLEX(r_vec), (double _Complex *) COMPLEX(res), COMPLEX(r_mat), dim, ldm);
  }

  UNPROTECT(1);
  return res;
}
