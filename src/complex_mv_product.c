#include "complex_mv_product.h"

#ifndef  USE_FC_LEN_T
# define USE_FC_LEN_T
#endif
#include <Rconfig.h>
#include <R_ext/Lapack.h>


#ifndef FCONE
# define FCONE
#endif


void zMatVecLa(double _Complex *x, double _Complex* y, Rcomplex* mat, int dim, int ldm) {

  Rcomplex alpha;
  alpha.r = 1.0; alpha.i = 0.0;
  Rcomplex beta;
  beta.r = 0.0; beta.i = 0.0;
  int inc = 1;

  F77_CALL(zgemv)("N", &dim, &dim, &alpha, mat,
                  &ldm,
                  (Rcomplex*) x, &inc, &beta,
                  (Rcomplex*) y, &inc FCONE);
}
