#ifndef COMPLEX_MV_PRODUCT_H_
#define COMPLEX_MV_PRODUCT_H_

#include <R.h>
#include <Rinternals.h>

void complex_mv_product(double _Complex *x, double _Complex* y,
               Rcomplex* mat, int dim, int ldm);

#endif // COMPLEX_MV_PRODUCT_H_
