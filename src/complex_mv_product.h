#ifndef COMPLEX_MV_PRODUCT_H_
#define COMPLEX_MV_PRODUCT_H_

#include <R.h>
#include <Rinternals.h>

void zMatVecLa(double _Complex *x, double _Complex* y,
               Rcomplex* mat, int dim, int ldm);


void zMatVec(double _Complex *x, double _Complex* y,
             Rcomplex* mat, int dim);


#endif // COMPLEX_MV_PRODUCT_H_
