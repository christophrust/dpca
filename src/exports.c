#include "dpca.h"


static const R_CallMethodDef CallEntries[] = {
  {"R_arnoldi_eigs", (DL_FUNC) &R_arnoldi_eigs, 6},
  {"R_zMatVec", (DL_FUNC) &R_zMatVec, 4},
  {"R_lagged_cov", (DL_FUNC) &R_lagged_cov, 7},
  {"R_lagged_covs", (DL_FUNC) &R_lagged_covs, 7},
  {"R_fourier_transform", (DL_FUNC) &R_fourier_transform, 7},
  {"R_fourier_inverse1", (DL_FUNC) &R_fourier_inverse1, 7},
  {"R_fourier_inverse", (DL_FUNC) &R_fourier_inverse, 7},

  {NULL, NULL, 0}
};



void R_init_dpca(DllInfo *dll) {
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
