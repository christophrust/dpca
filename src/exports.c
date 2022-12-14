#include "dpca.h"


static const R_CallMethodDef CallEntries[] = {
  {"R_arnoldi_eigs", (DL_FUNC) &R_arnoldi_eigs, 8},
  {"R_zMatVec", (DL_FUNC) &R_zMatVec, 4},
  {"R_lagged_cov", (DL_FUNC) &R_lagged_cov, 9},
  {"R_lagged_covs", (DL_FUNC) &R_lagged_covs, 9},
  {"R_fourier_transform", (DL_FUNC) &R_fourier_transform, 7},
  {"R_fourier_inverse1", (DL_FUNC) &R_fourier_inverse1, 7},
  {"R_fourier_inverse", (DL_FUNC) &R_fourier_inverse, 7},
  {"R_filter_process", (DL_FUNC) &R_filter_process, 11},
  {"R_dpca", (DL_FUNC) &R_dpca, 9},
  {"R_complex_crossprod", (DL_FUNC) &R_complex_crossprod, 2},
  {"R_recursive_filter", (DL_FUNC) &R_recursive_filter, 4},
  {"R_hl_ic", (DL_FUNC) &R_hl_ic, 7},
  {"R_hl_ic_n_path", (DL_FUNC) &R_hl_ic_n_path, 7},
  {"R_hl_q_path", (DL_FUNC) &R_hl_q_path, 4},
  {NULL, NULL, 0}
};



void R_init_dpca(DllInfo *dll) {
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
