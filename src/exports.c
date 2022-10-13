#include "dpca.h"


static const R_CallMethodDef CallEntries[] = {
  {"R_arnoldi_eigs", (DL_FUNC) &R_arnoldi_eigs, 4},
  {"R_zMatVec", (DL_FUNC) &R_zMatVec, 4},
  {NULL, NULL, 0}
};



void R_init_dpca(DllInfo *dll) {
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
