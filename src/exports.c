#include "dpca.h"


static const R_CallMethodDef CallEntries[] = {
    {"R_arnoldi_eigs", (DL_FUNC) &R_arnoldi_eigs, 3},
    {NULL, NULL, 0}
};



void R_init_dpca(DllInfo *dll) {
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
