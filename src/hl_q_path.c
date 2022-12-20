#include "R_ext/RS.h"
#include "Rinternals.h"
#include "dpca.h"

void hl_q_path(double *unpenalized_ic_vals, int ln, int max_q,
                      double penalty_scale,
                      double *penalties,
                      int * q_path) {

    int q;
    double ic, ic_temp;

    for (int i = 0; i < ln; i++) {
        q = 0;
        ic = unpenalized_ic_vals[i * (max_q + 1)] + penalties[i] * penalty_scale;
        for (int j = 1; j <= max_q; j++) {
            ic_temp = unpenalized_ic_vals[i * (max_q + 1) + j] +
                ((double) j) * penalties[i] * penalty_scale;
            if (ic_temp < ic) {
                ic = ic_temp;
                q = j;
            }
        }
        q_path[i] = q;
    }
}



SEXP R_hl_q_path(SEXP r_unpenalized_ic_vals, SEXP r_max_q,
                 SEXP r_penalty_scale,
                 SEXP r_penalties) {
    int ln = ncols(r_unpenalized_ic_vals);
    SEXP res = PROTECT(allocVector(INTSXP, ln));

    hl_q_path(REAL(r_unpenalized_ic_vals), ln, *INTEGER(r_max_q),
              *REAL(r_penalty_scale), REAL(r_penalties), INTEGER(res));

    UNPROTECT(1);
    return res;
}
