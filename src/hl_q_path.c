#include "R_ext/RS.h"
#include "dpca.h"

// TODO: start with q = 0!!
void hl_q_path(double *unpenalized_ic_vals, int ln, int max_q,
                      double penalty_scale,
                      double *penalties,
                      int * q_path) {

    int q;
    double ic, ic_temp;

    for (int i = 0; i < ln; i++) {
        q = 0;
        ic = unpenalized_ic_vals[i * max_q] + penalties[i] * penalty_scale;
        for (int j = 1; j < max_q; j++) {
            ic_temp = unpenalized_ic_vals[i * max_q + j] +
                ((double) j + 1.0) * penalties[i] * penalty_scale;
            if (ic_temp < ic) {
                ic = ic_temp;
                q = j + 1;
            }
        }
        q_path[i] = q;
    }
}
