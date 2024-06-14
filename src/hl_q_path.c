#include "hl_q_path.h"

void hl_q_path(double *unpenalized_ic_vals, int ln, int max_q,
                      double penalty_scale,
                      double *penalties,
                      int * q_path) {

    int q;
    double ic, ic_temp;

    for (int i = 0; i < ln; i++) {
        q = 0;
        ic = unpenalized_ic_vals[i * (max_q + 1)];
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
