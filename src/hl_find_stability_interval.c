#include "dpca.h"
#include <float.h>
#include <math.h>

/**
 * @brief find the stability interval of the sample variance of S_c
 * over the range of c
 * @param sample_var Array of length lsv with the sample variance of S_c
 * @param lsv the length of sample_var
 * @param ivl_idx On output this holds the indices (start, end) of the
 * first two stability intervals.
 * */
int hl_find_stability_intervals(double * sample_var, int lsv, int *ivl_idx) {

    ivl_idx[0] = 0;
    ivl_idx[1] = 0;

    int itvl_cnt = 0;
    int in_itvl = 1;
    int info = 2;
    for (int i = 1; i < lsv; i++) {

        if ((fabs(sample_var[i]) < 10* DBL_EPSILON) && (fabs(sample_var[i-1]) < 10 * DBL_EPSILON)) {
            in_itvl = 1;
            ivl_idx[itvl_cnt * 2 + 1] = i;

        } else {
            if (in_itvl) itvl_cnt++;
            in_itvl = 0;
            ivl_idx[itvl_cnt * 2] = i;
        }

        if (itvl_cnt > 1) {
            info=0;
            break;
        }
    }

    /* in case we have not succeed so far, we look stability
     * intervals where S_c is different from zero but still locally minimal (TODO) */

    if (info) {
        itvl_cnt = 0;
        in_itvl = 1;
        for (int i = 1; i < lsv; i++) {

            if ((fabs(sample_var[i] - sample_var[i-1]) < 10 * DBL_EPSILON)) {
                in_itvl = 1;
                ivl_idx[itvl_cnt * 2 + 1] = i;

            } else {
                if (in_itvl) itvl_cnt++;
                in_itvl = 0;
                ivl_idx[itvl_cnt * 2] = i;
            }

            if (itvl_cnt > 1) {
                info = 1;
                break;
            }
        }
    }

    if (info == 2) {
        ivl_idx[0] = 0;
        ivl_idx[1] = 0;
        ivl_idx[2] = 0;
        ivl_idx[3] = 0;

        double temp_min= sample_var[0];
        for (int i = 1; i < lsv; i++) {
            if (sample_var[i] < temp_min) {
                temp_min = sample_var[i];
                ivl_idx[0] = i;
            }
        }
    }
    /* int i = 1; */
    /* int first_min = 1; */
    /* while((sample_var[i] < sample_var[i - 1] || first_min)) { */
    /*     i++; */
    /* } */

    return info;
}



SEXP R_find_stability_intervals(SEXP r_sample_var) {
    SEXP res = PROTECT(allocVector(INTSXP, 4));

    int info = hl_find_stability_intervals(REAL(r_sample_var), length(r_sample_var) , INTEGER(res));

    if (info == 1)
        warning("No stability inverval where S_c is zero was found. Returning indices of second stability interval where S_c is different from zero.");
    if (info == 2)
        warning("No stability inverval found. Using index of global minimum.");
    UNPROTECT(1);
    return res;
}
