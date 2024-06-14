#include <float.h>
#include <math.h>

#include "hl_find_stability_interval.h"

int hl_find_stability_intervals(double * sample_var, int lsv, int *ivl_idx) {

    ivl_idx[0] = 0;
    ivl_idx[1] = 0;

    int itvl_cnt = 0;
    int in_itvl = 1;
    int info = 2;
    for (int i = 1; i < lsv; i++) {

        if ((fabs(sample_var[i]) < 10 * DBL_EPSILON) && (fabs(sample_var[i-1]) < 10 * DBL_EPSILON)) {
            in_itvl = 1;
            ivl_idx[itvl_cnt * 2 + 1] = i;

        } else {
            if (in_itvl) itvl_cnt++;
            in_itvl = 0;
            if (itvl_cnt <= 1)
                ivl_idx[itvl_cnt * 2] = i;
        }

        if (itvl_cnt > 1) {
            info=0;
            break;
        } else if ((i + 1) == lsv && itvl_cnt == 1 && fabs(sample_var[i]) < 10 * DBL_EPSILON) {
            info = 0;
        }

    }

    /* in case we have not succeeded so far, we look stability
     * intervals where S_c is different from zero but still locally minimal */

    if (info) {
        itvl_cnt = 0;
        in_itvl = 1;
        for (int i = 1; i < lsv; i++) {

            if ((fabs(sample_var[i] - sample_var[i-1]) < 10 * DBL_EPSILON)) {
                in_itvl = 1;
                ivl_idx[itvl_cnt * 2 + 1] = i;

            } else {
                if (in_itvl &&
                    ((itvl_cnt == 0 && sample_var[i] > sample_var[i-1]) ||
                     (itvl_cnt == 1 && (sample_var[i] > sample_var[i-1]) &&
                      sample_var[ivl_idx[itvl_cnt * 2]] < sample_var[ivl_idx[itvl_cnt * 2 ] - 1]))) {
                    itvl_cnt++;
                }
                in_itvl = 0;
                if (itvl_cnt <= 1)
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

        double temp_min = sample_var[0];
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
