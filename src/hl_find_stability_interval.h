#ifndef HL_FIND_STABILITY_INTERVAL_H_
#define HL_FIND_STABILITY_INTERVAL_H_

/**
 * @brief find the stability interval of the sample variance of S_c
 * over the range of c
 * @param sample_var Array of length lsv with the sample variance of S_c
 * @param lsv the length of sample_var
 * @param ivl_idx Pointer to an integer array of length 4.
 * On output this holds the indices (start, end) of the
 * first two stability intervals.
 * */
int hl_find_stability_intervals(double * sample_var, int lsv, int *ivl_idx);



#endif // HL_FIND_STABILITY_INTERVAL_H_
