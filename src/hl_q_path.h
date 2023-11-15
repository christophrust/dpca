#ifndef HL_Q_PATH_H_
#define HL_Q_PATH_H_


/**
 * @brief Compute q-path for unpenalized values of information criterium and
 * penalty information.
 *
 * @param unpenalized_ic_vals A max_q by ln array with the unpenalized information
 * criterium.
 * @param ln Number of columns of unpenalized_ic_vals (= size of n path).
 * @param max_q Number of rows of unpenalized_ic_vals (= maximum order of the factor space).
 * @param penalty_scale The scale of the penalty.
 * @param penalties An array of length ln holding the penalties for each n in n_path.
 * @param q_path Array of at least lenght ln holding the chosen q for each n in n_path.
 *
 * */
void hl_q_path
(
    double *unpenalized_ic_vals,
    int ln,
    int max_q,
    double penalty_scale,
    double *penalties,
    int * q_path
);

#endif // HL_Q_PATH_H_
