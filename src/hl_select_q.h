#ifndef HL_SELECT_Q_H_
#define HL_SELECT_Q_H_


/**
 * @param spec Estimated spectrum given as n by n by nfreqs array
 * @param evals An max_q by nfreqs array which will (on output) hold in
 * the first q by nfreqs entries the resulting eigenvalues.
 * @param evecs An n by max_q by nfreqs array which will (on output) hold
 * in the first n by q by nfreqs entries the resulting eigenvectors.
 * @param dim Dimension of the spectral density matrix.
 * @param nfreqs Number of frequencies where the spectrum is evaluated.
 * @param max_q The maximum number of dynamic factors.
 * @param select_q At this stage one of 1 or 2. 1 indicates that the IC1
 * criterion is used to select q, and 2 indicates IC2.
 * @param n_path Array holding the different values for subspectra to compute
 * the information criteria on to do fine tuning.
 * @param ln The length of n_path.
 * @param tol Tolerance used in the ARPACK routine
 * @param unpenalized_ic_vals An array of at least length ln * max_q containing
 * (on output) the unpenalized values of the selected information criterium.
 * @param penalties An array of length ln giving the penalty for every
 * element in n_path.
 * @param penalty_scales An array of length lps containing different scaling values
 * of the penalty.
 * @param lps The lenght of penalty_scale.
 * @param sample_var An array of length lps which (on output) holds the sample
 * variance of the q_choice over all entries in n_path for each element of penalty_scales.
 * @param info Single integer. If info = 0 then everything went fine. If info = 1,
 * no zero stability invervals were found. If info = 2, no stability was found, such that
 * the penalty scale which globally minimizes the sample variance is chosen.
 * @param q The number factors (on output).
 * @param q_path an integer array of length lps which will contain the q_path
 * over all values of penalty_scales for n.
 * */
void hl_select_q
(
    double _Complex * spec,
    double _Complex * evals,
    double _Complex *evecs,
    int dim,
    int nfreqs,
    int max_q,
    int select_q,
    int * n_path,
    int ln,
    double tol,
    double * unpenalized_ic_vals,
    double * penalties,
    double * penalty_scales,
    int lps,
    double * sample_var,
    int *info,
    int *q,
    int* q_path
);


#endif // HL_SELECT_Q_H_
