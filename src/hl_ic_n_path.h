#ifndef HL_IC_N_PATH_H_
#define HL_IC_N_PATH_H_


/**
 * @brief Selection of number of dynamic factors.
 *
 * @param spec The spectral density matrix (of dimension n by n by nfreqs).
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
 * @param ic_vals An array of at least length ln * (max_q + 1)
 *
 * */
void hl_ic_n_path(double _Complex * spec, double _Complex * evals, double _Complex *evecs,
                int dim, int nfreqs, int max_q, int select_q, int * n_path, int ln, double tol,
                double * ic_vals);


#endif // HL_IC_N_PATH_H_
