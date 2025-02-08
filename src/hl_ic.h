#ifndef HL_IC_H_
#define HL_IC_H_

/**
 * @brief Compute the unpenalized information criterium from Hallin & Liska
 * (2007)
 *
 * @param spec A dim by dim by nfreqs array. Note that only the
 *   first subarrays 0 to floor(nfreqs/2) fo spec are used for the
 *   computation, assuming symmetric frequencies around zero.
 * @param evals The first max_q eigenvalues of spec.
 * @param max_q Maximum number of factors to be considered.
 * @param nfreqs Number of frequencies along the spectrum
 * @param dim Dimension of spectrum.
 * @param ldm Leading dimension of spec.
 * @param select_q Which information criterion (IC1 or IC2).
 * @param ic_vals (On output) the computed criteria (must be at least of length
 * max_q + 1).
 * */
void hl_ic(double _Complex *spec, double _Complex *evals, int max_q, int nfreqs,
           int dim, int ldm, int select_q, double *ic_vals);

#endif // HL_IC_H_
