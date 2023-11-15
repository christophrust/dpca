#ifndef EIGS_H_
#define EIGS_H_

#include <Rinternals.h>


/**
 * @brief Compute truncated eigendecomposition of a hermitian matrix.
 *
 * @param mat An array of at least length ldm \times ldm holding the
 *   matrix to be decomposed (INPUT).
 * @param dim The dimension of the matrix mat.
 * @param ldm Leading dimension of mat.
 * @param q Number of eigenpairs to be computed.
 * @param eval Array of at least lenght q which on output will contain the
 *   computed eigenvalues.
 * @param evecs Array of at least length dim times q which will hold the
 *   eigenvectors on output.
 * @param tol Tolerance passed to the ARPACK routine.
 * @param normalize_evecs If set to a value different from zero, each of the resulting
 *   eigenvectors is normalized such that the complex part of the first entry is zero
 *   (the same is done in R's base::eigen()).
 * @param verbose If different from zero, we show some messages.
 * @param row_evecs If different from zero, the eigenvalues are roweigenvectors instead
 *   of column-eigenvectors
 * @param transpose_out If different from zero, the eigenvectors leading
 *  dimension is q, otherwise ldm
 */
void arnoldi_eigs(Rcomplex *mat, int dim, int ldm, int q,
                  Rcomplex *eval, Rcomplex *evecs,
                  double tol, int normalize_evecs, int verbose,
                  int row_evecs, int transpose_out);

#endif // EIGS_H_
