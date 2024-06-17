#include "hl_ic_n_path.h"
#include "Rinternals.h"
#include "eigs.h"
#include "hl_ic.h"
#include "complex.h"
#include <R_ext/Utils.h>

void hl_ic_n_path(double _Complex * spec, double _Complex * evals, double _Complex *evecs,
                int dim, int nfreqs, int max_q, int select_q, int * n_path, int ln,
                double tol, double * ic_vals) {


    int curr_dim;

    /*Compute spectral decomposition for any n in n_path */
    // NOTE: This is only done for the first nfreq/2 slices of the spec array
    // since every spectrum has to be symmetric around zero. If the passed
    // spec array is non-symmetric, this will result in undefined behaviour.
    for (int i = 0; i < ln; i++)  {

        curr_dim = n_path[i];

        if (nfreqs % 2) {
            int j = nfreqs / 2;
            arnoldi_eigs(spec + dim * dim * j, curr_dim, dim, max_q,
                         evals + j * max_q, evecs + dim * max_q * j, tol,
                         1, 0, 1, 1);
        }
        for (int k = 0; k < nfreqs / 2; k++) {

            R_CheckUserInterrupt();

            int j = nfreqs / 2 + k + nfreqs % 2;
            int j_neg = nfreqs / 2 - k - 1;


            arnoldi_eigs(spec + dim * dim * j, curr_dim, dim, max_q,
                         evals + j * max_q, evecs + dim * max_q * j, tol,
                         1, 0, 1, 1);

            for (int l = 0; l < max_q; l++) {
                *(evals + max_q * j_neg + l) = *(evals + max_q * j + l);
                for (int m = 0; m < dim ; m++) {
                    *(evecs + max_q * dim * j_neg + dim * l + m) = conj(*(evecs + max_q * dim * j + dim * l + m));
                }
            }

        }

        hl_ic(spec, evals, max_q, nfreqs, curr_dim, dim, select_q, ic_vals + i * (max_q + 1));
    }
}
