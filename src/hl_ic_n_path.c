#include "hl_ic_n_path.h"
#include "Rinternals.h"
#include "eigs.h"
#include "hl_ic.h"



void hl_ic_n_path(double _Complex * spec, double _Complex * evals, double _Complex *evecs,
                int dim, int nfreqs, int max_q, int select_q, int * n_path, int ln,
                double tol, double * ic_vals) {


    int curr_dim;

    /*Compute spectral decomposition for any n in n_path */
    for (int i = 0; i < ln; i++)  {

        curr_dim = n_path[i];

        for (int j = 0; j < nfreqs; j++) {
            arnoldi_eigs((Rcomplex *) spec + dim * dim * j, curr_dim, dim, max_q,
                         (Rcomplex *) evals + j * max_q, (Rcomplex *) evecs + dim * max_q * j, tol,
                         1, 0, 1, 1);
        }

        hl_ic(spec, evals, max_q, nfreqs, curr_dim, dim, select_q, ic_vals + i * (max_q + 1));
    }
}
