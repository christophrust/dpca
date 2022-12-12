#include "R_ext/RS.h"
#include <complex.h>
#include "dpca.h"

int hl_select_q(double _Complex * spec, double _Complex * evals, double _Complex *evecs,
                int dim, int nfreqs, int max_q, int select_q, int * n_path, int ln, double tol) {

    double total_trace = 0.0;
    double ic_val[max_q * ln];

    double _Complex *temp_evecs, *temp_evals;
    temp_evecs = (double _Complex *) R_Calloc(dim * nfreqs * (max_q -1) * max_q/2, double _Complex);
    temp_evals = (double _Complex *) R_Calloc(nfreqs * (max_q -1) * max_q/2, double _Complex);
    int curr_dim;

    /*Compute spectral decomposition for any n in n_path */
    int curr_idx = 0;
    for (int i = 0; i < ln; i++)  {
        curr_dim = n_path[i];
        arnoldi_eigs((Rcomplex *) spec, curr_dim, dim, max_q,
                     (Rcomplex *) temp_evals, (Rcomplex *) temp_evecs, tol,
                     1, 0, 1, 1);
        curr_idx += i;
    }

    /* Compute (unpenalized) Information Criterion for every */
    for (int i=0; i < nfreqs; i++) {
        for (int j = 0; j < dim; j++) {
            total_trace += creal(spec[i * dim * dim + j * (dim + 1)]);
        }
    }

    for (int i = 0; i < max_q; i++) {
        ic_val[i] = total_trace;
        for (int j = 0; j < nfreqs; j++) {
            ic_val[i] -= creal(evals[i * nfreqs + j]);
        }
        ic_val[i] *= 1.0/((double) (dim * nfreqs));

        if (select_q == 2) {
            ic_val[i] = log(ic_val[i]);
        }
        ic_val[i] += i*1;
    }

    R_Free(temp_evecs);
    R_Free(temp_evals);
    return 1;
}
