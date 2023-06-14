#include "R_ext/RS.h"
#include <complex.h>
#include "Rinternals.h"
#include "dpca.h"


void hl_ic_n_path(double _Complex * spec, double _Complex * evals, double _Complex *evecs,
                int dim, int nfreqs, int max_q, int select_q, int * n_path, int ln,
                double tol, double * ic_vals) {

    // double ic_vals[max_q * ln];

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

    /* for (int i = 0; i < max_q * ln; i++) { */
    /*     printf("IC[%i]: %f\n", i, ic_vals[i]); */
    /* } */

}


SEXP R_hl_ic_n_path(SEXP r_spec, SEXP r_n_path, SEXP r_max_q, SEXP r_dim,
                   SEXP r_nfreqs, SEXP r_select_q, SEXP r_tol) {

    int nfreqs = *INTEGER(r_nfreqs);
    int max_q = *INTEGER(r_max_q);
    int ln = length(r_n_path);
    int dim = *INTEGER(r_dim);
    int select_q = *INTEGER(r_select_q);

    SEXP res = PROTECT(allocMatrix(REALSXP, max_q + 1, ln));
    double _Complex *evals, *evecs;
    evals = (double _Complex *) R_Calloc(nfreqs * max_q, double _Complex);
    evecs = (double _Complex *) R_Calloc(nfreqs * max_q * dim, double _Complex);

    hl_ic_n_path((double _Complex *) COMPLEX(r_spec), evals, evecs, dim, nfreqs, max_q, select_q,
                INTEGER(r_n_path), ln, *REAL(r_tol), REAL(res));

    R_Free(evals);
    R_Free(evecs);

    // INTEGER(res)[0] = 1;
    UNPROTECT(1);
    return res;
}
