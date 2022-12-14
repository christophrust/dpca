#include "Rinternals.h"
#include "dpca.h"



int hl_select_q(double _Complex * spec, double _Complex * evals, double _Complex *evecs,
                int dim, int nfreqs, int max_q, int select_q, int * n_path, int ln,
                double tol, double * unpenalized_ic_vals, double * penalties,
                double * penalty_scales, int lps) {

    /* 1. obtain unpenalized ic vals for all q and n in n_path */
    hl_ic_n_path(spec, evals, evecs,
                 dim, nfreqs, max_q, select_q,  n_path, ln, tol,
                 unpenalized_ic_vals);

    /* 2. obtain for any c in the penalty_scales the q-path and compute its sample variability */
    int q_paths[ln * lps];
    double penalty_scale;
    double sample_var[lps];

    for (int i = 0; i < lps; i++) {
        penalty_scale = penalty_scales[i];
        hl_q_path(unpenalized_ic_vals, ln, max_q, penalty_scale, penalties, q_paths + i * ln);

    }

    /* 3.
     * find first minimum of sample variability  and choose q by majority vote over
     * the corresponding q_path
     */
    return 1;
}



SEXP R_hl_select_q(SEXP r_spec, SEXP r_n_path, SEXP r_max_q, SEXP r_dim,
                   SEXP r_nfreqs, SEXP r_select_q, SEXP r_tol) {

    SEXP res = PROTECT(allocVector(INTSXP, 1));

    INTEGER(res)[0] = 0;
    return res;
}
