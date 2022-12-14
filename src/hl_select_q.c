#include "Rinternals.h"
#include "dpca.h"



void hl_select_q(double _Complex * spec, double _Complex * evals, double _Complex *evecs,
                int dim, int nfreqs, int max_q, int select_q, int * n_path, int ln,
                double tol, double * unpenalized_ic_vals, double * penalties,
                double * penalty_scales, int lps, double * sample_var) {


    /* 1. obtain unpenalized ic vals for all q and n in n_path */
    hl_ic_n_path(spec, evals, evecs,
                 dim, nfreqs, max_q, select_q,  n_path, ln, tol,
                 unpenalized_ic_vals);


    /* 2. obtain for any c in the penalty_scales the q_path and compute its sample variability */
    int q_paths[ln * lps];
    double penalty_scale;
    // double sample_var[lps];
    double x2, x1;
    for (int i = 0; i < lps; i++) {
        penalty_scale = penalty_scales[i];
        hl_q_path(unpenalized_ic_vals, ln, max_q, penalty_scale, penalties, q_paths + i * ln);

        x2 = 0.0; x1 = 0.0;
        for (int j = 0; j<ln; j++) {
            x2 += (q_paths[i * ln + j] * q_paths[i * ln + j]);
            x1 += q_paths[i * ln + j];
        }
        sample_var[i] = x2 / ln - (x1 * x1)/ln;
    }

    /* 3.
     * find second 'stability interval' of sample variability and choose q either using q^T_c,n  or
     * by majority vote over the corresponding q_path
     */

}



SEXP R_hl_select_q(SEXP r_spec, SEXP r_n_path, SEXP r_max_q, SEXP r_dim,
                   SEXP r_nfreqs, SEXP r_select_q, SEXP r_tol) {

    SEXP res = PROTECT(allocVector(INTSXP, 1));

    INTEGER(res)[0] = 0;
    return res;
}
