#include "Rinternals.h"
#include "dpca.h"

// TODO: start with q = 0!!
void hl_ic(double _Complex * spec, double _Complex * evals, int max_q, int nfreqs,
           int dim, int ldm, int select_q, double * ic_vals) {

    double total_trace = 0.0;

     for (int i=0; i < nfreqs; i++) {
        for (int j = 0; j < dim; j++) {
            total_trace += creal(spec[i * ldm * ldm + j * (ldm + 1)]);
        }
     }

     for (int i = 0; i < max_q; i++) {
         for (int j = 0; j < nfreqs; j++) {
             total_trace -= creal(evals[j * max_q + i]);
         }
         ic_vals[i] = total_trace/((double) (dim * nfreqs));

         if (select_q == 2) {
             ic_vals[i] = log(ic_vals[i]);
         }
     }
}





SEXP R_hl_ic(SEXP r_spec, SEXP r_evals, SEXP r_max_q, SEXP r_nfreqs, SEXP r_dim,
             SEXP r_ldm, SEXP r_select_q) {

    SEXP res = PROTECT(allocVector(REALSXP, *INTEGER(r_max_q)));

    hl_ic((double _Complex *) COMPLEX(r_spec),
          (double _Complex *) COMPLEX(r_evals),
          *INTEGER(r_max_q),
          *INTEGER(r_nfreqs),
          *INTEGER(r_dim),
          *INTEGER(r_ldm),
          *INTEGER(r_select_q),
          REAL(res));

    UNPROTECT(1);
    return res;
}
