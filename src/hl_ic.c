#include "hl_ic.h"
#include <complex.h>
#include <math.h>

void hl_ic(double _Complex * spec, double _Complex * evals, int max_q, int nfreqs,
           int dim, int ldm, int select_q, double * ic_vals) {

    double total_trace = 0.0;

     for (int i=0; i < nfreqs; i++) {
        for (int j = 0; j < dim; j++) {
            total_trace += creal(spec[i * ldm * ldm + j * (ldm + 1)]);
        }
     }

     for (int i = 0; i <= max_q; i++) {

         if (i > 0){
             for (int j = 0; j < nfreqs; j++) {
                 total_trace -= creal(evals[j * max_q + (i - 1)]);
             }
         }
         ic_vals[i] = total_trace/((double) (dim * nfreqs));

         if (select_q == 2) {
             ic_vals[i] = log(ic_vals[i]);
         }
     }
}
