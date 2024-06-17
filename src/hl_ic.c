#include "hl_ic.h"
#include <complex.h>
#include <math.h>

void hl_ic(double _Complex * spec, double _Complex * evals, int max_q, int nfreqs,
           int dim, int ldm, int select_q, double * ic_vals) {

    // NOTE: we only use the first nfreqs/2 slices of the spectrum array
    // to compute the Hallin & Liska criterion. If the passed spec array is
    // non-symmetric, this will result in undefined behaviour.
    double total_trace = 0.0;

    if (nfreqs % 2) {
        int i = nfreqs / 2;
        for (int j = 0; j < dim; j++) {
            total_trace += creal(spec[i * ldm * ldm + j * (ldm + 1)]);
        }
    }
    for (int i=0; i < nfreqs / 2; i++) {
        for (int j = 0; j < dim; j++) {
            total_trace += 2.0 * creal(spec[i * ldm * ldm + j * (ldm + 1)]);
        }
    }

    for (int i = 0; i <= max_q; i++) {

         if (i > 0){
             if (nfreqs % 2) {
                 int j = nfreqs / 2;
                 total_trace -= creal(evals[j * max_q + (i - 1)]);
             }
             for (int j = 0; j < nfreqs/2; j++) {
                 total_trace -= 2.0 * creal(evals[j * max_q + (i - 1)]);
             }
         }
         ic_vals[i] = total_trace/((double) (dim * nfreqs));

         if (select_q == 2) {
             ic_vals[i] = log(ic_vals[i]);
         }
     }
}
