#include "dpca.h"


void hl_ic(double _Complex * spec, double _Complex * evals, int max_q, int nfreqs,
           int dim, int select_q, double * ic_vals) {

    double total_trace = 0.0;

     for (int i=0; i < nfreqs; i++) {
        for (int j = 0; j < dim; j++) {
            total_trace += creal(spec[i * dim * dim + j * (dim + 1)]);
        }
     }

     for (int i = 0; i < max_q; i++) {
        ic_vals[i] = total_trace;
        for (int j = 0; j < nfreqs; j++) {
            ic_vals[i] -= creal(evals[i * nfreqs + j]);
        }
        ic_vals[i] *= 1.0/((double) (dim * nfreqs));

        if (select_q == 2) {
            ic_vals[i] = log(ic_vals[i]);
        }
        ic_vals[i] += i*1;
     }

}
