#include "R_ext/RS.h"
#include <complex.h>
#include "dpca.h"

int hl_select_q(double _Complex * spec, double _Complex * evals, int dim, int nfreqs, int max_q, int select_q){

    double total_trace = 0.0;
    double ic_val[max_q];

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

    return 1;
}
