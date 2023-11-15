#ifndef DPCA_H_
#define DPCA_H_

#include <R.h>
#include <Rinternals.h>
#include <complex.h>

#include "complex_crossprod.h"
#include "eigs.h"
#include "filter_process.h"
#include "fourier_inverse.h"
#include "fourier_transform.h"
#include "complex_mv_product.h"
#include "hl_find_stability_interval.h"
#include "hl_ic.h"
#include "hl_ic_n_path.h"
#include "hl_q_path.h"
#include "hl_select_q.h"
#include "lagged_cov.h"


SEXP R_dpca(SEXP r_x, SEXP r_q, SEXP r_freqs, SEXP r_bandwidth,
            SEXP r_tol, SEXP kernel, SEXP r_max_q, SEXP r_select_q,
            SEXP r_n_path, SEXP r_t_path, SEXP r_penalties, SEXP r_penalty_scales);




#endif // DPCA_H_
