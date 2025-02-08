#ifndef DPCA_H_
#define DPCA_H_

#include <Rinternals.h>

SEXP R_dpca(SEXP r_x, SEXP r_q, SEXP r_freqs, SEXP r_bandwidth, SEXP r_tol,
            SEXP kernel, SEXP r_max_q, SEXP r_select_q, SEXP r_n_path,
            SEXP r_t_path, SEXP r_penalties, SEXP r_penalty_scales);

#endif // DPCA_H_
