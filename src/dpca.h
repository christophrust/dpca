#ifndef DPCA_H_
#define DPCA_H_

#include <R.h>
#include <Rinternals.h>
#include <complex.h>
#include <arpack.h>
#include <stat_c.h>
#include <debug_c.h>


SEXP R_arnoldi_eigs(SEXP r_mat, SEXP r_dim, SEXP r_q);
#endif // DPCA_H_