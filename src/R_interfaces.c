#include "R_interfaces.h"

#ifndef  USE_FC_LEN_T
# define USE_FC_LEN_T
#endif

#include <Rconfig.h>

#ifndef FCONE
# define FCONE
#endif

#include "complex_mv_product.h"
#include "complex_crossprod.h"
#include "eigs.h"
#include "filter_process.h"
#include "fourier_inverse.h"
#include "fourier_transform.h"
#include "hl_find_stability_interval.h"
#include "hl_ic.h"
#include "hl_ic_n_path.h"
#include "hl_q_path.h"
#include "hl_select_q.h"
#include "lagged_cov.h"

SEXP R_arnoldi_eigs(SEXP r_mat, SEXP r_dim, SEXP r_q, SEXP r_tol,
                    SEXP r_normalize_evecs, SEXP r_verbose,
                    SEXP r_row_evecs, SEXP r_transpose_out) {

  Rcomplex *mat = COMPLEX(r_mat);
  int ldm = nrows(r_mat);
  int dim = *INTEGER(r_dim);
  int q = *INTEGER(r_q);
  double tol = *REAL(r_tol);
  int normalize_evecs = *INTEGER(r_normalize_evecs);
  int verbose = *INTEGER(r_verbose);
  int row_evecs = *INTEGER(r_row_evecs);
  int transpose_out = *INTEGER(r_transpose_out);

  // result objects
  SEXP res = PROTECT(allocVector(VECSXP, 2));;
  SEXP evecs;
  if (transpose_out){
    evecs = PROTECT(allocMatrix(CPLXSXP, q, dim));
  } else {
    evecs = PROTECT(allocMatrix(CPLXSXP, dim, q));
  }
  SEXP evals = PROTECT(allocVector(CPLXSXP, q));

  arnoldi_eigs(
      (double _Complex *) mat, dim, ldm, q,
      (double _Complex *) COMPLEX(evals),
      (double _Complex *) COMPLEX(evecs), tol,
      normalize_evecs, verbose, row_evecs,
      transpose_out);

  SET_VECTOR_ELT(res, 0, evals);
  SET_VECTOR_ELT(res, 1, evecs);

  SEXP nms = PROTECT(allocVector(STRSXP, 2));
  SET_STRING_ELT(nms, 0, mkChar("values"));
  SET_STRING_ELT(nms, 1, mkChar("vectors"));

  setAttrib(res, R_NamesSymbol, nms);

  UNPROTECT(4);
  return res;
}


SEXP R_zMatVec(SEXP r_mat, SEXP r_vec, SEXP r_dim, SEXP version) {
  int dim = *INTEGER(r_dim);
  SEXP res = PROTECT(allocVector(CPLXSXP, dim));
  int ldm = nrows(r_mat);

  if (ldm < dim) {
    error("ldm must be larger than dim!");
  }

  complex_mv_product(
      (double _Complex *) COMPLEX(r_vec),
      (double _Complex *) COMPLEX(res),
      (double _Complex *) COMPLEX(r_mat), dim, ldm);

  UNPROTECT(1);
  return res;
}



SEXP R_complex_crossprod(SEXP r_x, SEXP r_trans_conj) {

    int nrx = nrows(r_x);
    int ncx = ncols(r_x);
    int trans_conj = *INTEGER(r_trans_conj);

    int dim;
    if (trans_conj) {
        dim = nrx;
    } else {
        dim = ncx;
    }

    SEXP res = PROTECT(allocMatrix(CPLXSXP, dim, dim));

    complex_crossprod((_Complex double *) COMPLEX(r_x), nrx, ncx, (_Complex double *) COMPLEX(res), trans_conj);

    UNPROTECT(1);
    return res;
}


SEXP R_filter_process(SEXP r_f, SEXP r_x, SEXP r_lags,
                      SEXP r_nrf, SEXP r_ncf, SEXP r_nrx, SEXP r_ncx,
                      SEXP r_nlags, SEXP r_inx, SEXP r_transf, SEXP r_rev) {

    int nrf = *INTEGER(r_nrf);
    int ncf = *INTEGER(r_ncf);
    int ncx = *INTEGER(r_ncx);
    int transf = *INTEGER(r_transf);
    int rev = *INTEGER(r_rev);
    int nry = nrf;
    if (transf) nry = ncf;
    SEXP res = PROTECT(allocMatrix(REALSXP, nry , ncx));

    filter_process(REAL(r_f), REAL(r_x), INTEGER(r_lags), nrf, *INTEGER(r_ncf),
                   *INTEGER(r_nrx), *INTEGER(r_ncx), *INTEGER(r_nlags), REAL(res),
                   *INTEGER(r_inx), transf, rev);
    UNPROTECT(1);
    return res;
}


SEXP R_fourier_inverse(SEXP r_f, SEXP r_nrf, SEXP r_ncf,
                     SEXP r_lags, SEXP r_nlags,
                     SEXP r_freqs, SEXP r_nfreqs) {

    int nrf = *INTEGER(r_nrf);
    int ncf = *INTEGER(r_ncf);
    int nlags = *INTEGER(r_nlags);
    int nfreqs = *INTEGER(r_nfreqs);
    double accum;


    SEXP res = PROTECT(alloc3DArray(REALSXP, nrf, ncf, nlags));

    fourier_inverse((double _Complex *) COMPLEX(r_f), nrf, ncf,
                    INTEGER(r_lags), nlags, REAL(r_freqs), nfreqs, REAL(res), &accum);



    if (accum > 1e-9) {
        Rprintf("Accum %f\n", accum);
        warning("The imaginary part of the coefficients was not zero, probably due to an assymmetric spectrum!");
    }

    UNPROTECT(1);
    return res;
}


SEXP R_fourier_transform(SEXP z, SEXP dim_z1, SEXP dim_z2,
                       SEXP freq, SEXP n_freq,
                         SEXP lags, SEXP
                         n_lags)
{


    SEXP res = PROTECT(alloc3DArray(CPLXSXP, *INTEGER(dim_z1), *INTEGER(dim_z2), *INTEGER(n_freq)));

    fourier_transform(REAL(z), *INTEGER(dim_z1), *INTEGER(dim_z2),
                       REAL(freq), *INTEGER(n_freq),
                       INTEGER(lags), *INTEGER(n_lags),
                       (double _Complex *) COMPLEX(res));

    UNPROTECT(1);
    return res;
}



SEXP R_find_stability_intervals(SEXP r_sample_var) {

    SEXP res = PROTECT(allocVector(INTSXP, 4));

    int info = hl_find_stability_intervals(REAL(r_sample_var), length(r_sample_var) , INTEGER(res));

    if (info == 1)
        warning("No stability inverval where S_c is zero was found. Returning indices of second stability interval where S_c is different from zero.");
    if (info == 2)
        warning("No stability inverval found. Using index of global minimum.");

    UNPROTECT(1);
    return res;
}


SEXP R_hl_ic(SEXP r_spec, SEXP r_evals, SEXP r_max_q, SEXP r_nfreqs, SEXP r_dim,
             SEXP r_ldm, SEXP r_select_q) {

    SEXP res = PROTECT(allocVector(REALSXP, *INTEGER(r_max_q) + 1));

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




SEXP R_hl_ic_n_path(SEXP r_spec, SEXP r_n_path, SEXP r_max_q, SEXP r_dim,
                   SEXP r_nfreqs, SEXP r_select_q, SEXP r_tol) {

    int nfreqs = *INTEGER(r_nfreqs);
    int max_q = *INTEGER(r_max_q);
    int ln = length(r_n_path);
    int dim = *INTEGER(r_dim);
    int select_q = *INTEGER(r_select_q);

    SEXP res = PROTECT(allocMatrix(REALSXP, max_q + 1, ln));
    double _Complex *evals, *evecs;
    evals = (double _Complex *) R_Calloc(nfreqs * max_q, double _Complex);
    evecs = (double _Complex *) R_Calloc(nfreqs * max_q * dim, double _Complex);

    hl_ic_n_path((double _Complex *) COMPLEX(r_spec), evals, evecs, dim, nfreqs, max_q, select_q,
                INTEGER(r_n_path), ln, *REAL(r_tol), REAL(res));

    R_Free(evals);
    R_Free(evecs);

    // INTEGER(res)[0] = 1;
    UNPROTECT(1);
    return res;
}


SEXP R_hl_q_path(SEXP r_unpenalized_ic_vals, SEXP r_max_q,
                 SEXP r_penalty_scale,
                 SEXP r_penalties) {
    int ln = ncols(r_unpenalized_ic_vals);
    SEXP res = PROTECT(allocVector(INTSXP, ln));

    hl_q_path(REAL(r_unpenalized_ic_vals), ln, *INTEGER(r_max_q),
              *REAL(r_penalty_scale), REAL(r_penalties), INTEGER(res));

    UNPROTECT(1);
    return res;
}



SEXP R_hl_select_q(SEXP r_spec, SEXP r_n_path, SEXP r_max_q, SEXP r_dim,
                   SEXP r_nfreqs, SEXP r_select_q, SEXP r_tol, SEXP r_penalties,
                   SEXP r_penalty_scales) {


    int max_q = *INTEGER(r_max_q);
    int dim = *INTEGER(r_dim);
    int nfreqs = *INTEGER(r_nfreqs);
    int ln = length(r_n_path);
    int lps = length(r_penalty_scales);

    SEXP evals = PROTECT(allocMatrix(CPLXSXP, max_q, nfreqs));
    SEXP evecs = PROTECT(alloc3DArray(CPLXSXP, max_q, dim, nfreqs));
    SEXP unpenalized_ic_vals = PROTECT(allocMatrix(REALSXP, max_q + 1, ln));
    SEXP sample_var = PROTECT(allocVector(REALSXP, lps));
    SEXP info = PROTECT(allocVector(INTSXP, 1));
    SEXP q = PROTECT(allocVector(INTSXP, 1));


    hl_select_q((_Complex double *) COMPLEX(r_spec),
                (_Complex double *) COMPLEX(evals),
                (_Complex double *) COMPLEX(evecs),
                dim,
                *INTEGER(r_nfreqs), max_q, *INTEGER(r_select_q),
                INTEGER(r_n_path), length(r_n_path),
                *REAL(r_tol),
                REAL(unpenalized_ic_vals),
                REAL(r_penalties),
                REAL(r_penalty_scales), length(r_penalty_scales),
                REAL(sample_var),
                INTEGER(info),
                INTEGER(q));


    /* for (int i=0; i < dim * max_q * nfreqs; i++) { */
    /*     printf("%f+%fi, ", COMPLEX(evecs)[i].r, COMPLEX(evecs)[i].i); */
    /* } */


    SEXP res = PROTECT(allocVector(VECSXP, 6));
    SET_VECTOR_ELT(res, 0, evals);
    SET_VECTOR_ELT(res, 1, evecs);
    SET_VECTOR_ELT(res, 2, unpenalized_ic_vals);
    SET_VECTOR_ELT(res, 3, sample_var);
    SET_VECTOR_ELT(res, 4, info);
    SET_VECTOR_ELT(res, 5, q);

    SEXP nms = PROTECT(allocVector(STRSXP, 6));
    SET_STRING_ELT(nms, 0, mkChar("evals"));
    SET_STRING_ELT(nms, 1, mkChar("evecs"));
    SET_STRING_ELT(nms, 2, mkChar("unpenalized_ic_vals"));
    SET_STRING_ELT(nms, 3, mkChar("sample_var_criterion"));
    SET_STRING_ELT(nms, 4, mkChar("info"));
    SET_STRING_ELT(nms, 5, mkChar("q"));
    setAttrib(res, R_NamesSymbol, nms);

    UNPROTECT(8);
    return res;
}



// here r_x and r_y are {n,m} by T arrays
SEXP R_lagged_cov(SEXP r_x, SEXP r_y, SEXP r_lag, SEXP r_nrx,
                  SEXP r_ncx, SEXP r_nry, SEXP r_ncy,
                  SEXP r_weight, SEXP r_center) {
    int ncx = *INTEGER(r_ncx);
    int nrx = *INTEGER(r_nrx);
    int ncy = *INTEGER(r_ncy);
    int nry = *INTEGER(r_nry);
    int lag = *INTEGER(r_lag);
    int center = *INTEGER(r_center);

    if (length(r_x) != nrx * ncx)
        error("Supplied dimension of X do not fit to supplied array length!");
    if (length(r_y) != nry * ncy)
        error("Supplied dimension of Y do not fit to supplied array length!");
    if (ncx != ncy)
        error("X and Y must have same number of columns!");
    if (ncx - 1 <= lag)
        error("Not enough observations to compute covariance with lag %d!", lag);

    SEXP res = PROTECT(allocMatrix(REALSXP, nrx, nry));

    lagged_cov(REAL(r_x), REAL(r_y), REAL(res), lag, nrx, ncx, nry, center, *REAL(r_weight));

    UNPROTECT(1);
    return res;
}


SEXP R_lagged_covs(SEXP r_x, SEXP r_y, SEXP r_lags, SEXP r_nrx, SEXP r_ncx,
                   SEXP r_nry, SEXP r_ncy, SEXP r_weights, SEXP r_center) {

    int ncx = *INTEGER(r_ncx);
    int nrx = *INTEGER(r_nrx);
    int ncy = *INTEGER(r_ncy);
    int nry = *INTEGER(r_nry);
    int center = *INTEGER(r_center);
    int *lags = INTEGER(r_lags);
    int nlags = length(r_lags);
    int maxlag = 0;

    for (int i = 0; i < nlags; i++)
        if (abs(lags[i]) > maxlag) maxlag = abs(lags[i]);

    if (length(r_x) != nrx * ncx)
        error("Supplied dimension of X do not fit to supplied array length!");
    if (length(r_y) != nry * ncy)
        error("Supplied dimension of Y do not fit to supplied array length!");
    if (ncx != ncy)
        error("X and Y must have same number of columns!");
    if (ncx - 1 <= maxlag)
        error("Not enough observations to compute covariance with maximum lag %d!", maxlag);

    SEXP res = PROTECT(alloc3DArray(REALSXP, nrx, nry, nlags));

    lagged_covs(REAL(r_x), REAL(r_y), REAL(res), INTEGER(r_lags), nlags, nrx, ncx, nry, ncy, REAL(r_weights), center);

    UNPROTECT(1);
    return res;
}



