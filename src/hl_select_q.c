#include "Rinternals.h"
#include "dpca.h"


/**
 * @brief Implementation of the the Hallin & Liska (2007) procedure to determine
 * the numbers of dynamic factors for a spectrum.
 *
 * @param spec Estimated spectrum given as n by n by nfreqs array
 * @param evals An max_q by nfreqs array which will (on output) hold in
 * the first q by nfreqs entries the resulting eigenvalues.
 * @param evecs An n by max_q by nfreqs array which will (on output) hold
 * in the first n by q by nfreqs entries the resulting eigenvectors.
 * @param dim Dimension of the spectral density matrix.
 * @param nfreqs Number of frequencies where the spectrum is evaluated.
 * @param max_q The maximum number of dynamic factors.
 * @param select_q At this stage one of 1 or 2. 1 indicates that the IC1
 * criterion is used to select q, and 2 indicates IC2.
 * @param n_path Array holding the different values for subspectra to compute
 * the information criteria on to do fine tuning.
 * @param ln The length of n_path.
 * @param tol Tolerance used in the ARPACK routine
 * @param unpenalized_ic_vals An array of at least length ln * max_q containing
 * (on output) the unpenalized values of the selected information criterium.
 * @param penalties An array of length ln giving the penalty for every
 * element in n_path.
 * @param penalty_scales An array of length lps containing different scaling values
 * of the penalty.
 * @param lps The lenght of penalty_scale.
 * @param sample_var An array of length lps which (on output) holds the sample
 * variance of the q_choice over all entries in n_path for each element of penalty_scales.
 * @param info Single integer. If info = 0 then everything went fine. If info = 1,
 * no zero stability invervals were found. If info = 2, no stability was found, such that
 * the penalty scale which globally minimizes the sample variance is chosen.
 * @param q The number factors (on output).
 * */
void hl_select_q(double _Complex * spec, double _Complex * evals, double _Complex *evecs,
                 int dim, int nfreqs, int max_q, int select_q, int * n_path, int ln,
                 double tol, double * unpenalized_ic_vals, double * penalties,
                 double * penalty_scales, int lps, double * sample_var, int *info, int *q) {


    /* 1. obtain unpenalized ic vals for all q and n in n_path */
    hl_ic_n_path(spec, evals, evecs,
                 dim, nfreqs, max_q, select_q,  n_path, ln, tol,
                 unpenalized_ic_vals);


    /* 2. obtain for any c in the penalty_scales the q_path and compute its sample variability */
    int q_paths[ln * lps];
    double penalty_scale;

    double x2, x1;
    for (int i = 0; i < lps; i++) {
        penalty_scale = penalty_scales[i];
        hl_q_path(unpenalized_ic_vals, ln, max_q, penalty_scale, penalties, q_paths + i * ln);

        x2 = 0.0; x1 = 0.0;
        for (int j = 0; j<ln; j++) {
            x2 += (q_paths[i * ln + j] * q_paths[i * ln + j]);
            x1 += q_paths[i * ln + j];
        }

        sample_var[i] = x2 / ln - (x1 * x1)/(ln * ln);

        // printf("sample_var[%i]: %f; x2: %f, x1: %f\n", i, sample_var[i], x2, x1);
    }

    /* 3.
     * find second 'stability interval' of sample variability and choose q either using q^T_c,n  or
     * by majority vote over the corresponding q_path
     */
    int stability_interval[4];
    *info = hl_find_stability_intervals(sample_var, lps, stability_interval);

    *q = q_paths[stability_interval[2] * ln - 1];
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
    SEXP evecs = PROTECT(alloc3DArray(CPLXSXP, dim, max_q, nfreqs));
    SEXP unpenalized_ic_vals = PROTECT(allocMatrix(REALSXP, max_q, ln));
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
                REAL(sample_var), INTEGER(info), INTEGER(q));


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



