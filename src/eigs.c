
#include <R.h>
#include <Rinternals.h>
#include <complex.h>
#include <arpack.h>
#include <stat_c.h>
#include <debug_c.h>
#include <R_ext/Lapack.h>
#include <math.h>

#include "dpca.h"


void get_rank(double *values, int *rank, int n) {

  for (int i = 0; i < n; i++) {
    int curRank=0;
    for (int j = 0; j < i; j++) {
      if (values[i] < values[j]) {
        curRank++;
      } else rank[j]++;
    }
    rank[i] = curRank;
  }
}


void arnoldi_eigs(Rcomplex *mat, int dim, int ldm, int q,
                  Rcomplex *eval, Rcomplex *evecs, double tol,
                  int normalize_evecs, int verbose, int row_evecs,
                  int transpose_out) {


  // znaupd parameters
  a_int ido = 0;                // reverse communication flag, (handled internally) must be zero at start
  char bmat[] = "I";            // B matrix ( I -> standard EV problem )
  a_int N = (a_int) dim;        // number of rows (dimension of the eigenproblem)
  char which[] = "LM";          // LM -> largest eigenvalues are of interest
  a_int nev = (a_int) q;        // Number of eigenvalues
  double _Complex resid[N];     // residual vector

  a_int ncv = 2 * nev + 1; //
  if (ncv < 20) ncv = 20;  // usage consistent to octave
  if (ncv > N) ncv = N;

  double _Complex V[ncv * N];
  a_int ldv = N;
  a_int iparam[11] = {0};
  a_int ipntr[14];
  for (int i=0; i< 14; i++) ipntr[i] = 0;

  double _Complex workd[3 * N];
  a_int rvec = 1;
  char howmny[] = "A";
  double _Complex d[nev+1];

  a_int select[ncv];
  for (int i = 0; i < ncv; i++) select[i] = 1;

  double _Complex z[N  * nev];
  a_int ldz = N;
  double _Complex sigma = 0. + I * 0.;
  int k;
  for (k = 0; k < 3 * N; ++k) workd[k] = 0;

  a_int lworkl = ncv *  (3 * ncv + 5);
  double _Complex workl[lworkl];
  for (k = 0; k < lworkl; ++k) workl[k] = 0;

  double rwork[ncv];
  double _Complex workev[2 * ncv];
  a_int info = 0;

  iparam[0] = 1;
  iparam[2] = 10 * N;
  iparam[3] = 1;
  iparam[4] = 0;  // number of ev found by arpack.
  iparam[5] = 0;
  iparam[6] = 1;
  iparam[7] = 0;
  iparam[8] = 0;
  iparam[9] = 0;
  iparam[10] = 0;

  if (verbose) printf("starting znaupd iteration\n");

  int cnt = 0;
  while (ido != 99) {
    /* call arpack like you would have, but, use znaupd_c instead of znaupd_ */
    znaupd_c(&ido, bmat, N, which, nev, tol, resid, ncv, V, ldv, iparam, ipntr,
             workd, workl, lworkl, rwork, &info);

    zMatVecLa(&(workd[ipntr[0] - 1]), &(workd[ipntr[1] - 1]), mat, dim, ldm);
    cnt++;
  }

  if (verbose) printf("finished znaupd iteration, info: %d, numer of iterations: %d\n", info, cnt);


  if (iparam[4] != nev) {
    printf("Error: iparam[4] %d, nev %d\n", iparam[4], nev); // check number of ev found by arpack.
  }


  /* call arpack like you would have, but, use zneupd_c instead of zneupd_ */
  zneupd_c(rvec, howmny, select, d, z, ldz, sigma, workev, bmat, N, which, nev,
           tol, resid, ncv, V, ldv, iparam, ipntr, workd, workl, lworkl, rwork,
           &info);

  if (verbose) printf("finished zneupd call, info: %d\n", info);

  if (verbose) printf("copying results\n");

  // sort eigenvalues and eigenvectors
  int rank[q];
  double abs_vals[q];
  for (int i = 0; i < q; i++)
    abs_vals[i] = sqrt( pow(creal(d[i]), 2) +  pow(cimag(d[i]), 2));

  get_rank(abs_vals, rank, q);


  // rotate according to eigen result:
  double z1 = 1.0, z2 = 0.0;
  int didx, ridx, cidx;

  int evecs_rdim = dim;
  if (transpose_out) evecs_rdim = q;

  // copy results -> how to avoid?
  for (int i = 0; i < q; i++) {

    didx = rank[i];
    eval[i].r = creal(d[didx]);
    eval[i].i = cimag(d[didx]);

    if (normalize_evecs) {

      z2 = sqrt(1.0 / (1.0 + pow(creal(z[didx * dim + 0]) / cimag(z[didx * dim + 0]), 2.0)));
      z1 = - creal(z[didx * dim + 0]) / cimag(z[didx * dim + 0]) * z2;
      if ((creal(z[didx * dim + 0]) * z1 - cimag(z[didx * dim + 0]) * z2) < 0.0) {
        z1 = -z1;
        z2 = -z2;
      }
    }

    for (int j = 0; j < dim; j++){
      if (transpose_out) {ridx = j; cidx = i;} else {ridx = i; cidx = j;}

      evecs[ridx * evecs_rdim + cidx].r = creal(z[didx * dim + j]) * z1 - cimag(z[didx * dim + j]) * z2;
      if (row_evecs) {
        // complex conjugate if row eigenvectors are requested (assumes hermitian matrix)
        evecs[ridx * evecs_rdim + cidx].i = - cimag(z[didx * dim + j]) * z1 - creal(z[didx * dim + j]) * z2;
      } else {
        evecs[ridx * evecs_rdim + cidx].i = cimag(z[didx * dim + j]) * z1 + creal(z[didx * dim + j]) * z2;
      }
    }
  }

}



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

  arnoldi_eigs(mat, dim, ldm, q, COMPLEX(evals), COMPLEX(evecs), tol,
               normalize_evecs, verbose, row_evecs, transpose_out);

  SET_VECTOR_ELT(res, 0, evals);
  SET_VECTOR_ELT(res, 1, evecs);

  SEXP nms = PROTECT(allocVector(STRSXP, 2));
  SET_STRING_ELT(nms, 0, mkChar("values"));
  SET_STRING_ELT(nms, 1, mkChar("vectors"));

  setAttrib(res, R_NamesSymbol, nms);

  UNPROTECT(4);
  return res;
}

