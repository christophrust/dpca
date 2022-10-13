// LAPACK stegr

#include "dpca.h"


void zMatVec(double _Complex* x, double _Complex* y, double _Complex* mat, int dim) {

  int i, j;

  for (i = 0; i < dim; i++) {
    double _Complex accum = 0.0 + 0.0 * _Complex_I;
    for (j = 0; j < dim; j++) {
      accum += x[j] * mat[j * dim + i];
      /* printf("Curr mat[%d,%d]: %f%+f\n", i,j, creal(mat[i * dim + j]), cimag(mat[i * dim + j])); */
      /* printf("Curr v[%d]: %f%+f\n", j, creal(x[j]), cimag(x[j])); */
      /* printf("Curr sum: %f%+f\n", creal(x[j] * mat[i * dim + j]), cimag(x[j] * mat[i * dim + j])); */
      /* printf("Curr accum: %f%+f\n", creal(accum), cimag(accum)); */
    }
    y[i] = accum;
  }

}

void zMatVecLa(double _Complex *x, double _Complex* y, Rcomplex* mat, int dim) {

  Rcomplex alpha;
  alpha.r = 1.0; alpha.i = 0.0;
  Rcomplex beta;
  beta.r = 0.0; beta.i = 0.0;
  int inc =1;

  F77_CALL(zgemv)("N", &dim, &dim, &alpha, mat,
                  &dim,
                  (Rcomplex*) x, &inc, &beta,
                  (Rcomplex*) y, &inc);
}


void arnoldi_eigs(Rcomplex *mat, int dim, int q,
  Rcomplex *eval, Rcomplex *evecs, double tol) {


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
  a_int iparam[11];
  a_int ipntr[14];
  for (int i=0; i< 14; i++) ipntr[i] = 0;

  double _Complex workd[3 * N];
  a_int rvec = 1;
  char howmny[] = "A";
  double _Complex* d =
      (double _Complex*)malloc((nev + 1) * sizeof(double _Complex));
  a_int select[ncv];
  for (int i = 0; i < ncv; i++) select[i] = 1;
  //double _Complex z[(N + 1) * (nev + 1)];
  double _Complex z[N  * nev];
  a_int ldz = N + 1;
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

  // we still copy the array
  double _Complex cmplx_mat[dim * dim];
  for (int i= 0; i < dim * dim; i++) cmplx_mat[i] = mat[i].r + _Complex_I * mat[i].i;

  int cnt = 0;
  while (ido != 99) {
    /* call arpack like you would have, but, use znaupd_c instead of znaupd_ */
    znaupd_c(&ido, bmat, N, which, nev, tol, resid, ncv, V, ldv, iparam, ipntr,
             workd, workl, lworkl, rwork, &info);

    zMatVec(&(workd[ipntr[0] - 1]), &(workd[ipntr[1] - 1]), cmplx_mat, dim);
    //for (int i=0; i<N; i++) printf("xVec[%d]: %f%+fi\n", i, creal(workd[ipntr[0] - 1 + i]), cimag(workd[ipntr[0] - 1 + i]));
    //for (int i=0; i<N; i++) printf("yVec[%d]: %f%+fi\n", i, creal(workd[ipntr[1] - 1 + i]), cimag(workd[ipntr[1] - 1 + i]));
    cnt++;
  }

  //printf("Info: %d\n", info);
  // printf("Number of iterations: %d or %i\n", iparam[2], cnt);

  if (iparam[4] != nev) {
    printf("Error: iparam[4] %d, nev %d\n", iparam[4], nev); // check number of ev found by arpack.
  }

  /* printf("Tol: %e\n", tol); */
  /* for (int i=0; i < N; i++) printf("presid[%i,1]: %f%+fi\n", i, creal(resid[i]), cimag(resid[i])); */
  /* for (int i=0; i < N; i++) printf("V[%i,1]: %f%+fi\n", i, creal(V[i]), cimag(V[i])); */
  /* for (int i=0; i <11; i++) printf("iparam[%i]: %i\n", i, iparam[i]); */
  /* for (int i=0; i <14; i++) printf("ipntr[%i]: %i\n", i, ipntr[i]); */
  /* for (int i=0; i < 3*N; i++) printf("workd[%i,1]: %f%+fi\n", i, creal(workd[i]), cimag(workd[i])); */
  /* for (int i=0; i < lworkl; i++) printf("workl[%i,1]: %f%+fi\n", i, creal(workl[i]), cimag(workl[i])); */
  /* printf("lworkl: %i\n", lworkl); */
  /* for (int i=0; i < ncv; i++) printf("rwork[%i,1]: %f\n", i, rwork[i]); */

  /* call arpack like you would have, but, use zneupd_c instead of zneupd_ */
  zneupd_c(rvec, howmny, select, d, z, ldz, sigma, workev, bmat, N, which, nev,
           tol, resid, ncv, V, ldv, iparam, ipntr, workd, workl, lworkl, rwork,
           &info);

  // copy results
  for (int i = 0; i < q; i++) {
    eval[i].r = creal(d[i]);
    eval[i].i = cimag(d[i]);
    for (int j = 0; j < dim; j++){
      evecs[i * dim + j].r = creal(z[i * dim + j]);
      evecs[i * dim + j].i = cimag(z[i * dim + j]);
    }
  }

}


SEXP R_arnoldi_eigs(SEXP r_mat, SEXP r_dim, SEXP r_q, SEXP r_tol) {

  Rcomplex *mat = COMPLEX(r_mat);
  int dim = *INTEGER(r_dim);
  int q = *INTEGER(r_q);
  double tol = *REAL(r_tol);

  // result objects
  SEXP res = PROTECT(allocVector(VECSXP, 2));;
  SEXP evecs = PROTECT(allocMatrix(CPLXSXP, dim, q));
  SEXP evals = PROTECT(allocVector(CPLXSXP, q));

  arnoldi_eigs(mat, dim, q, COMPLEX(evals), COMPLEX(evecs), tol);

  SET_VECTOR_ELT(res, 0, evals);
  SET_VECTOR_ELT(res, 1, evecs);

  SEXP nms = PROTECT(allocVector(STRSXP, 2));
  SET_STRING_ELT(nms, 0, mkChar("values"));
  SET_STRING_ELT(nms, 1, mkChar("vectors"));

  setAttrib(res, R_NamesSymbol, nms);

  UNPROTECT(4);
  return res;
}


SEXP R_zMatVec(SEXP r_mat, SEXP r_vec, SEXP r_dim) {
  int dim = *INTEGER(r_dim);
  SEXP res = PROTECT(allocVector(CPLXSXP, dim));
  double _Complex m[dim * dim];
  double _Complex x[dim];
  double _Complex y[dim];

  for (int i = 0; i < dim*dim; i ++)
    m[i] = COMPLEX(r_mat)[i].r + COMPLEX(r_mat)[i].i * _Complex_I;
  for (int i = 0; i < dim; i ++)
    x[i] = COMPLEX(r_vec)[i].r + COMPLEX(r_vec)[i].i * _Complex_I;

  zMatVec(x, y, m, dim);

  for (int i = 0; i < dim; i ++){
    COMPLEX(res)[i].r = creal(y[i]);
    COMPLEX(res)[i].i = cimag(y[i]);
  }

  UNPROTECT(1);
  return res;
}


