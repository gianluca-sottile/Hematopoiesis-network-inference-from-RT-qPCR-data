#include <R_ext/RS.h>

#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>

void
F77_SUB(setup)(int *n, int *p, double *Y, double *lo, double *up,
               double *Yna, int *R, int *startmis, int *order);

void
F77_SUB(fitmcgm)(int *n, int *p, double *Y, double *lo, double *up,
                 int *R, int *nstp, double *eps, double *ym, double *yv,
                 int *conv);

void
F77_SUB(cglasso_v1)(int *n, int *p, double *Y, int *Id, int *nP,
                    int *InfoP, double *lo, double *up, int *pendiag, double *wTht,
                    double *ym, double *yv, int *nrho, double *rhoratio, double *rho,
                    int *maxit_em, double *thr_em, int *maxit_bcd, double *thr_bcd, double *Yipt,
                    double *B, double *mu, double *R, double *S, double *Sgm,
                    double *Tht, int *Adj_yy, int *dfB, int *dfTht, int *ncomp,
                    int *Ck, int *pk, int *nit, int *conv, int *subrout,
                    int *trace);

void
F77_SUB(cglasso_v2)(int *n, int *q, double *X, int *p, double *Y,
                    int *Id, int *nP, int *InfoP, double *lo, double *up,
                    double *wB, int *pendiag, double *wTht, double *ym, double *yv,
                    int *nlambda, double *lambdaratio, double *lambda, int *nrho, double *rhoratio,
                    double *rho, int *maxit_em, double *thr_em, int *maxit_bcd, double *thr_bcd,
                    double *Yipt, double *B, double *mu, double *R, double *S,
                    double *Sgm, double *Tht, int *Adj_yy, int *dfB, int *dfTht,
                    int *ncomp, int *Ck, int *pk, int *Adj_xy, int *nit,
                    int *conv, int *subrout, int *trace);

void
F77_SUB(predict)(double *newrho, double *newlambda, int *nrho, double *rho, int *nlambda,
                 double *lambda, int *d1, int *d2, double *Min, double *Mout);

void
F77_SUB(impute)(double *newrho, double *newlambda, int *nrho, double *rho, int *nlambda,
                double *lambda, int *n, int *p, double *Yin, int *Id,
                double *Yout);

void
F77_SUB(cggm_v1)(int *n, int *p, double *Y, int *Id, int *nP,
                 int *InfoP, double *lo, double *up, double *ym, double *yv,
                 int *pendiag, double *wTht, int *ntp, double *rho, int *maxit_em,
                 double *thr_em, int *maxit_bcd, double *thr_bcd, double *Yipt, double *B,
                 double *mu, double *R, double *S, double *Sgm, double *Tht,
                 int *nit, int *conv, int *subrout, int *trace);

void
F77_SUB(cggm_v2)(int *n, int *q, double *X, int *p, double *Y,
                 int *Id, int *nP, int *InfoP, double *lo, double *up,
                 double *ym, double *yv, double *wB, int *pendiag, double *wTht,
                 int *ntp, double *lambda, double *rho, int *maxit_em, double *thr_em,
                 int *maxit_bcd, double *thr_bcd, double *Yipt, double *B, double *mu,
                 double *R, double *S, double *Sgm, double *Tht, int *nit,
                 int *conv, int *subrout, int *trace);

void
  F77_SUB(cglasso_v3)(int *n, int *p, double *Y, int *Id, int *nP,
          int *InfoP, double *lo, double *up, int *pendiag, double *wTht,
          double *ym, double *yv, int *nrho, double *rhoratio, double *rho,
          int *maxit_em, double *thr_em, int *maxit_bcd, double *thr_bcd, double *Yipt,
          double *B, double *mu, double *R, double *S, double *Sgm,
          double *Tht, int *Adj_yy, int *dfB, int *dfTht, int *ncomp,
          int *Ck, int *pk, int *nit, int *conv, int *subrout,
          int *trace);


void
  F77_SUB(cglasso_v4)(int *n, int *q, double *X, int *p, double *Y,
          int *Id, int *nP, int *InfoP, double *lo, double *up,
          double *wB, int *pendiag, double *wTht, double *ym, double *yv,
          int *nlambda, double *lambdaratio, double *lambda, int *nrho, double *rhoratio,
          double *rho, int *maxit_em, double *thr_em, int *maxit_bcd, double *thr_bcd,
          double *Yipt, double *B, double *mu, double *R, double *S,
          double *Sgm, double *Tht, int *Adj_yy, int *dfB, int *dfTht,
          int *ncomp, int *Ck, int *pk, int *Adj_xy, int *nit,
          int *conv, int *subrout, int *trace);

void
F77_SUB(cggm_v3)(int *n, int *p, double *Y, int *Id, int *nP,
                 int *InfoP, double *lo, double *up, double *ym, double *yv,
                 int *pendiag, double *wTht, int *ntp, double *rho, int *maxit_em,
                 double *thr_em, int *maxit_bcd, double *thr_bcd, double *Yipt, double *B,
                 double *mu, double *R, double *S, double *Sgm, double *Tht,
                 int *nit, int *conv, int *subrout, int *trace);

void
F77_SUB(cggm_v4)(int *n, int *q, double *X, int *p, double *Y,
                 int *Id, int *nP, int *InfoP, double *lo, double *up,
                 double *ym, double *yv, double *wB, int *pendiag, double *wTht,
                 int *ntp, double *lambda, double *rho, int *maxit_em, double *thr_em,
                 int *maxit_bcd, double *thr_bcd, double *Yipt, double *B, double *mu,
                 double *R, double *S, double *Sgm, double *Tht, int *nit,
                 int *conv, int *subrout, int *trace);

void
  F77_SUB(admm)(int *p, double *S, double *initOmega, double *initZ, double *initY,
          double *lam, double *thr, int *maxit, double *Sgm, int *iter, int *trace);

void
  F77_SUB(glassosub_v2)(int *p, double *S, int *pendiag, double *rho, double *wTht, int *maxit, 
          double *thr, double *Tht, int *k, int *Ck, int *pk, 
          int *nit, int *conv, int *trace);

void
F77_SUB(e_step_v2)(int *n, int *p, double *Y, double *lo, double *up, int *nP, int *InfoP,
                   double *T1o, double *T2o, double *mu, double *Sgm, double *Tht, double *Yipt_lo,
                   double *Yipt_up, double *Yipt, double *T1, double *T2, int *conv);

void
F77_SUB(multilasso)(int *p, int *q, double *ym, double *xm, double *xtx_n, double *xtr_n, double *B,
                    double *Tht, double *W, double *lmb, int *maxit, double *thr, int *conv,
                    int *subrout, int *nit, int *df, int *trace);

void
  F77_SUB(admm_multilasso_2)(int *p, int *q, int *K, int *qp, double *ym, double *xm, double *xtx_n, 
                            double *xtr_n, double *B, double *Tht, double *W, double *lmb1, double *lmb2, 
                            int *maxit, double *thr, double *rho, int *conv, int *subrout, int *nit, 
                            int *df, int *trace);

void
  F77_SUB(admm_b_single)(int *p, int *q, int *qp, double *xtx_n,
                            double *xtr_n, double *B, double *Tht, double *W, double *lmb,
                            int *maxit, double *thr, double *rho, int *trace, int *conv, int *subrout, int *nit,
                            int *df);

void
  F77_SUB(admm_tht_single)(int *p, double *S,
                            double *Theta, double *wTht, double *lam,
                            int *maxit, double *thr, double *rho, int *trace, int *conv, int *nit,
                            int *df);

void
  F77_SUB(admm_b_joint)(int *p, int *q, int *qp, int *K, double *fk, double *xtx_n,
                            double *xtr_n, double *B, double *Tht, double *W, double *lmb, double *alpha,
                            int *maxit, double *thr, double *rho, int *trace, int *conv, int *subrout, int *nit,
                            int *df);

void
  F77_SUB(admm_tht_joint)(int *p, int *K, double *fk, double *S,
                            double *Theta, double *wTht, double *lam, double *alpha,
                            int *maxit, double *thr, double *rho, int *trace, int *conv, int *nit,
                            int *df);

void
  F77_SUB(graph_adjacency)(int *p, int *N, double *fk, double *S,
                            double *W, double *rho, double *alpha, int *k,
                            int *Ck, int *pk);

void
  F77_SUB(admm_tht_sub)(int *p, int *N, double *fk, double *S, double *wTht, int *pendiag, double *rho, double *alpha, int *maxit, double *thr, double *Tht, int *k, int *Ck, int *pk, int *nit, int *df, int *conv, int *trace);

void F77_SUB(grad_b)(int *q, int *p, int *n, double *m, double *f, double *A, double *b, double *x, double *Tht, double *g);
void F77_SUB(prox_lasso)(int *p, int *q, int *k, double *x, double *lam);
void F77_SUB(prox_grouplasso)(int *p, int *q, int *k, double *x, double *lam);
void F77_SUB(prox_sparsegrouplasso)(int *p, int *q, int *k, double *x, double *lam, double *alpha);
void F77_SUB(apg)(int *p, int *q, int *n, int *k, double *nk, double *fk, double *A, double *b, double *Tht, double *weights, double *lambda, double *alpha, double *xm, double *ym, double *x, double *beta,int *maxit, double *thr, int *trace, int *df, int *nit, int *conv);

void F77_SUB(softmatrix_b)(double *S, double *Tau, int *p, int *q);
void F77_SUB(softmatrix_b_joint)(double *S, double *Tau, double *norm, int *p, int *q);
