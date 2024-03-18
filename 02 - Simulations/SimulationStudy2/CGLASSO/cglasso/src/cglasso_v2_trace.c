#include <R.h>

void F77_SUB(trace_cglasso_v2_1)(int *k, double *rho, double *lambda, int *em_nit, int *nit) {
    Rprintf("\ncglasso model number %d\n\t     rho = %11.6f\n\t  lambda = %11.6f\n\tEM-steps = %4d\n\t   steps = %4d\n", *k, *rho, *lambda, *em_nit, *nit);
}

void F77_SUB(trace_cglasso_v2_2_1)(int *k, double *rho, double * lambda) {
    Rprintf("\n*************************************************************************\n");
    Rprintf("\ncglasso model number %d\n\t     rho = %11.6f\n\t  lambda = %11.6f", *k, *rho, *lambda);
}

void F77_SUB(trace_cglasso_v2_2_2)(int *i) {
    Rprintf("\t\nstarting EM algorithm step = %d\n", *i);
}

void F77_SUB(trace_cglasso_v2_2_3)(void) {
    Rprintf("\n\tE-step: completed!\n");
}

void F77_SUB(trace_cglasso_v2_2_4)(double * lambda) {
    Rprintf("\tM-step:\n\t       fitting multivariate lasso model with lambda = %11.6f\n\n", *lambda);
}

void F77_SUB(trace_cglasso_v2_2_4_1)(void) {
    Rprintf("\t\tbcd step\tlasso steps\t||B_old - B_new||_inf\n");
}

void F77_SUB(trace_cglasso_v2_2_4_2)(int *nit, int *lnit, double *dB) {
    Rprintf("\t\t%8d\t%11d\t%25.10f\n", *nit, *lnit, *dB);
}

void F77_SUB(trace_cglasso_v2_2_4_3)(double * thr) {
    Rprintf("\tConvergence is met (threshold = %f)\n\n", *thr);
}

void F77_SUB(trace_cglasso_v2_2_5)(double * rho) {
    Rprintf("\tM-step:\n\t       fitting glasso model with rho = %11.6f\n\n", *rho);
}

void F77_SUB(trace_cglasso_v2_2_6)(void) {
    Rprintf("\tM-step: completed!\n");
}

void F77_SUB(trace_cglasso_v2_2_7)(double *thr_em, double *dB, double *dSgm, double *dBdSgm) {
    Rprintf("\n\tChecking convergence criterion (threshold = %f)\n", *thr_em);
    Rprintf("\t||B_old - B_new||_2  = %f\n", *dB);
    Rprintf("\t||Sgm_old - Sgm_new||_2 = %f\n", *dSgm);
    if(*dBdSgm <= *thr_em) Rprintf("\tConvergence criterion is met!\n");
}
