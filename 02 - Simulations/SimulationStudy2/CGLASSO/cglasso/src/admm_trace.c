#include <R.h>


void F77_SUB(admm_trace_1_2)(void) {
    Rprintf("\t\tadmm step\trho value\t||Tht_old - Tht_new||_2\n");
}
void F77_SUB(admm_trace_1)(int *nit, double *Db, double *rho) {
    Rprintf("\t\t%8d\t%11.5f\t%25.10f\n", *nit, *rho, *Db);
}

void F77_SUB(admm_trace_2)(int *nit, double *r, double *eps1, double *s, double *eps2, double *rho) {
    Rprintf("\titer num %6d\trho: %6.4f\tr: %7.6f\teps1: %7.6f\ts: %7.6f\teps2: %7.6f\n", *nit, *rho, *r, *eps1, *s, *eps2);
}

void F77_SUB(apg_trace)(int *i, double *err1, double *t) {
    Rprintf("\titer num %6d\tnorm(tGk): %11.8f\tstep-size: %7.4f\n", *i, *err1, *t);
}
