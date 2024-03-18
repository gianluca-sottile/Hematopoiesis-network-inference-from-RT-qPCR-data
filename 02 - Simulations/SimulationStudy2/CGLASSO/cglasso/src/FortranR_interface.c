/*  wrappers for calling R's functions from the  FORTRAN code */

#include "cglasso.h"

double F77_SUB(rdnorm)(double const *x, double const *mu, double const *sigma, int const *give_log){
    return dnorm4( *x, *mu, *sigma, *give_log);
}

double F77_SUB(rpnorm)(double const *x, double const *mu, double const *sigma, int const *lower_tail, int const *log_p){
    return pnorm5( *x, *mu, *sigma, *lower_tail, *log_p);
}
