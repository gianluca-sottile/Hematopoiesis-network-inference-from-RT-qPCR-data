#include <stdlib.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Visibility.h>

#include "cglasso.h"

static const R_FortranMethodDef FortEntries[] = {
    {"setup", (DL_FUNC) &F77_SUB(setup), 9},
    {"fitmcgm", (DL_FUNC) &F77_SUB(fitmcgm), 11},
    {"cglasso_v1", (DL_FUNC) &F77_SUB(cglasso_v1), 36},
    {"cglasso_v2", (DL_FUNC) &F77_SUB(cglasso_v2), 43},
    {"predict", (DL_FUNC) &F77_SUB(predict), 10},
    {"impute", (DL_FUNC) &F77_SUB(impute), 11},
    {"cggm_v1", (DL_FUNC) &F77_SUB(cggm_v1), 29},
    {"cggm_v2", (DL_FUNC) &F77_SUB(cggm_v2), 33},
    {"multilasso", (DL_FUNC) &F77_SUB(multilasso), 17},
    
    {"e_step_v2", (DL_FUNC) &F77_SUB(e_step_v2), 18},
    {"cglasso_v3", (DL_FUNC) &F77_SUB(cglasso_v3), 36},
    {"cglasso_v4", (DL_FUNC) &F77_SUB(cglasso_v4), 43},
    {"cggm_v3", (DL_FUNC) &F77_SUB(cggm_v3), 29},
    {"cggm_v4", (DL_FUNC) &F77_SUB(cggm_v4), 33},
    {"admm", (DL_FUNC) &F77_SUB(admm), 11},
    {"glassosub_v2", (DL_FUNC) &F77_SUB(glassosub_v2), 14},
    {"admm_multilasso_2", (DL_FUNC) &F77_SUB(admm_multilasso_2), 21},
    {"admm_b_single", (DL_FUNC) &F77_SUB(admm_b_single), 17},
    {"admm_tht_single", (DL_FUNC) &F77_SUB(admm_tht_single), 12},
    {"admm_b_joint", (DL_FUNC) &F77_SUB(admm_b_joint), 20},
    {"admm_tht_joint", (DL_FUNC) &F77_SUB(admm_tht_joint), 15},
    {"graph_adjacency", (DL_FUNC) &F77_SUB(graph_adjacency), 10},
    {"admm_tht_sub", (DL_FUNC) &F77_SUB(admm_tht_sub), 18},
    {"grad_b", (DL_FUNC) &F77_SUB(grad_b), 10},
    {"prox_lasso", (DL_FUNC) &F77_SUB(prox_lasso), 5},
    {"prox_grouplasso", (DL_FUNC) &F77_SUB(prox_grouplasso), 5},
    {"prox_sparsegrouplasso", (DL_FUNC) &F77_SUB(prox_sparsegrouplasso), 6},

    {"softmatrix_b", (DL_FUNC) &F77_SUB(softmatrix_b), 4},
    {"softmatrix_b_joint", (DL_FUNC) &F77_SUB(softmatrix_b_joint), 5},

    {"apg", (DL_FUNC) &F77_SUB(apg), 22},
    {NULL, NULL, 0}
};

void attribute_visible R_init_cglasso(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, NULL, FortEntries, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
