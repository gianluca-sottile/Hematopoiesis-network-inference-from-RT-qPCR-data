library(lattice)

rm(list = ls())

percorso <- "results"

mse_tbl_df <- NULL

# for(i in 1:length(files)){
  load(paste(percorso, "mse_tbl.RData", sep = "/"))
  
  mse_tbl_df <- rbind(mse_tbl_df, 
                      data.frame("mse" = c(sapply(1:2, function(j) sapply(1:5, function(r) mean(mse_tbl_Tht[1,j,r,])))),
                                 "ratio" = rep(c(1, .75, .50, .25, .1), 2L),
                                 "method" = rep(c("cjglasso", "jgl"), each = 5L)))
# }

mse_mat <- matrix(format(round(mse_tbl_df[, 1], 2), justify = "right"), 
                  nrow = 2L, ncol = 5L, byrow = TRUE, 
                  dimnames = list("method" = c("jcglasso", "jgl"),  
                                  "rho/rho[max]" = format(round(unique(mse_tbl_df[, 2]), 2), justify = "right")))[, 5:1]
mse_mat
# xtable:::xtable(mse_mat)


auc_tbl_df <- NULL
# for(i in 1:length(files)){
  load(paste(percorso, "auc_tbl.RData", sep = "/"))
  
  auc_tbl_df <- rbind(auc_tbl_df, 
                      data.frame("auc" = sapply(1:2, function(j) mean(auc_Tht[1,j,1,])),
                                 "sd" = sapply(1:2, function(j) mean(auc_Tht[2,j,1,])),
                                 "method" = c("cjglasso", "jgl")))
# }

auc_mat <- matrix(apply(format(round(auc_tbl_df[, 1:2], 2), justify = "right"), 1, function(x) paste0(x[1], " (", x[2], ")")), 
                  nrow = 2L, ncol = 1L, dimnames = list("method" = c("jcglasso", "jgl"),  "AUC (sd)"))
auc_mat
# xtable:::xtable(auc_mat)

cbind(mse_mat, auc_mat)
