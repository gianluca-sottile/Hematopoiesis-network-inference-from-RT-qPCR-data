library(lattice)

rm(list = ls())

percorsoNu0 <- "results/01_nu_0"
percorsoNuBic <- "results/02_nu_bic"
percorsoNuMax <- "results/03_nu_max"

ratio <- c(1, .75, .5, .25, .1)
N <- length(ratio)

mse_tbl_omega_df <- mse_tbl_tht_df <- mse_tbl_b_df <- NULL
auc_tbl_omega_df <- auc_tbl_tht_df <- auc_tbl_b_df <- NULL

# for(i in 1:length(files1)){
  ##### MSE #####
  ##### Theta
  load(paste(percorsoNu0, "mse_tbl_Tht.RData", sep = "/"))
  
  mse_tbl_tht_df <- rbind(mse_tbl_tht_df, 
                          data.frame("mse" = sapply(1:N, function(r) mean(mse_tbl_Tht[1,1,r,])),
                                     "ratio" = ratio,
                                     "method" = rep("nu0", N)))
  
  load(paste(percorsoNuBic, "mse_tbl_Tht.RData", sep = "/"))
  
  mse_tbl_tht_df <- rbind(mse_tbl_tht_df, 
                          data.frame("mse" = sapply(1:N, function(r) mean(mse_tbl_Tht[1,1,r,])),
                                     "ratio" = ratio,
                                     "method" = rep("nuBic", N)))
  
  load(paste(percorsoNuMax, "mse_tbl_Tht.RData", sep = "/"))
  
  mse_tbl_tht_df <- rbind(mse_tbl_tht_df, 
                          data.frame("mse" = sapply(1:N, function(r) mean(mse_tbl_Tht[1,1,r,])),
                                     "ratio" = ratio,
                                     "method" = rep("nuMax", N)))
  
  ##### Omega
  load(paste(percorsoNu0, "mse_tbl_Thtx.RData", sep = "/"))
  
  mse_tbl_omega_df <- rbind(mse_tbl_omega_df, 
                            data.frame("mse" = sapply(1:N, function(r) mean(mse_tbl_Thtx[1,1,r,])),
                                       "ratio" = ratio,
                                       "method" = rep("nu0", N)))
  
  load(paste(percorsoNuBic, "mse_tbl_Thtx.RData", sep = "/"))
  
  mse_tbl_omega_df <- rbind(mse_tbl_omega_df, 
                            data.frame("mse" = sapply(1:N, function(r) mean(mse_tbl_Thtx[1,1,r,])),
                                       "ratio" = ratio,
                                       "method" = rep("nuBic", N)))
  
  load(paste(percorsoNuMax, "mse_tbl_Thtx.RData", sep = "/"))
  
  mse_tbl_omega_df <- rbind(mse_tbl_omega_df, 
                            data.frame("mse" = sapply(1:N, function(r) mean(mse_tbl_Thtx[1,1,r,])),
                                       "ratio" = ratio,
                                       "method" = rep("nuMax", N)))
  
  ##### B
  load(paste(percorsoNu0, "mse_tbl_B.RData", sep = "/"))
  
  mse_tbl_b_df <- rbind(mse_tbl_b_df, 
                        data.frame("mse" = sapply(1:N, function(r) mean(mse_tbl_B[1,1,r,])),
                                   "ratio" = ratio,
                                   "method" = rep("nu0", N)))
  
  load(paste(percorsoNuBic, "mse_tbl_B.RData", sep = "/"))
  
  mse_tbl_b_df <- rbind(mse_tbl_b_df, 
                        data.frame("mse" = sapply(1:N, function(r) mean(mse_tbl_B[1,1,r,])),
                                   "ratio" = ratio,
                                   "method" = rep("nuBic", N)))
  
  load(paste(percorsoNuMax, "mse_tbl_B.RData", sep = "/"))
  
  mse_tbl_b_df <- rbind(mse_tbl_b_df, 
                        data.frame("mse" = sapply(1:N, function(r) mean(mse_tbl_B[1,1,r,])),
                                   "ratio" = ratio,
                                   "method" = rep("nuMax", N)))
  
  ##### AUC #####
  ##### Theta
  load(paste(percorsoNu0, "auc_tbl_Tht.RData", sep = "/"))
  
  auc_tbl_tht_df <- rbind(auc_tbl_tht_df, 
                          data.frame("auc" = sapply(1:N, function(r) mean(rowMeans(auc_tbl_Tht[,1,,r]))),
                                     "ratio" = ratio,
                                     "method" = rep("nu0", N)))
  
  load(paste(percorsoNuBic, "auc_tbl_Tht.RData", sep = "/"))
  
  auc_tbl_tht_df <- rbind(auc_tbl_tht_df, 
                          data.frame("auc" = sapply(1:N, function(r) mean(rowMeans(auc_tbl_Tht[,1,,r]))),
                                     "ratio" = ratio,
                                     "method" = rep("nuBic", N)))
  
  load(paste(percorsoNuMax, "auc_tbl_Tht.RData", sep = "/"))
  
  auc_tbl_tht_df <- rbind(auc_tbl_tht_df, 
                          data.frame("auc" = sapply(1:N, function(r) mean(rowMeans(auc_tbl_Tht[,1,,r]))),
                                     "ratio" = ratio,
                                     "method" = rep("nuMax", N)))
  
  ##### Omega
  load(paste(percorsoNu0, "auc_tbl_Thtx.RData", sep = "/"))
  
  auc_tbl_omega_df <- rbind(auc_tbl_omega_df, 
                            data.frame("auc" = sapply(1:N, function(r) mean(rowMeans(auc_tbl_Thtx[,1,,r]))),
                                       "ratio" = ratio,
                                       "method" = rep("nu0", N)))
  
  load(paste(percorsoNuBic, "auc_tbl_Thtx.RData", sep = "/"))
  
  auc_tbl_omega_df <- rbind(auc_tbl_omega_df, 
                            data.frame("auc" = sapply(1:N, function(r) mean(rowMeans(auc_tbl_Thtx[,1,,r]))),
                                       "ratio" = ratio,
                                       "method" = rep("nuBic", N)))
  
  load(paste(percorsoNuMax, "auc_tbl_Thtx.RData", sep = "/"))
  
  auc_tbl_omega_df <- rbind(auc_tbl_omega_df, 
                            data.frame("auc" = sapply(1:N, function(r) mean(rowMeans(auc_tbl_Thtx[,1,,r]))),
                                       "ratio" = ratio,
                                       "method" = rep("nuMax", N)))
  
  ##### B
  load(paste(percorsoNu0, "auc_tbl_B.RData", sep = "/"))
  
  auc_tbl_b_df <- rbind(auc_tbl_b_df, 
                        data.frame("auc" = sapply(1:N, function(r) mean(rowMeans(auc_tbl_B[,1,,r]))),
                                   "ratio" = ratio,
                                   "method" = rep("nu0", N)))
  
  load(paste(percorsoNuBic, "auc_tbl_B.RData", sep = "/"))
  
  auc_tbl_b_df <- rbind(auc_tbl_b_df, 
                        data.frame("auc" = sapply(1:N, function(r) mean(rowMeans(auc_tbl_B[,1,,r]))),
                                   "ratio" = ratio,
                                   "method" = rep("nuBic", N)))
  
  load(paste(percorsoNuMax, "auc_tbl_B.RData", sep = "/"))
  
  auc_tbl_b_df <- rbind(auc_tbl_b_df, 
                        data.frame("auc" = sapply(1:N, function(r) mean(rowMeans(auc_tbl_B[,1,,r]))),
                                   "ratio" = ratio,
                                   "method" = rep("nuMax", N)))
  
# }

mse_tbl_tht_df$method <- factor(mse_tbl_tht_df$method, levels = c("nu0", "nuBic", "nuMax"))
mse_tbl_omega_df$method <- factor(mse_tbl_omega_df$method, levels = c("nu0", "nuBic", "nuMax"))
mse_tbl_b_df$method <- factor(mse_tbl_b_df$method, levels = c("nu0", "nuBic", "nuMax"))

auc_tbl_tht_df$method <- factor(auc_tbl_tht_df$method, levels = c("nu0", "nuBic", "nuMax"))
auc_tbl_omega_df$method <- factor(auc_tbl_omega_df$method, levels = c("nu0", "nuBic", "nuMax"))
auc_tbl_b_df$method <- factor(auc_tbl_b_df$method, levels = c("nu0", "nuBic", "nuMax"))


MSE <- rbind(
  cbind.data.frame(mse_tbl_tht_df, "measure" = "mse_tht"),
  cbind.data.frame(mse_tbl_omega_df, "measure" = "mse_omega"),
  cbind.data.frame(mse_tbl_b_df, "measure" = "mse_b")
)
AUC <- rbind(
  cbind.data.frame(auc_tbl_tht_df, "measure" = "auc_tht"),
  cbind.data.frame(auc_tbl_omega_df, "measure" = "auc_omega"),
  cbind.data.frame(auc_tbl_b_df, "measure" = "auc_b")
)
colnames(MSE)[1] <- "value"
colnames(AUC)[1] <- "value"
ALL <- rbind(MSE, AUC)
ALL$measure <- factor(ALL$measure, levels = c("auc_tht", "auc_omega", "auc_b",
                                              "mse_tht", "mse_omega", "mse_b"))

flevels <- c(
  expression(paste(AUC[rho], "(", hat(Theta), ")")),
  expression(paste(AUC[nu], "(", hat(Omega), ")")),
  expression(paste(AUC[lambda], "(", hat(B), ")")),
  expression(paste(min[rho], "MSE(", hat(Theta), ")")),
  expression(paste(min[nu], "MSE(", hat(Omega), ")")),
  expression(paste(min[lambda], "MSE(", hat(B), ")"))
)

flevels2 <- c(
  expression(lambda / lambda[max]),
  expression(lambda / lambda[max]),
  expression(rho / rho[max])
)

key.method <- list(text = list(c(expression(paste(nu, " = ", nu[BIC])), 
                                 expression(paste(nu, " = ", 0)), 
                                 expression(paste(nu, " = ", nu[max])))),
                   points = list(pch = 1:3, col = "black"), column = 3, row = 1)

ALL2 <- ALL
ALL2[grep("mse_", ALL2$measure), 1] <- c(unlist(tapply(ALL[grep("mse_tht", ALL$measure), "value"], 
                                                          ALL[grep("mse_tht", ALL$measure), "method"], 
                                                          function(x) x / max(x))),
                                            unlist(tapply(ALL[grep("mse_omega", ALL$measure), "value"], 
                                                          ALL[grep("mse_omega", ALL$measure), "method"], 
                                                          function(x) x / max(x))),
                                            unlist(tapply(ALL[grep("mse_b", ALL$measure), "value"], 
                                                          ALL[grep("mse_b", ALL$measure), "method"], 
                                                          function(x) x / max(x))))


# pdf(file = "figs/fig1supp_mse_auc_nu.pdf", width = .8*9, height = .8*7)
xyplot(value ~ ratio | measure, data = ALL2, type = "o", groups = method, layout = c(3,2),
       lty = 1, pch = 1:3, col.line = "darkgrey", col.symbol = "black", 
       ylab = "", xlab = flevels2,#expression(lambda / lambda[max]), 
       strip = strip.custom(bg = "grey92", factor.levels = flevels, 
                            par.strip.text = list(cex = 0.9)), between = list(x = 1, y = 1), 
       scales = list(alternating = 2, relation = "free", tick.number = 4), key = key.method, 
       # ylim = list(c(.35,.95), c(.35,.95), c(.35,.95), c(0,35), c(0,35), c(0,35)))#mod6
       ylim = list(c(.20,.90), c(.20,.90), c(.20,.90), c(.1,1), c(.1,1), c(.1,1)))#mod8
# dev.off()

