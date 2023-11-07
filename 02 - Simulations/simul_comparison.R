rm(list = ls())

# libraries
library(MASS)
library(cglasso)
library(JGL)
library(huge)
library(ggplot2)
library(caret)
setwd('~/Dropbox/group-cglasso/03 - R/01 - RCode/')

source("datajcggm.R")
source("jcglasso.R")
source("jcggm.R")
source("gof.R")
source("to_graph.R")

simulTheta <- function(p, ncomp, unconnected = NULL, mintht = 0.3, maxtht = 0.5){
  if(p %% ncomp != 0) warning("problem in dimensions p and ncomp")
  frac <- p / ncomp
  block <- matrix(0.0, ncomp, ncomp)
  block[1L, ] <- block[, 1L] <- sort(runif(ncomp, min = mintht, maxtht), TRUE)
  diag(block) <- 1.0
  subTht <- diag(frac)
  Tht <- kronecker(subTht, block)
  
  if(!is.null(unconnected)) {
    tmp <- which(Tht[upper.tri(Tht)] != 0)
    id <- sample(tmp, size = floor(unconnected * length(tmp)))
    Tht[upper.tri(Tht)][id] <- 0.0
    Tht <- t(Tht)
    Tht[upper.tri(Tht)][id] <- 0.0
  }
  
  diag(Tht) <- 1.0
  Tht
}
simulOmega <- function(p, ncomp, unconnected = NULL, mintht = 0.3, maxtht = 0.5){
  if(p %% ncomp != 0) warning("problem in dimensions p and ncomp")
  frac <- p / ncomp
  block <- matrix(0.0, ncomp, ncomp)
  block[1L, ] <- block[, 1L] <- sort(runif(ncomp, min = mintht, maxtht), TRUE)
  diag(block) <- 1.0
  subTht <- diag(frac)
  Tht <- kronecker(subTht, block)
  
  if(!is.null(unconnected)) {
    tmp <- which(Tht[upper.tri(Tht)] != 0)
    id <- sample(tmp, size = floor(unconnected * length(tmp)))
    Tht[upper.tri(Tht)][id] <- 0.0
    Tht <- t(Tht)
    Tht[upper.tri(Tht)][id] <- 0.0
  }
  
  diag(Tht) <- 1.0
  Tht
}
simulB <- function(X, perc, up, q, p, Bmin = .5, Bmax = .7, indexB, unconnectedB) {
  B <- matrix(0.0, nrow = q + 1L, ncol = p)
  for(j in seq_len(p)) B[1L + indexB[, j], j] <- runif(2, Bmin, Bmax)
  
  eta <- X %*% B[-1L, ]
  cutoff <- function(b0, i, up, perc = perc) (mean(pnorm(up - b0 - eta[, i], lower.tail =  FALSE))) - perc
  B[1L, ] <- sapply(1:p, function(.i) uniroot(cutoff, interval = c(0, 100), i = .i, up = up, perc = perc[.i])$root)
  
  if(!is.null(unconnectedB)) {
    tmp <- which(B[-1L, ] != 0)
    id <- sample(tmp, size = floor(unconnectedB * length(tmp)), replace = FALSE)
    B[-1L, ][id] <- 0.0
  }
  
  return(B)
}

lmbMax_fun <- function(obj, lambda.id = 1L, rho.id = 1L) {
  K <- length(obj$Z)
  seq_k <- seq_len(K)
  n <- obj$nobs
  q <- obj$npred
  p <- obj$nresp
  fk <- proportions(n)
  
  id_Y <- obj$InfoStructure$id_Y
  id_X <- obj$InfoStructure$id_X
  
  lambda_max <- max(sapply(seq_k, function(k) {
    X <- obj$Zipt[seq_len(n[k]), id_X, k, lambda.id, rho.id]
    Y <- obj$Zipt[seq_len(n[k]), id_Y, k, lambda.id, rho.id]
    B <- obj$B[, , k, lambda.id, rho.id]
    Ycenter <- Y - cbind(1, X) %*% B
    Sxy <- crossprod(X, Ycenter) / n[k]
    Thtyy <- obj$Tht[id_Y, id_Y, k, lambda.id, rho.id]
    # max(fk[k] * abs(Sxy %*% Thtyy / matrix(diag(Thtyy), nrow = q, ncol = p)))
    max(fk[k] * abs(Sxy %*% Thtyy))
  }))
  lambda_max
}
rhoMax_fun <- function(obj, lambda.id = 1L, rho.id = 1L) {
  K <- length(obj$Z)
  seq_k <- seq_len(K)
  n <- obj$nobs
  fk <- proportions(n)
  p <- obj$nresp
  U <- outer(1:p, 1:p, "<")
  
  id_Y <- obj$InfoStructure$id_Y
  id_X <- obj$InfoStructure$id_X
  
  rho_max <- max(sqrt(rowSums(sapply(seq_k, function(k) {
    X <- obj$Zipt[seq_len(n[k]), id_X, k, lambda.id, rho.id]
    Y <- obj$Zipt[seq_len(n[k]), id_Y, k, lambda.id, rho.id]
    B <- obj$B[, , k, lambda.id, rho.id]
    Ycenter <- Y - cbind(1, X) %*% B
    Syy <- crossprod(Ycenter) / n[k]
    pmax(fk[k] * abs(Syy)[U], 0)})^2)))
  
  rho_max
}

PrecRecAUC <- function(A, B, plot.it = FALSE) {
  dimA <- dim(A)
  n <- dimA[3]
  p <- dimA[2]
  q <- dimA[1]
  U <- if(q == p) outer(1:p, 1:p, "<=") else matrix(TRUE, q, p)
  tmp <- sapply(1:n, \(i) {
    data <- factor(1*c(A[, , i][U] != 0), levels = c(1,0))
    reference <- factor(1 * c(B[U] != 0), levels = c(1,0))
    confusionMatrix(data = data, reference = reference, mode = "prec_recall")$byClass[c("Precision", "Recall")]
  })
  tmp <- cbind(c(1, sum(ncol(B) != 0) / sum(B[U] != 0)), 
               tmp, 
               c(sum(B[U] != 0) / length(B[U]), 1))
  splFun <- suppressWarnings(splinefun(tmp["Recall", ], tmp["Precision", ], method = "monoH.FC"))
  if(plot.it) {
    plot(tmp["Recall", ], tmp["Precision", ], type = "b", xlim = c(0, 1), ylim = c(0, 1))
    x.seq <- seq(max(tmp["Recall", ]), min(tmp["Recall", ]), length = 100L)
    lines(x.seq, splFun(x.seq), lty = 2, col = "red")
  }
  auc <- integrate(splFun, lower = min(tmp["Recall", ]), upper = max(tmp["Recall", ]))$value
  list(precision = tmp["Precision", ], recall = tmp["Recall", ], auc = auc)
}

impute <- function(x) {
  require(missForest)
  tmp <- missForest(x)
  tmp$ximp
  # idNA <- is.na(x)
  # nNA <- sum(idNA)
  # if(nNA > 0) x[idNA] <- mean(x, na.rm = TRUE)
  # x
}
JSEMcust <- function(Y, Theta.Y, lambda.vec = NULL) {
  source('JMMLE/JMLE.R')
  source('JMMLE/l1LS_Main.R')
  source('JMMLE/jsem.R')
  source('JMMLE/Objval.R')
  require(glmnet)
  require(grpreg)
  require(glasso)
  require(parallel)
  
  K <- length(Y)
  n <- nrow(Y[[1]])
  p <- ncol(Y[[1]])
  
  Y.indices <- lapply(seq_len(K), \(k) rep.int(k, times = n))
  Y.groups <- lapply(seq_len(p), \(i) matrix(1.0, nrow = K, ncol = p - 1))
  
  Theta.group.array <- list2array(lapply(seq_len(K), \(k) {
    tmp <- matrix(1.0, p, p)
    diag(tmp) <- 0.0
    tmp
  }))
  
  if(is.null(lambda.vec)){
    lambda.vec <- sqrt(log(p)/n) * seq(1.8, 0.4, -0.2)
  }
  nlambda <- length(lambda.vec)
  
  Y2 <- lapply(Y, scale, center = TRUE, scale = TRUE)
  
  Theta.array <- array(0, c(p, p, K, nlambda))
  
  cl <- makeCluster(4L)
  clusterExport(cl, c("bdiag", "sort.variables", "symmetrize",
                      "JSEM", "grpreg",
                      "zeroInd", "list2array",
                      "multi.glasso", "glasso"))
  tmp <- parLapply(cl = cl, X = 1:nlambda, fun = function(l) {
    init.jsem.model <- JSEM(do.call(rbind, Y2), unlist(Y.indices),
                            Y.groups, lambda = lambda.vec[l])
    Ahat <- init.jsem.model$Ahat
    Info <- list()
    for (k in 1:K) Info[[k]] <- zeroInd(Ahat[[k]], 1)$zeroArr
    Theta_refit <- multi.glasso(do.call(rbind, Y2), unlist(Y.indices), 1E-8, Info)
    list2array(Theta_refit$Theta)
  })
  stopCluster(cl)
  
  for(l in 1:nlambda) Theta.array[, , , l] <- tmp[[l]]
  
  ## get all models
  # for(l in 1:nlambda){
  #   init.jsem.model <- JSEM(do.call(rbind, Y2), unlist(Y.indices),
  #                           Y.groups, lambda = lambda.vec[l])
  #   Ahat <- init.jsem.model$Ahat
  #   Info <- list()
  #   for (k in 1:K) Info[[k]] <- zeroInd(Ahat[[k]], 1)$zeroArr
  #   Theta_refit <- multi.glasso(do.call(rbind, Y2), unlist(Y.indices), 1E-8, Info)
  #   Theta.array[, , , l] <- list2array(Theta_refit$Theta)
  # } 
  
  return(Theta.array)
}
JMMLEcust <- function(X, Theta.X, Y, Theta.Y, B0, lambda.vec = NULL, gamma.vec = NULL) {
  source('JMMLE/JMLE.R')
  source('JMMLE/l1LS_Main.R')
  source('JMMLE/jsem.R')
  source('JMMLE/Objval.R')
  require(glmnet)
  require(grpreg)
  require(glasso)
  require(parallel)
  
  K <- length(X)
  n <- nrow(X[[1]])
  q <- ncol(X[[1]])
  p <- ncol(Y[[1]])
  
  X.indices <- lapply(seq_len(K), \(k) rep.int(k, times = n))
  X.groups <- lapply(seq_len(q), \(i) matrix(1.0, nrow = K, ncol = q - 1))
  Y.indices <- lapply(seq_len(K), \(k) rep.int(k, times = n))
  Y.groups <- lapply(seq_len(p), \(i) matrix(1.0, nrow = K, ncol = p - 1))
  
  B0.group.array <- list2array(lapply(seq_len(K), \(k) 
                                      matrix(seq_len(p * q), nrow = q, ncol = p, byrow = TRUE)))
  Theta.group.array <- list2array(lapply(seq_len(K), \(k) {
    tmp <- matrix(1.0, p, p)
    diag(tmp) <- 0.0
    tmp
  }))
  
  if(is.null(lambda.vec)){
    lambda.vec <- sqrt(log(p)/n) * seq(1.8, 0.4, -0.2)
  }
  if(is.null(gamma.vec)){
    gamma.vec <- sqrt(log(q)/n) * seq(1, 0.4, -0.1)
  }
  nlambda <- length(lambda.vec)
  ngamma <- length(gamma.vec)
  model.list<- vector("list", nlambda)
  
  X2 <- lapply(X, scale, center = TRUE, scale = TRUE)
  Y2 <- lapply(Y, scale, center = TRUE, scale = TRUE)
  
  cl <- makeCluster(4L)
  clusterExport(cl, c("bdiag", "sort.variables", "symmetrize",
                      "JSEM", "grpreg", 
                      "sel.lambda.jsem", "matTr",
                      "l1LS_Main", "glmnet", "huge",
                      "jmmle.1step", "Obj", "squaredError",
                      "zeroInd", "list2array",
                      "multi.glasso", "glasso"))
  
  for(l in 1:nlambda){
    model.list[[l]] <- parLapply(cl = cl, X = 1:ngamma, fun = function(g) {
      jmmle.1step(Y.list = Y2, Y.indices = Y.indices, X.list = X2, 
                  B.group.array = B0.group.array, Theta.groups = Y.groups,
                  lambda = lambda.vec[l], gamma = gamma.vec[g], 
                  init.option = 1, tol = 1e-3, VERBOSE = FALSE, refit.B = TRUE)
    })
  }
  stopCluster(cl)
  
  # ## get all models
  # for(l in 1:nlambda){
  #   for(g in 1:ngamma) {
  #     if(g == 1) {
  #       model.list[[l]][[g]] <- jmmle.1step(Y.list = Y2, Y.indices = Y.indices, X.list = X2, 
  #                                           B.group.array = B0.group.array, Theta.groups = Y.groups,
  #                                           lambda = lambda.vec[l], gamma = gamma.vec[g], 
  #                                           init.option = 1, tol = 1e-3, VERBOSE = FALSE, refit.B = TRUE)
  #     } 
  #     else {
  #       model.list[[l]][[g]] <- jmmle.1step(Y.list = Y2, Y.indices = Y.indices, X.list = X2, 
  #                                           B.group.array = B0.group.array, Theta.groups = Y.groups,
  #                                           lambda = lambda.vec[l], gamma = gamma.vec[g], 
  #                                           init.gamma = gamma.vec[g], B_init.array = model.list[[l]][[g - 1L]]$B.refit, 
  #                                           Theta_init.array = list2array(model.list[[l]][[g - 1L]]$Theta_refit$Theta),
  #                                           init.option = 0, tol = 1e-3, VERBOSE = FALSE, refit.B = TRUE)
  #     }
  #   }
  # } 
  # 
  B.refit <- array(0.0, dim = c(q + 1L, p, K, nlambda, ngamma))
  Theta.refit <- array(0.0, dim = c(p, p, K, nlambda, ngamma))
  
  for(l in 1:nlambda) {
    for(g in 1:ngamma) {
      B.refit[-1, , , l, g] <- model.list[[l]][[g]]$B.refit
      for(k in 1:K) {
        B.refit[1, , k, l, g] <- attr(Y2[[k]], "scaled:center") - attr(X2[[k]], "scaled:center") %*% B.refit[-1, , k, l, g]
        Theta.refit[ , , k, l, g] <- model.list[[l]][[g]]$Theta_refit$Theta[[k]]
      }
    }
  }
  
  return(list("B" = B.refit, "Theta" = Theta.refit))
}

n <- 100L                # sample size
p <- 25L               # number of response variables
q <- 10L                 # number of predictors
ncomp <- 5L
n.sim <- 50L

c2 <- c(1, .75, .5, .25)
c3 <- seq(1, .1, l = 10)
Bmin <- 0.5
Bmax <- 0.7
Thtmin <- 0.3
Thtmax <- 0.5
Omgmin <- 0.3
Omgmax <- 0.5

alpha.seq <- c(.25, .50, .75)
alpha.true <- 0.50

K <- 2L
perc.NA <- 0.25
perc.X.na <- 0.25
perc.Y.na <- 0.25

up <- 40
perc <- rep(1E-6, p)

#############################################
# starting simulation study
#############################################
set.seed(1234)

X <- Ximp <- Theta.x <- Sigma.x <- B <- mu <- E <- 
  Y <- Yimp <- Theta.y <- Sigma.y <- Z <- Zimp <- ZX <- vector(mode = "list", length = K)
unconnectedTheta <- list(NULL, alpha.true)
unconnectedOmg <- list(NULL, alpha.true)

indexB <- sapply(seq_len(p), \(i) sort(sample(q, size = 2L, replace = FALSE)))
unconnectedB <- list(NULL, alpha.true)

for(k in 1:K) {
  Theta.y[[k]] <- simulTheta(p = p, unconnected = unconnectedTheta[[k]],
                             ncomp = ncomp, mintht = Thtmin, maxtht = Thtmax)
  Sigma.y[[k]] <- solve(Theta.y[[k]])
  
  Theta.x[[k]] <- simulOmega(p = q, unconnected = unconnectedTheta[[k]],
                             ncomp = ncomp, mintht = Omgmin, maxtht = Omgmax)
  Sigma.x[[k]] <- solve(Theta.x[[k]])
  
  X[[k]] <- mvrnorm(n = n, mu = rep(0.0, q), Sigma = Sigma.x[[k]])
  B[[k]] <- simulB(X[[k]], perc, up, q, p, Bmin, Bmax, indexB, unconnectedB[[k]])
}

# pdf("~/Downloads/simul_new/structureK_Example.pdf", width = 16*.75, height = 10*.75)
par(mfrow = c(2, 3))
for(k in seq_len(K)) {
  image(Theta.y[[k]] != 0, col = c("white", "black"), 
        ylab = "", cex.lab = 2, cex.main = 3, 
        main = ifelse(k == 1, expression(Theta), ""), axes = FALSE)
  mtext(text = paste0("k = ", k), side = 2, line = 2.5, cex = 1.5)
  axis(1, at = c(0, .2, .4, .6, .8, 1), labels = c(0, 10, 20, 30, 40, 50), cex.axis = 2)
  axis(2, at = c(0, .2, .4, .6, .8, 1), labels = c(0, 10, 20, 30, 40, 50), cex.axis = 2)
  image(Theta.x[[k]] != 0, col = c("white", "black"), 
        ylab = "", cex.lab = 2, cex.main = 3, 
        main = ifelse(k == 1, expression(Omega), ""), axes = FALSE)
  axis(1, at = c(0, .2, .4, .6, .8, 1), labels = c(0, 5, 10, 15, 20, 25), cex.axis = 2)
  axis(2, at = c(0, .2, .4, .6, .8, 1), labels = c(0, 5, 10, 15, 20, 25), cex.axis = 2)
  image(B[[k]] != 0, col = c("white", "black"), 
        ylab = "", cex.lab = 2, cex.main = 3, 
        main = ifelse(k == 1, expression(B), ""), axes = FALSE)
  axis(1, at = c(0, .2, .4, .6, .8, 1), labels = c(0, 5, 10, 15, 20, 25), cex.axis = 2)
  axis(2, at = c(0, .2, .4, .6, .8, 1), labels = c(0, 10, 20, 30, 40, 50), cex.axis = 2)
}
# dev.off()

U <- outer(1:p, 1:p, "<=")
U2 <- outer(1:p, 1:p, "<")
U.2 <- outer(1:q, 1:q, "<=")
U2.2 <- outer(1:q, 1:q, "<")

# 1 - sum((B[[2]][-1L, ] != 0) & (B[[1]][-1L, ] != 0)) / sum(B[[1]][-1L, ] != 0)
# 1 - sum((Theta.y[[2]][U2] != 0) & (Theta.y[[1]][U2] != 0)) / sum(Theta.y[[1]][U2] != 0)
# 1 - sum((Theta.x[[2]][U2.2] != 0) & (Theta.x[[1]][U2.2] != 0)) / sum(Theta.x[[1]][U2.2] != 0)

# Output objects
MSE_JMMLE <- AUC_JMMLE <- array(0.0, dim = c(n.sim, 3L, K + 1L), 
                                dimnames = list(sim = seq_len(n.sim), 
                                                statistic = c("B", "Theta", "Omega"),
                                                class.id = c(seq_len(K), "mean")))

MSE_CGLASSO <- AUC_CGLASSO <- array(0.0, dim = c(n.sim, 3L, K + 1L), 
                                    dimnames = list(sim = seq_len(n.sim), 
                                                    statistic = c("B", "Theta", "Omega"),
                                                    class.id = c(seq_len(K), "mean")))

MSE_JCGLASSO <- AUC_JCGLASSO <- array(0.0, dim = c(length(alpha.seq), n.sim, 3L, K + 1L), 
                                      dimnames = list(alpha = alpha.seq, 
                                                      sim = seq_len(n.sim), 
                                                      statistic = c("B", "Theta", "Omega"),
                                                      class.id = c(seq_len(K), "mean")))

# Main simulation
for(i in 1:n.sim) {
  for(k in 1:K) {
    X[[k]] <- mvrnorm(n = n, mu = rep(0.0, q), Sigma = Sigma.x[[k]])
    B[[k]] <- simulB(X[[k]], perc, up, q, p, Bmin, Bmax, indexB, unconnectedB[[k]])
    eta <- cbind(1.0, X[[k]]) %*% B[[k]]
    mu[[k]] <- eta
    
    X[[k]][sample(n, floor(n * perc.NA)), sample(q, floor(q * perc.X.na))] <- NA
    Ximp[[k]] <- impute(X[[k]])
    
    E[[k]] <- mvrnorm(n = n, mu = rep(0, p), Sigma = Sigma.y[[k]])
    Y[[k]] <- mu[[k]] + E[[k]]
    
    Y[[k]][sample(n, floor(n * perc.NA)), sample(p, floor(p * perc.Y.na))] <- NA
    Yimp[[k]] <- impute(Y[[k]])
  }
  ZX <- datajcggm(Y = X) 
  Z <- datajcggm(Y = Y, X = X) 
  Zimp <- datajcggm(Y = Yimp, X = Ximp)
  
  ##### JMMLE SECTION #####
  tmpJESM <- JSEMcust(Ximp, Theta.X)
  tmpJMMLE <- JMMLEcust(Ximp, Theta.x, Yimp, Theta.y, B)
  for(k in 1:K) {
    MSE_JMMLE[i, "Omega", k] <- min(sapply(1:8, \(r) sum((tmpJESM[, , k, r] - Theta.x[[k]])[U.2]^2) / sum(Theta.x[[k]][U.2] != 0)))
    omegaJSEM <- PrecRecAUC(A = tmpJESM[, , k, ], B = Theta.x[[k]])
    AUC_JMMLE[i, "Omega", k] <- omegaJSEM$auc
    
    mse_b <- array(0, dim = c(8L, 7L), dimnames = list(lambda = 1:8, rho = 1:7))
    for(l in 1:8)
      for(r in 1:7)
        mse_b[l, r] <- sum((tmpJMMLE$B[, , k, l, r] - B[[k]])^2) / sum(B[[k]] != 0)
    MSE_JMMLE[i, "B", k] <- min(mse_b)
    id.b <- which(mse_b == min(mse_b), arr.ind = T)
    AUC_JMMLE[i, "B", k] <- sapply(1:7, \(r) PrecRecAUC(A = tmpJMMLE$B[, , k, , r], B = B[[k]])$auc)[id.b[1, 2]]
    
    mse_tht <- array(0, dim = c(8L, 7L), dimnames = list(lambda = 1:8, rho = 1:7))
    for(l in 1:8)
      for(r in 1:7)
        mse_tht[l, r] <- sum((tmpJMMLE$Theta[, , k, l, r] - Theta.y[[k]])[U]^2) / sum(Theta.y[[k]][U] != 0)
    MSE_JMMLE[i, "Theta", k] <- min(mse_tht)
    id.tht <- which(mse_tht == min(mse_tht), arr.ind = T)
    AUC_JMMLE[i, "Theta", k] <- sapply(1:8, \(l) PrecRecAUC(A = tmpJMMLE$Theta[, , k, l, ], B = Theta.y[[k]])$auc)[id.tht[1, 1]]
  }
  AUC_JMMLE[i, , K + 1] <- rowMeans(AUC_JMMLE[i, , 1:K])
  MSE_JMMLE[i, , K + 1] <- rowMeans(MSE_JMMLE[i, , 1:K])
  
  ##### CGLASSO SECTION #####
  for(k in 1:K){
    tmpCGLASSO <- cglasso(data = datacggm(Y = X[[k]]) , nrho = 10L, rho.min.ratio = 0.01,
                          trace = 0L, thr.bcd = 1E-5, thr.em = 1E-4)  
    MSE_CGLASSO[i, "Omega", k] <- min(sapply(1:10L, \(r) 
                                             sum((coef(tmpCGLASSO, type = "Theta", lambda.id = 1L, rho.id = r) - Theta.x[[k]])[U.2]^2) / sum(Theta.x[[k]][U.2] != 0)))
    omegaCGLASSO <- PrecRecAUC(A = tmpCGLASSO$Tht[, , 1L, ], B = Theta.x[[k]])
    AUC_CGLASSO[i, "Omega", k] <- omegaCGLASSO$auc
    
    id.min <- which.min(sapply(1:10L, \(r) 
                               sum((coef(tmpCGLASSO, type = "Theta", lambda.id = 1L, rho.id = r) - Theta.x[[k]])[U.2]^2) / sum(Theta.x[[k]][U.2] != 0)))
    tmpXimp <- tmpCGLASSO$Yipt[, , 1L, id.min]
    colnames(tmpXimp) <- paste0("X", 1:q)
    tmpZ <- datacggm(Y = Y[[k]], X = tmpXimp)
    
    # mse_B <- array(0.0, dim = c(4L, 10L), dimnames = list(rho = c2, lambda = c3))
    # auc_B <- double(4L)
    # tmpCGLASSO2 <- cglasso(data = tmpZ, nlambda = 1L, nrho = 4L, rho.min.ratio = .25, thr.bcd = 1E-5, thr.em = 1E-4)
    # 
    # for(r in 1:4L) {
    #   kkt <- abs(t(tmpXimp) %*% tmpCGLASSO2$R[, , 1L, r] %*% tmpCGLASSO2$Tht[, , 1L, r]) / matrix(diag(tmpCGLASSO2$Tht[, , 1L, r]), nrow = q, ncol = p, byrow = TRUE)
    #   lmb_max <- max(kkt) / n
    #   lmb <- lmb_max * c3
    #   tmpCGLASSO2.2 <- cglasso(data = tmpZ, rho = tmpCGLASSO2$rho[r], lambda = lmb, thr.em = 1e-4, thr.bcd = 1e-5, trace = 0L)
    #   mse_B[r, ] <-  apply(tmpCGLASSO2.2$B[, , , 1L], 3L, \(Bh) sum((Bh - B[[k]])^2) / sum(B[[k]] != 0))
    #   auc_B[r] <- PrecRecAUC(A = tmpCGLASSO2.2$B[, , , 1L], B = B[[k]])$auc
    # }
    # id.b <- which(mse_B == min(mse_B), arr.ind = TRUE)
    # MSE_CGLASSO[i, "B", k] <- mse_B[id.b[1, 1], id.b[1, 2]]
    # AUC_CGLASSO[i, "B", k] <- auc_B[id.b[1, 1]]
    # 
    # mse_tht <- array(0.0, dim = c(4L, 10L), dimnames = list(lambda = c2, rho = c3))
    # auc_tht <- double(4L)
    # tmpCGLASSO3 <- cglasso(data = tmpZ, nlambda = 4L, lambda.min.ratio = .25, nrho = 1L, thr.bcd = 1E-5, thr.em = 1E-4)
    # 
    # for(l in 1:4L) {
    #   rho_max <- max(abs(tmpCGLASSO3$S[, , l, 1L][U]))
    #   rho <- rho_max * c3
    #   tmpCGLASSO3.2 <- cglasso(data = tmpZ, lambda = tmpCGLASSO3$lambda[l], rho = rho, thr.em = 1e-4, thr.bcd = 1e-5, trace = 0L)
    #   mse_tht[l, ] <-  apply(tmpCGLASSO3.2$Tht[, , 1L, ], 3L, \(Tht) sum((Tht[U] - Theta.y[[k]][U])^2) / sum(Theta.y[[k]][U] != 0))
    #   auc_tht[l] <- PrecRecAUC(A = tmpCGLASSO3.2$Tht[, , 1L, ], B = Theta.y[[k]])$auc
    # }
    # id.tht <- which(mse_tht == min(mse_tht), arr.ind = TRUE)
    # MSE_CGLASSO[i, "Theta", k] <- mse_B[id.b[1, 1], id.b[1, 2]]
    # AUC_CGLASSO[i, "Theta", k] <- auc_B[id.b[1, 1]]
    
    tmpCGLASSOmax <- cglasso(data = tmpZ, nlambda = 10L, lambda.min.ratio = .1,
                             nrho = 10L, rho.min.ratio = .1, thr.bcd = 1E-5, thr.em = 1E-4)
    
    mse_b <- array(0, dim = c(10L, 10L), dimnames = list(lambda = 1:10, rho = 1:10))
    for(l in 1:10)
      for(r in 1:10)
        mse_b[l, r] <- sum((tmpCGLASSOmax$B[, , l, r] - B[[k]])^2) / sum(B[[k]] != 0)
    MSE_CGLASSO[i, "B", k] <- min(mse_b)
    id.b <- which(mse_b == min(mse_b), arr.ind = T)
    AUC_CGLASSO[i, "B", k] <- sapply(1:10, \(r) PrecRecAUC(A = tmpCGLASSOmax$B[, , , r], B = B[[k]])$auc)[id.b[1, 2]]
    
    mse_tht <- array(0, dim = c(10L, 10L), dimnames = list(lambda = 1:10, rho = 1:10))
    for(l in 1:10)
      for(r in 1:10)
        mse_tht[l, r] <- sum((tmpCGLASSOmax$Tht[, , l, r] - Theta.y[[k]])[U]^2) / sum(Theta.y[[k]][U] != 0)
    MSE_CGLASSO[i, "Theta", k] <- min(mse_tht)
    id.tht <- which(mse_tht == min(mse_tht), arr.ind = T)
    AUC_CGLASSO[i, "Theta", k] <- sapply(1:10, \(l) PrecRecAUC(A = tmpCGLASSOmax$Tht[, , l, ], B = Theta.y[[k]])$auc)[id.tht[1, 1]]
  }
  AUC_CGLASSO[i, , K + 1] <- rowMeans(AUC_CGLASSO[i, , 1:K])
  MSE_CGLASSO[i, , K + 1] <- rowMeans(MSE_CGLASSO[i, , 1:K])
  
  ##### JCGLASSO SECTION #####
  # Simulation varying alpha
  for(a in 1:length(alpha.seq)) {
    tmpX <- jcglasso(data = ZX, nrho = 10L, rho.min.ratio = 0.01, 
                     trace = 0L, thr.bcd = 1E-5, thr.em = 1E-3, alpha1 = alpha.seq[a])
    tmpMSE_JCGLASSO <- rowSums(sapply(1:K, \(k) 
                                      sapply(1:10L, \(r) 
                                             sum((coef(tmpX, type = "Theta", lambda.id = 1L, rho.id = r, class.id = k) - Theta.x[[k]])[U.2]^2) / sum(Theta.x[[k]][U.2] != 0))))
    AUC_JCGLASSO[a, i, "Omega", 1:K] <- sapply(1:K, \(k) 
                                               PrecRecAUC(A = tmpX$Tht[, , k, 1L, ], B = Theta.x[[k]])$auc)
    
    nu.opt <- tmpX$rho[which.min(tmpMSE_JCGLASSO)]
    
    # mse_omg <- array(0.0, dim = c(4L, 10L, 2L), dimnames = list(rho = c2, lambda = c3, K = 1:K))
    # mse_B <- array(0.0, dim = c(4L, 10L, 2L), dimnames = list(rho = c2, lambda = c3, K = 1:K))
    # auc_B <- matrix(0.0, nrow = 4L, ncol = 2L, dimnames = list(rho = c2, K = 1:K))
    # tmp2 <- jcglasso(data = Z, nlambda = 1L, nrho = 4L, rho.min.ratio = 0.25, nu = nu.opt,
    #                  trace = 0L, thr.bcd = 1E-4, thr.em = 1E-3, 
    #                  alpha1 = alpha.seq[a], alpha2 = alpha.seq[a], alpha3 = alpha.seq[a])
    # for(r in 1:length(c2)) {
    #   lMax2 <- lmbMax_fun(tmp2, lambda.id = 1L, rho.id = r)
    #   lambda.seq <- lMax2 * c3
    #   tmp2.2 <- jcglasso(data = Z, lambda = lambda.seq, rho = tmp2$rho[r], nu = nu.opt,
    #                      trace = 0L, thr.bcd = 1E-4, thr.em = 1E-3, 
    #                      alpha1 = alpha.seq[a], alpha2 = alpha.seq[a], alpha3 = alpha.seq[a])
    #   
    #   mse_omg[r, , ] <- sapply(seq_len(K), \(k) apply(tmp2.2$Omega[, , k, , 1L], 3, function(Omg) sum((Omg[U.2] - Theta.x[[k]][U.2])^2) / sum(Theta.x[[k]][U.2] != 0)))
    #   mse_B[r, , ] <- sapply(seq_len(K), \(k) apply(tmp2.2$B[, , k, , 1L], 3, function(Bh) sum((Bh - B[[k]])^2) / sum(B[[k]] != 0)))
    #   auc_B[r, ] <- sapply(1:K, \(k) PrecRecAUC(A = tmp2.2$B[, , k, , 1L], B = B[[k]])$auc)
    # }
    # id.b <- which(apply(mse_B, 1, rowSums) == min(apply(mse_B, 1, rowSums)), arr.ind = TRUE)
    # MSE_JCGLASSO[a, i, "Omega", 1:K] <- mse_omg[id.b[1, 2], id.b[1, 1], ]
    # MSE_JCGLASSO[a, i, "B", 1:K] <- mse_B[id.b[1, 2], id.b[1, 1], ]
    # AUC_JCGLASSO[a, i, "B", 1:K] <- auc_B[id.b[1, 2],]
    # 
    # mse_omg2 <- array(0.0, dim = c(4L, 10L, 2L), dimnames = list(lambda = c2, rho = c3, K = 1:K))
    # mse_tht <- array(0.0, dim = c(4L, 10L, 2L), dimnames = list(lambda = c2, rho = c3, K = 1:K))
    # auc_tht <- matrix(0.0, nrow = 4L, ncol = 2L, dimnames = list(lambda = c2, K = 1:K))
    # tmp3 <- jcglasso(data = Z, nlambda = 4L, nrho = 1L, lambda.min.ratio = 0.25, nu = nu.opt,
    #                  trace = 0L, thr.bcd = 1E-4, thr.em = 1E-3, 
    #                  alpha1 = alpha.seq[a], alpha2 = alpha.seq[a], alpha3 = alpha.seq[a])
    # 
    # for(l in 1:length(c2)) {
    #   rMax2 <- rhoMax_fun(tmp3, lambda.id = l, rho.id = 1L)
    #   # rMax2 <- max(abs(tmp3$S[, , l, 1L][outer(1:p, 1:p, "<")]))
    #   rho.seq <- rMax2 * c3
    #   tmp3.2 <- jcglasso(data = Z, lambda = tmp3$lambda[l], rho = rho.seq, nu = nu.opt,
    #                      trace = 0L, thr.bcd = 1E-4, thr.em = 1E-3, 
    #                      alpha1 = alpha.seq[a], alpha2 = alpha.seq[a], alpha3 = alpha.seq[a])
    #   
    #   mse_omg2[l, , ] <- sapply(seq_len(K), \(k) apply(tmp3.2$Omega[, , k, 1L, ], 3, function(Omg) sum((Omg[U.2] - Theta.x[[k]][U.2])^2) / sum(Theta.x[[k]][U.2] != 0)))
    #   mse_tht[l, , ] <-  sapply(seq_len(K), \(k) apply(tmp3.2$Tht[1:p, 1:p, k, 1L, ], 3L, function(Th) sum((Th[U] - Theta.y[[k]][U])^2) / sum(Theta.y[[k]][U] != 0)))
    #   auc_tht[l, ] <- sapply(1:K, \(k) PrecRecAUC(A = tmp3.2$Tht[1:p, 1:p, k, 1L, ], B = Theta.y[[k]])$auc)
    # }
    # id.tht <- which(apply(mse_tht, 1, rowSums) == min(apply(mse_tht, 1, rowSums)), arr.ind = TRUE)
    # id.omg <- which.min(rowSums(rbind(MSE_JCGLASSO[a, i, "Omega", 1:K], mse_omg2[id.tht[1, 2], id.tht[1, 1], ])))
    # MSE_JCGLASSO[a, i, "Omega", 1:K] <- rbind(MSE_JCGLASSO[a, i, "Omega", 1:K], mse_omg2[id.tht[1, 2], id.tht[1, 1], ])[id.omg, ]
    # MSE_JCGLASSO[a, i, "Theta", 1:K] <- mse_tht[id.tht[1, 2], id.tht[1, 1], ]
    # AUC_JCGLASSO[a, i, "Theta", 1:K] <- auc_tht[id.tht[1, 2],]
    
    tmp2 <- jcglasso(data = Z, nlambda = 10L, lambda.min.ratio = .1,
                     nrho = 10L, rho.min.ratio = .1, nu = nu.opt,
                     trace = 0L, thr.bcd = 1E-4, thr.em = 1E-3,
                     alpha1 = alpha.seq[a], alpha2 = alpha.seq[a], alpha3 = alpha.seq[a])

    mse_B <- array(0, dim = c(10L, 10L), dimnames = list(lambda = 1:10, rho = 1:10))
    for(l in 1:10)
      for(r in 1:10)
        mse_B[l, r] <- sum(sapply(seq_len(K), \(k) sum((tmp2$B[, , k, l, r] - B[[k]])^2) / sum(B[[k]] != 0)))
    id.b <- which(mse_B == min(mse_B), arr.ind = TRUE)
    MSE_JCGLASSO[a, i, "B", 1:K] <- sapply(seq_len(K), \(k) sum((tmp2$B[, , k, id.b[1, 1], id.b[1, 2]] - B[[k]])^2) / sum(B[[k]] != 0))
    AUC_JCGLASSO[a, i, "B", 1:K] <- sapply(1:10L, \(r) sapply(1:K, \(k) PrecRecAUC(A = tmp2$B[, , k, , r], B = B[[k]])$auc))[, id.b[1, 2]]

    mse_tht <- array(0, dim = c(10L, 10L), dimnames = list(lambda = 1:10, rho = 1:10))
    for(l in 1:10)
      for(r in 1:10)
        mse_tht[l, r] <- sum(sapply(seq_len(K), \(k) sum((tmp2$Tht[1:p, 1:p, k, l, r] - Theta.y[[k]])[U]^2) / sum(Theta.y[[k]][U] != 0)))
    id.tht <- which(mse_tht == min(mse_tht), arr.ind = TRUE)
    MSE_JCGLASSO[a, i, "Theta", 1:K] <- sapply(seq_len(K), \(k) sum((tmp2$Tht[1:p, 1:p, k, id.tht[1, 1], id.tht[1, 2]] - Theta.y[[k]])[U]^2) / sum(Theta.y[[k]][U] != 0))
    AUC_JCGLASSO[a, i, "Theta", 1:K] <- sapply(1:10L, \(l) sapply(1:K, \(k) PrecRecAUC(A = tmp2$Tht[1:p, 1:p, k, l, ], B = Theta.y[[k]])$auc))[, id.tht[1, 1]]

    mse_omg <- array(0, dim = c(10L, 10L), dimnames = list(lambda = 1:10, rho = 1:10))
    for(l in 1:10)
      for(r in 1:10)
        mse_omg[l, r] <- sum(sapply(seq_len(K), \(k) sum((tmp2$Omega[, , k, l, r] - Theta.x[[k]])[U.2]^2) / sum(Theta.x[[k]][U.2] != 0)))
    id.omg <- which(mse_omg == min(mse_omg), arr.ind = TRUE)
    MSE_JCGLASSO[a, i, "Omega", 1:K] <- sapply(seq_len(K), \(k) sum((tmp2$Omega[, , k, id.omg[1, 1], id.omg[1, 2]] - Theta.x[[k]])[U.2]^2) / sum(Theta.x[[k]][U.2] != 0))
    
    AUC_JCGLASSO[a, i, , K + 1] <- rowMeans(AUC_JCGLASSO[a, i, , 1:K])
    MSE_JCGLASSO[a, i, , K + 1] <- rowMeans(MSE_JCGLASSO[a, i, , 1:K])
  }
  
  ##### PRINTING SECTION #####
  cat("iter n. ", i)
  
  cat("\nMSE JMMLE\n")
  print(MSE_JMMLE[i, , K + 1])
  cat("\nMSE of CGLASSO\n")
  print(MSE_CGLASSO[i, , K + 1])
  cat("\nMSE of JCGLASSO")
  print(MSE_JCGLASSO[, i, , K + 1])
  
  cat("\nAUC JMMLE\n")
  print(AUC_JMMLE[i, , K + 1])
  cat("\nAUC of CGLASSO\n")
  print(AUC_CGLASSO[i, , K + 1])
  cat("\nAUC of JCGLASSO")
  print(AUC_JCGLASSO[, i, , K + 1])
  
  cat("\n")
  
  # if(i %% 5 == 0) save.image(paste0("~/Downloads/simul_new/SimulJRSSC_31102023_n", n, "p", p, "q", q,
  #                                   "NA", perc.NA*100, "Xna", perc.X.na*100, "Yna", perc.Y.na*100,
  #                                   "alpha", alpha.true, "_comparison.RData"))
}

MSE1 <- sapply(1:K, \(k) colMeans(MSE_JMMLE[, , k]))
MSE2 <- sapply(1:K, \(k) colMeans(MSE_CGLASSO[, , k]))

MSE3 <- lapply(1:length(alpha.seq), \(a) 
               sapply(1:K, \(k) colMeans(MSE_JCGLASSO[a, , , k])))

AUC1 <- sapply(1:K, \(k) colMeans(AUC_JMMLE[, , k]))
AUC2 <- sapply(1:K, \(k) colMeans(AUC_CGLASSO[, , k]))

AUC3 <- lapply(1:length(alpha.seq), \(a) 
               sapply(1:K, \(k) colMeans(AUC_JCGLASSO[a, , , k])))

round(cbind(t(MSE1), t(MSE2), t(MSE3[[1]]), t(MSE3[[2]]), t(MSE3[[3]])), 3L)
round(cbind(t(AUC1), t(AUC2), t(AUC3[[1]]), t(AUC3[[2]]), t(AUC3[[3]])), 3L)

# clip <- pipe("pbcopy", "w")
# write.table(round(cbind(t(MSE1), t(MSE2), t(MSE3[[1]]), t(MSE3[[2]]), t(MSE3[[3]])), 3L), file = clip, sep = "\t", row.names = TRUE)
# close(clip)
# 
# clip <- pipe("pbcopy", "w")
# write.table(round(cbind(t(AUC1), t(AUC2), t(AUC3[[1]]), t(AUC3[[2]]), t(AUC3[[3]])), 3L), file = clip, sep = "\t", row.names = TRUE)
# close(clip)

