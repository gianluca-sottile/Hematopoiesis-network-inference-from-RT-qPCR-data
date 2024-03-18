rm(list = ls())

# libraries
library(MASS)
library(cglasso)
library(JGL)
library(huge)
library(ggplot2)
library(caret)
library(BDgraph)

source("../../01 - RCode/datajcggm.R")
source("../../01 - RCode/jcglasso.R")
source("../../01 - RCode/jcggm.R")
source("../../01 - RCode/gof.R")
source("../../01 - RCode/to_graph.R")

source('JMMLE/JMLE.R')
source('JMMLE/l1LS_Main.R')
source('JMMLE/jsem.R')
source('JMMLE/Objval.R')
require(glmnet)
require(grpreg)
require(glasso)
require(parallel)

simulTheta <- function(n, p, prob, alpha, seed = 1234, vis = TRUE, graph = c("random", "smallworld", "star", "circle", "cluster", "scale-free", "lattice", "hub", "AR(1)", "AR(2)"), ...){
  graph <- match.arg(graph)
  # full.size <- p * (p-1) / 2
  # size <- round(alpha * prob * full.size)
  set.seed(1234)
  data.sim <- bdgraph.sim(p = p, n = n, prob = prob * alpha, vis = vis, graph = graph, ...)
  data.sim
}
simulOmega <- function(n, q, prob, alpha, seed = 1234, vis = TRUE, graph = c("random", "smallworld", "star", "circle", "cluster", "scale-free", "lattice", "hub", "AR(1)", "AR(2)"), ...){
  graph <- match.arg(graph)
  # full.size <- q * (q-1) / 2
  # size <- round(alpha * prob * full.size)
  set.seed(1234)
  data.sim <- bdgraph.sim(p = q, n = n, prob = prob * alpha, vis = vis, graph = graph, ...)
  data.sim
}
simulB <- function(X, perc, up, q, p, Bmin = .5, Bmax = .7, indexB, unconnectedB) {
  B <- matrix(0.0, nrow = q + 1L, ncol = p)
  for(j in seq_len(p)) B[1L + indexB[, j], j] <- runif(nrow(indexB), Bmin, Bmax)
  
  eta <- X %*% B[-1L, ]
  if(up > 0) {
    cutoff <- function(b0, i, up, perc = perc) (mean(pnorm(up - b0 - eta[, i], lower.tail =  FALSE))) - perc
    B[1L, ] <- sapply(1:p, function(.i) uniroot(cutoff, interval = c(0, 100), i = .i, up = up, perc = perc[.i])$root)
  }
  
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

PrecRecAUC <- function(A, B, plot.it = FALSE, unpenalized = TRUE) {
  dimA <- dim(A)
  n <- dimA[3]
  p <- dimA[2]
  q <- dimA[1]
  if(unpenalized) 
    U <- if(q == p) outer(1:p, 1:p, "<=") else matrix(TRUE, q, p)
  else U <- if(q == p) outer(1:p, 1:p, "<") else rbind(rep(FALSE, p), matrix(TRUE, q-1, p))
  
  tmp <- sapply(1:n, \(i) {
    data <- factor(1*c(A[, , i][U] != 0), levels = c(1,0))
    reference <- factor(1 * c(B[U] != 0), levels = c(1,0))
    confusionMatrix(data = data, reference = reference, mode = "prec_recall")$byClass[c("Precision", "Recall")]
  })
  tmp <- cbind(c(1, 0), 
               if(unpenalized) c(1, p / sum(B[U] != 0)) else c(1, 0), 
               tmp, 
               c(sum(B[U] != 0) / sum(U), 1),
               c(0, 1))
  splFun <- suppressWarnings(splinefun(tmp["Recall", ], tmp["Precision", ], method = "monoH.FC"))
  # splFun <- approxfun(tmp["Recall", ], tmp["Precision", ], rule = 2)
  if(plot.it) {
    plot(tmp["Recall", ], tmp["Precision", ], type = "b", xlim = c(0, 1), ylim = c(0, 1))
    curve(splFun, from = 0, to = 1, add = TRUE, col = 2, lty = 2)
    # x.seq <- seq(max(tmp["Recall", ]), min(tmp["Recall", ]), length = 100L)
    # lines(x.seq, splFun(x.seq), lty = 2, col = "red")
  }
  auc <- integrate(splFun, lower = 0, upper = 1)$value
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
    # lambda.vec <- sqrt(log(p)/n) * seq(1.8, 0.4, -0.2)
    lambda.vec <- sqrt(log(p)/n) * exp(seq(log(2), log(.1), l = 8))
  }
  nlambda <- length(lambda.vec)
  
  Y2 <- lapply(Y, scale, center = TRUE, scale = FALSE)
  
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
    if(length(Info) == 0) Info <- NULL
    Theta_refit <- multi.glasso(do.call(rbind, Y2), unlist(Y.indices), lambda.vec[l], Info)
    list2array(Theta_refit$Theta)
  })
  stopCluster(cl)
  for(l in 1:nlambda) Theta.array[, , , l] <- tmp[[l]]
  
  # for(l in 1:nlambda) {
  #   init.jsem.model <- JSEM(do.call(rbind, Y2), unlist(Y.indices),
  #                           Y.groups, lambda = lambda.vec[l])
  #   Ahat <- init.jsem.model$Ahat
  #   Info <- list()
  #   for (k in 1:K) Info[[k]] <- zeroInd(Ahat[[k]], 1)$zeroArr
  #   if(length(Info) == 0) Info <- NULL
  #   Theta_refit <- multi.glasso(do.call(rbind, Y2), unlist(Y.indices), 1E-8, Info)
  #   Theta.array[, , , l] <- list2array(Theta_refit$Theta)
  # }
  
  return(Theta.array)
}
JMMLEcust <- function(X, Theta.X, Y, Theta.Y, B0, lambda.vec = NULL, gamma.vec = NULL) {
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
    # lambda.vec <- sqrt(log(q)/n) * seq(1.8, 0.4, -0.2)
    lambda.vec <- sqrt(log(q)/n) * exp(seq(log(2), log(.1), l = 8))
  }
  if(is.null(gamma.vec)){
    # gamma.vec <- sqrt(log(p)/n) * seq(1, 0.4, -0.1)
    gamma.vec <- sqrt(log(p)/n) * exp(seq(log(2), log(.1), l = 7))
  }
  nlambda <- length(lambda.vec)
  ngamma <- length(gamma.vec)
  model.list<- vector("list", nlambda)
  
  X2 <- lapply(X, scale, center = TRUE, scale = FALSE)
  Y2 <- lapply(Y, scale, center = TRUE, scale = FALSE)
  
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
                  init.option = 1, tol = 1e-3, VERBOSE = FALSE, refit.B = FALSE)
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

n <- 100L              # sample size
p <- 200L              # number of response variables
q <- 50L               # number of predictors
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
unconnectedTheta <- list(1.0, alpha.true)
unconnectedOmg <- list(1.0, alpha.true)

indexB <- sapply(seq_len(p), \(i) sort(sample(q, size = 2L, replace = FALSE)))
unconnectedB <- list(NULL, alpha.true)

for(k in 1:K) {
  tempTheta <- simulTheta(n, p, prob = .005, alpha = unconnectedTheta[[k]], seed = 1234, # prob = c(0.005, 0.02)
                          vis = FALSE, graph = "cluster", rewire = .5, class = 1)
  Theta.y[[k]] <- round(tempTheta$K, 5L)
  Sigma.y[[k]] <- tempTheta$sigma
  
  tempOmega <- simulOmega(n, q, prob = .02, alpha = unconnectedOmg[[k]], seed = 1234, # prob = c(0.02, 0.08)
                          vis = FALSE, graph = "cluster", rewire = .5, class = 1)
  Theta.x[[k]] <- round(tempOmega$K, 5)
  Sigma.x[[k]] <- tempOmega$sigma
  
  X[[k]] <- tempOmega$data
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

summary(Theta.y[[1]][U2][Theta.y[[1]][U2]!=0])
summary(Theta.y[[2]][U2][Theta.y[[2]][U2]!=0])

mean(Theta.y[[1]][U]!=0)
mean(Theta.y[[2]][U]!=0)

summary(Theta.x[[1]][U2.2][Theta.x[[1]][U2.2]!=0])
summary(Theta.x[[2]][U2.2][Theta.x[[2]][U2.2]!=0])

mean(Theta.x[[1]][U.2]!=0)
mean(Theta.x[[2]][U.2]!=0)

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
  tmpJESM <- JSEMcust(Ximp, Theta.x, lambda.vec = sqrt(log(q) / n) * exp(seq(log(2), log(.1), l = 5L)))
  
  tmpMSE <- sapply(1:K, \(k) sapply(1:5, \(r) norm(tmpJESM[, , k, r] - Theta.x[[k]], "F")))
  # plot(sqrt(log(q)/n) * exp(seq(log(1), log(.1), l = 5)), rowSums(tmpMSE), type = "b")
  
  MSE_JMMLE[i, "Omega", 1:K] <- tmpMSE[which.min(rowSums(tmpMSE)), ]
  AUC_JMMLE[i, "Omega", 1:K] <- sapply(1:K, \(k) PrecRecAUC(A = tmpJESM[, , k, ], B = Theta.x[[k]])$auc)
  
  tmpJMMLE <- JMMLEcust(Ximp, Theta.x, Yimp, Theta.y, B, lambda.vec = sqrt(log(q)/n) * exp(seq(log(2), log(.1), l = 5)), gamma.vec = 1E6)
  tmpMSE <- sapply(1:K, \(k) sapply(1:5, \(l) norm(tmpJMMLE$B[, , k, l, 1L] - B[[k]], "F")))
  # plot(sqrt(log(q)/n) * exp(seq(log(1), log(.1), l = 5)), rowSums(tmpMSE), type="b")
  
  MSE_JMMLE[i, "B", 1:K] <- tmpMSE[which.min(rowSums(tmpMSE)), ]
  AUC_JMMLE[i, "B", 1:K] <- sapply(1:K, \(k) PrecRecAUC(A = tmpJMMLE$B[, , k, , 1L], B = B[[k]], unpenalized = (up!=0))$auc)
  
  lmb.opt <- sqrt(log(q)/n) * exp(seq(log(1), log(.1), l = 5))[which.min(rowSums(tmpMSE))]
  
  tmpJMMLE2 <- JMMLEcust(Ximp, Theta.x, Yimp, Theta.y, B, lambda.vec = lmb.opt, gamma.vec = sqrt(log(p)/n) * exp(seq(log(2), log(.1), l = 5)))
  tmpMSE <- sapply(1:K, \(k) sapply(1:5, \(r) norm(tmpJMMLE2$Theta[, , k, 1L, r] - Theta.y[[k]], "F")))
  # plot(sqrt(log(p)/n) * exp(seq(log(1), log(.1), l = 5)), rowSums(tmpMSE), type="b")
  
  MSE_JMMLE[i, "Theta", 1:K] <- tmpMSE[which.min(rowSums(tmpMSE)), ]
  AUC_JMMLE[i, "Theta", 1:K] <- sapply(1:K, \(k) PrecRecAUC(A = tmpJMMLE2$Theta[, , k, 1L, ], B = Theta.y[[k]])$auc)
  
  rho.opt <- sqrt(log(p)/n) * exp(seq(log(2), log(.1), l = 8))[which.min(rowSums(tmpMSE))]
  
  AUC_JMMLE[i, , K + 1] <- rowMeans(AUC_JMMLE[i, , 1:K])
  MSE_JMMLE[i, , K + 1] <- rowMeans(MSE_JMMLE[i, , 1:K])
  
  ##### CGLASSO SECTION #####
  for(k in 1:K){
    tmpCGLASSO <- cglasso(data = datacggm(Y = X[[k]]), nrho = 10L, rho.min.ratio = .01,
                          trace = 0L, thr.bcd = 1E-5, thr.em = 1E-3)
    # tmpMSE <- sapply(1:10L, \(r) 
    #                   norm(cggm(tmpCGLASSO, mle = TRUE, lambda.id = 1L, rho.id = r)$Tht[,,1L,1L] - Theta.x[[k]], "F") / sum(U.2))
    tmpMSE <- sapply(1:10L, \(r)
                      norm(tmpCGLASSO$Tht[,,1L,r] - Theta.x[[k]], "F"))
    # plot(tmpCGLASSO$rho, tmpMSE, type="b")
    
    id.min <- which.min(tmpMSE)
    MSE_CGLASSO[i, "Omega", k] <- tmpMSE[id.min]
    AUC_CGLASSO[i, "Omega", k] <- PrecRecAUC(A = tmpCGLASSO$Tht[, , 1L, ], B = Theta.x[[k]])$auc
    
    tmpXimp <- cggm(tmpCGLASSO, mle = TRUE, lambda.id = 1L, rho.id = id.min)$Yipt[, , 1L, 1L]
    colnames(tmpXimp) <- paste0("X", 1:q)
    tmpZ <- datacggm(Y = Y[[k]], X = tmpXimp)
    
    tmpCGLASSOmax <- cglasso(data = tmpZ, nlambda = 1L, rho = 1E6, thr.bcd = 1E-5, thr.em = 1E-3, trace = 0)
    lmb.seq <- exp(seq(log(tmpCGLASSOmax$lambda), log(tmpCGLASSOmax$lambda*.01), l = 10))
    
    tmpCGLASSOmax2 <- cglasso(data = tmpZ, lambda = lmb.seq, rho = 1E6,
                               trace = 0L, thr.bcd = 1E-5, thr.em = 1E-3)
    # tmpMSE <- sapply(1:10, \(l)
    #                  norm(cggm(tmpCGLASSOmax2, lambda.id = l, rho.id = 1)$B[, , 1L, 1L] - B[[k]], "F") / prod(dim(B[[k]])))
    tmpMSE <- sapply(1:10, \(l)
                     norm(tmpCGLASSOmax2$B[, , l, 1L] - B[[k]], "F"))
    # plot(tmpCGLASSOmax2$lambda, tmpMSE, type="b")
    
    MSE_CGLASSO[i, "B", k] <- tmpMSE[which.min(tmpMSE)]
    AUC_CGLASSO[i, "B", k] <- PrecRecAUC(A = tmpCGLASSOmax2$B[, , , 1L], B = B[[k]], unpenalized = (up!=0))$auc
    
    lmb.opt <- tmpCGLASSOmax2$lambda[which.min(tmpMSE)]
    rmax <- max(abs((crossprod(tmpCGLASSOmax2$R[, , which.min(tmpMSE), 1L]) / n)[U2]))
    rho.seq <- exp(seq(log(rmax), log(rmax * .01), l = 10))
    
    tmpCGLASSOmax3 <- cglasso(data = tmpZ, lambda = lmb.opt, rho = rho.seq, 
                         trace = 0L, thr.bcd = 1E-5, thr.em = 1E-3)
    # tmpMSE <- sapply(1:10, \(r)
    #                  norm(cggm(tmpCGLASSOmax3, lambda.id = 1L, rho.id = r)$Tht[, , 1L, 1L] - Theta.y[[k]], "F") / sum(U))
    tmpMSE <- sapply(1:10, \(r)
                     norm(tmpCGLASSOmax3$Tht[, , 1L, r] - Theta.y[[k]], "F"))
    # plot(tmpCGLASSOmax3$rho, tmpMSE, type="b")
    
    MSE_CGLASSO[i, "Theta", k] <- tmpMSE[which.min(tmpMSE)]
    AUC_CGLASSO[i, "Theta", k] <- PrecRecAUC(A = tmpCGLASSOmax3$Tht[, , 1L, ], B = Theta.y[[k]])$auc
    
    rho.opt <- tmpCGLASSOmax3$rho[which.min(tmpMSE)]
    # lmax <- max(abs(crossprod(as.matrix(getMatrix(tmpZ, "X")), tmpCGLASSOmax3$R[, , 1L, which.min(tmpMSE)]) %*% 
    #   tmpCGLASSOmax3$Tht[, , 1L, which.min(tmpMSE)] / (n * matrix(diag(tmpCGLASSOmax3$Tht[, , 1L, which.min(tmpMSE)]), q, p))))
    # lmb.seq <- exp(seq(log(lmax), log(lmax * .10), l = 10))
  }
  AUC_CGLASSO[i, , K + 1] <- rowMeans(AUC_CGLASSO[i, , 1:K])
  MSE_CGLASSO[i, , K + 1] <- rowMeans(MSE_CGLASSO[i, , 1:K])
  
  ##### JCGLASSO SECTION #####
  # Simulation varying alpha
  for(a in 1:length(alpha.seq)) {
    tmpX <- jcglasso(data = ZX, nrho = 10L, rho.min.ratio = 0.01, 
                     trace = 0L, thr.bcd = 1E-4, thr.em = 1E-3, alpha1 = alpha.seq[a])
    # tmpMSE <- t(sapply(1:10, \(r) {
    #   temp <- jcggm(tmpX, lambda.id = 1L, rho.id = r)
    #   sapply(1:K, \(k) norm(temp$Tht[, , k, 1L, 1L] - Theta.x[[k]], "F") / sum(U.2))
    # }))
    tmpMSE <- sapply(1:K, \(k)
                     sapply(1:10, \(r)
                            norm(tmpX$Tht[, , k, 1L, r] - Theta.x[[k]], "F")))
    # plot(tmpX$rho, rowSums(tmpMSE), type="b")
    
    AUC_JCGLASSO[a, i, "Omega", 1:K] <- sapply(1:K, \(k) PrecRecAUC(A = tmpX$Tht[, , k, 1L, ], B = Theta.x[[k]])$auc)
    
    nu.opt <- tmpX$rho[which.min(rowSums(tmpMSE))]
    
    tmp2max <- jcglasso(data = Z, nlambda = 1L, rho = 1E6, nu = nu.opt,
                     trace = 0L, thr.bcd = 1E-4, thr.em = 1E-3, 
                     alpha1 = alpha.seq[a], alpha2 = alpha.seq[a], alpha3 = alpha.seq[a])
    
    lmb.seq <- exp(seq(log(tmp2max$lambda), log(tmp2max$lambda*.01), l = 10))
    
    tmp2max2 <- jcglasso(data = Z, lambda = lmb.seq, rho = 1E6, nu = nu.opt,
                        trace = 0L, thr.bcd = 1E-4, thr.em = 1E-3,
                        alpha1 = alpha.seq[a], alpha2 = alpha.seq[a], alpha3 = alpha.seq[a])
    # tmpMSE <- t(sapply(1:10, \(l) {
    #   temp <- jcggm(tmp2max2, lambda.id = l, rho.id = 1)
    #   sapply(1:K, \(k) norm(temp$B[, , k, 1L, 1L] - B[[k]], "F") / prod(dim(B[[k]])))
    # }))
    tmpMSE <- sapply(1:K, \(k)
                     sapply(1:10, \(l)
                            norm(tmp2max2$B[, , k, l, 1L] - B[[k]], "F")))
    # plot(tmp2max2$lambda, rowSums(tmpMSE), type="b")
    
    MSE_JCGLASSO[a, i, "B", 1:K] <- tmpMSE[which.min(rowSums(tmpMSE)), ]
    AUC_JCGLASSO[a, i, "B", 1:K] <- sapply(1:K, \(k) PrecRecAUC(A = tmp2max2$B[, , k, , 1L], B = B[[k]], unpenalized = (up!=0))$auc)
    
    lmb.opt <- tmp2max2$lambda[which.min(rowSums(tmpMSE))]
    rmax <- rhoMax_fun(tmp2max2, lambda.id = which.min(rowSums(tmpMSE)), rho.id = 1L)
    rho.seq <- exp(seq(log(rmax), log(rmax * .01), l = 10))
    
    tmp2max3 <- jcglasso(data = Z, lambda = lmb.opt, rho = rho.seq, nu = nu.opt,
                         trace = 0L, thr.bcd = 1E-4, thr.em = 1E-3,
                         alpha1 = alpha.seq[a], alpha2 = alpha.seq[a], alpha3 = alpha.seq[a])
    # tmpMSE <- t(sapply(1:10, \(r) {
    #   temp <- jcggm(tmp2max3, lambda.id = 1L, rho.id = r)
    #   sapply(1:K, \(k) norm(temp$Tht[1:p, 1:p, k, 1L, 1L] - Theta.y[[k]], "F") / sum(U))
    # }))
    tmpMSE <- sapply(1:K, \(k)
                     sapply(1:10, \(r)
                            norm(tmp2max3$Tht[1:p, 1:p, k, 1L, r] - Theta.y[[k]], "F")))
    # plot(tmp2max3$rho, rowSums(tmpMSE), type="b")
    
    MSE_JCGLASSO[a, i, "Theta", 1:K] <- tmpMSE[which.min(rowSums(tmpMSE)), ]
    AUC_JCGLASSO[a, i, "Theta", 1:K] <- sapply(1:K, \(k) PrecRecAUC(A = tmp2max3$Tht[1:p, 1:p, k, 1L, ], B = Theta.y[[k]])$auc)
    
    MSE_JCGLASSO[a, i, "Omega", 1:K] <- sapply(1:K, \(k) norm(tmp2max3$Omega[, , k, 1L, which.min(rowSums(tmpMSE))] - Theta.x[[k]], "F"))
    
    rho.opt <- tmp2max3$rho[which.min(rowSums(tmpMSE))]
    # lmax <- lmbMax_fun(tmp2max3, lambda.id = 1L, rho.id = which.min(rowSums(tmpMSE))) * 2
    # lmb.seq <- exp(seq(log(lmax), log(lmax * .01), l = 10))
    
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
  
  # if(i %% 5 == 0) save.image(paste0("~/Downloads/Simul2.2_14032024_n", n, "p", p, "q", q,
  #                                   "NA", perc.NA*100, "Xna", perc.X.na*100, "Yna", perc.Y.na*100,
  #                                   "alpha", alpha.true, "_comparison.RData"))
}

# MSE_JMMLE[,,1] <- MSE_JMMLE[,,1] * rep(c((p*(q+1)), sum(U), sum(U.2)), each = n.sim)
# MSE_JMMLE[,,2] <- MSE_JMMLE[,,2] * rep(c((p*(q+1)), sum(U), sum(U.2)), each = n.sim)
# 
# MSE_CGLASSO[,,1] <- MSE_CGLASSO[,,1] * rep(c((p*(q+1)), sum(U), sum(U.2)), each = n.sim)
# MSE_CGLASSO[,,2] <- MSE_CGLASSO[,,2] * rep(c((p*(q+1)), sum(U), sum(U.2)), each = n.sim)
# 
# MSE_JCGLASSO[1,,,1] <- MSE_JCGLASSO[1,,,1] * rep(c((p*(q+1)), sum(U), sum(U.2)), each = n.sim)
# MSE_JCGLASSO[1,,,2] <- MSE_JCGLASSO[1,,,2] * rep(c((p*(q+1)), sum(U), sum(U.2)), each = n.sim)
# 
# MSE_JCGLASSO[2,,,1] <- MSE_JCGLASSO[2,,,1] * rep(c((p*(q+1)), sum(U), sum(U.2)), each = n.sim)
# MSE_JCGLASSO[2,,,2] <- MSE_JCGLASSO[2,,,2] * rep(c((p*(q+1)), sum(U), sum(U.2)), each = n.sim)
# 
# MSE_JCGLASSO[3,,,1] <- MSE_JCGLASSO[3,,,1] * rep(c((p*(q+1)), sum(U), sum(U.2)), each = n.sim)
# MSE_JCGLASSO[3,,,2] <- MSE_JCGLASSO[3,,,2] * rep(c((p*(q+1)), sum(U), sum(U.2)), each = n.sim)


MSE1 <- sapply(1:K, \(k) apply(MSE_JMMLE[, , k], 2, mean))
MSE2 <- sapply(1:K, \(k) apply(MSE_CGLASSO[, , k], 2, mean))

MSE3 <- lapply(1:length(alpha.seq), \(a) 
               sapply(1:K, \(k) apply(MSE_JCGLASSO[a, , , k], 2, mean)))

AUC1 <- sapply(1:K, \(k) apply(AUC_JMMLE[, , k], 2, mean))
AUC2 <- sapply(1:K, \(k) apply(AUC_CGLASSO[, , k], 2, mean))

AUC3 <- lapply(1:length(alpha.seq), \(a) 
               sapply(1:K, \(k) apply(AUC_JCGLASSO[a, , , k], 2, mean)))

TAB <- rbind(
  cbind(AUC1[3,], AUC2[3,], AUC3[[1]][3,], AUC3[[2]][3,], AUC3[[3]][3,]),
  cbind(MSE1[3,], MSE2[3,], MSE3[[1]][3,], MSE3[[2]][3,], MSE3[[3]][3,]),
  cbind(AUC1[1,], AUC2[1,], AUC3[[1]][1,], AUC3[[2]][1,], AUC3[[3]][1,]),
  cbind(MSE1[1,], MSE2[1,], MSE3[[1]][1,], MSE3[[2]][1,], MSE3[[3]][1,]),
  cbind(AUC1[2,], AUC2[2,], AUC3[[1]][2,], AUC3[[2]][2,], AUC3[[3]][2,]),
  cbind(MSE1[2,], MSE2[2,], MSE3[[1]][2,], MSE3[[2]][2,], MSE3[[3]][2,])
)

# clip <- pipe("pbcopy", "w")
# write.table(round(TAB, 3L), file = clip, sep = "\t", row.names = TRUE)
# close(clip)
# 
# clip <- pipe("pbcopy", "w")
# write.table(round(cbind(t(AUC1), t(AUC2), t(AUC3[[1]]), t(AUC3[[2]]), t(AUC3[[3]])), 3L), file = clip, sep = "\t", row.names = TRUE)
# close(clip)

