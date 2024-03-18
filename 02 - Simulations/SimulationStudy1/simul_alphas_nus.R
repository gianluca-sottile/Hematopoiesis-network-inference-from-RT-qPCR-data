.nu.vec <- c("0", "BIC", "MAX")
rm(list = ls())

# libraries
library(MASS)
library(cglasso)
library(JGL)
library(huge)
library(ggplot2)
library(caret)

source("../../01 - RCode/datajcggm.R")
source("../../01 - RCode/jcglasso.R")
source("../../01 - RCode/jcggm.R")
source("../../01 - RCode/gof.R")
source("../../01 - RCode/to_graph.R")

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
  tmp <- tmp[, -c(1, ncol(tmp))]
  list(precision = tmp["Precision", ], recall = tmp["Recall", ], auc = auc)
}

n <- 100L                # sample size
p <- 200L               # number of response variables
q <- 50L                 # number of predictors
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

# Output objects
Recall <- Precision <- array(0.0, dim = c(length(.nu.vec), length(alpha.seq), n.sim, 2L, 
                                          length(c2), length(c3) + 2L, K + 1L), 
                             dimnames = list(nu = .nu.vec, 
                                             alpha = alpha.seq, 
                                             sim = seq_len(n.sim), 
                                             statistic = c("B", "Theta"),
                                             proportionMax = c2,
                                             propCurve = c("MAX", c3, "MIN"), 
                                             class.id = c(seq_len(K), "mean")))

MSE <- array(0.0, dim = c(length(.nu.vec), length(alpha.seq), n.sim, 2L, 
                          length(c2), length(c3), K + 1L), 
             dimnames = list(nu = .nu.vec, 
                             alpha = alpha.seq, 
                             sim = seq_len(n.sim), 
                             statistic = c("B", "Theta"),
                             proportionMax = c2,
                             propCurve = c3, 
                             class.id = c(seq_len(K), "mean")))

AUC <- array(0.0, dim = c(length(.nu.vec), length(alpha.seq), n.sim, 2L, length(c2), K + 1L), 
             dimnames = list(nu = .nu.vec,
                             alpha = alpha.seq, 
                             sim = seq_len(n.sim), 
                             statistic = c("B", "Theta"),
                             proportionMax = c2,
                             class.id = c(seq_len(K), "mean")))

for(.nu in .nu.vec) {
  cat("\nSetting for nu =", .nu, "\n")
  
  #############################################
  # starting simulation study
  #############################################
  set.seed(1234)
  
  X <- Theta.x <- Sigma.x <- B <- mu <- E <- 
    Y <- Theta.y <- Sigma.y <- Z <- ZX <- vector(mode = "list", length = K)
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
  
  # Main simulation
  for(i in 1:n.sim) {
    for(k in 1:K) {
      X[[k]] <- mvrnorm(n = n, mu = rep(0.0, q), Sigma = Sigma.x[[k]])
      B[[k]] <- simulB(X[[k]], perc, up, q, p, Bmin, Bmax, indexB, unconnectedB[[k]])
      eta <- cbind(1.0, X[[k]]) %*% B[[k]]
      mu[[k]] <- eta
      
      X[[k]][sample(n, floor(n * perc.NA)), sample(q, floor(q * perc.X.na))] <- NA
      
      E[[k]] <- mvrnorm(n = n, mu = rep(0, p), Sigma = Sigma.y[[k]])
      Y[[k]] <- mu[[k]] + E[[k]]
      
      Y[[k]][sample(n, floor(n * perc.NA)), sample(p, floor(p * perc.Y.na))] <- NA
    }
    ZX <- datajcggm(Y = X) 
    Z <- datajcggm(Y = Y, X = X) 
    
    # Simulation varying alpha
    for(a in 1:length(alpha.seq)) {
      
      # Simulation for B fixing a grid of rho: c(1, .75, .5, .25, .1)*rMAx
      if(.nu == "0") nu <- 0.0
      if(.nu == "BIC") {
        tmpX <- jcglasso(data = ZX, nrho = 25L, rho.min.ratio = 0.01,
                         trace = 0L, thr.bcd = 1E-5, thr.em = 1E-3, alpha1 = alpha.seq[a])
        plot(BIC(tmpX))
        plot(select.jcglasso(tmpX, BIC))
        nu <- select.jcglasso(tmpX, BIC)$rho
      }
      if(.nu == "MAX") {
        tmpX <- jcglasso(data = ZX, nrho = 1L, trace = 0L, thr.bcd = 1E-5, thr.em = 1E-3, alpha1 = alpha.seq[a])
        nu <- tmpX$rho
      }
      
      tmp2 <- jcglasso(data = Z, nlambda = 1L, nrho = 4L, rho.min.ratio = 0.25, nu = nu,
                       trace = 0L, thr.bcd = 1E-4, thr.em = 1E-3, 
                       alpha1 = alpha.seq[a], alpha2 = alpha.seq[a], alpha3 = alpha.seq[a])
      for(r in 1:length(c2)) {
        lMax2 <- lmbMax_fun(tmp2, lambda.id = 1L, rho.id = r)
        # lMax2 <- max(abs(t(X[[1]]) %*% tmp2$R[, , 1L, r] %*% tmp2$Tht[, , 1L, r]) / matrix(diag(tmp2$Tht[, , 1L, r]), nrow = q, ncol = p, byrow = TRUE)) / n
        lambda.seq <- lMax2 * c3
        tmp2.2 <- jcglasso(data = Z, lambda = lambda.seq, rho = tmp2$rho[r], nu = nu,
                           trace = 0L, thr.bcd = 1E-4, thr.em = 1E-3, 
                           alpha1 = alpha.seq[a], alpha2 = alpha.seq[a], alpha3 = alpha.seq[a])
        
        #### MSE of B (marginalising on rho)
        MSE[match(.nu, .nu.vec), a, i, "B", r, , seq_len(K)] <- sapply(seq_len(K), \(k) apply(tmp2.2$B[, , k, , 1L], 3, function(Bh) norm(Bh - B[[k]], "F")))
        MSE[match(.nu, .nu.vec), a, i, "B", r, , K + 1L] <- if(K == 1L) MSE[match(.nu, .nu.vec), a, i, "B", r, , seq_len(K)] else rowMeans(MSE[match(.nu, .nu.vec), a, i, "B", r, , seq_len(K)])
        
        #### PR curve and AUC of B (marginalising on rho)
        for(k in seq_len(K)) {
          out_prcurve <- PrecRecAUC(A = tmp2.2$B[, , k, , 1L], B = B[[k]])
          Precision[match(.nu, .nu.vec), a, i, "B", r, , k] <- out_prcurve$precision
          Recall[match(.nu, .nu.vec), a, i, "B", r, , k] <- out_prcurve$recall
          AUC[match(.nu, .nu.vec), a, i, "B", r, k] <- out_prcurve$auc
        }
        
        Precision[match(.nu, .nu.vec), a, i, "B", r, , K + 1L] <- if(K == 1L) Precision[match(.nu, .nu.vec), a, i, "B", r, , seq_len(K)] else rowMeans(Precision[match(.nu, .nu.vec), a, i, "B", r, , seq_len(K)])
        Recall[match(.nu, .nu.vec), a, i, "B", r, , K + 1L] <- if(K == 1L) Recall[match(.nu, .nu.vec), a, i, "B", r, , seq_len(K)] else rowMeans(Recall[match(.nu, .nu.vec), a, i, "B", r, , seq_len(K)])
        AUC[match(.nu, .nu.vec), a, i, "B", r, K + 1L] <- if(K == 1L) AUC[match(.nu, .nu.vec), a, i, "B", r, seq_len(K)] else mean(AUC[match(.nu, .nu.vec), a, i, "B", r, seq_len(K)])
      }
      
      # Simulation for Theta fixing a grid of lambda: c(1, .75, .5, .25, .1)*lMAx
      tmp3 <- jcglasso(data = Z, nlambda = 4L, nrho = 1L, lambda.min.ratio = 0.25, nu = nu,
                       trace = 0L, thr.bcd = 1E-4, thr.em = 1E-3, 
                       alpha1 = alpha.seq[a], alpha2 = alpha.seq[a], alpha3 = alpha.seq[a])
      
      for(l in 1:length(c2)) {
        rMax2 <- rhoMax_fun(tmp3, lambda.id = l, rho.id = 1L)
        # rMax2 <- max(abs(tmp3$S[, , l, 1L][outer(1:p, 1:p, "<")]))
        rho.seq <- rMax2 * c3
        tmp3.2 <- jcglasso(data = Z, lambda = tmp3$lambda[l], rho = rho.seq, nu = nu,
                           trace = 0L, thr.bcd = 1E-4, thr.em = 1E-3, 
                           alpha1 = alpha.seq[a], alpha2 = alpha.seq[a], alpha3 = alpha.seq[a])
        
        #### MSE of Theta (marginalising on lambda)
        MSE[match(.nu, .nu.vec), a, i, "Theta", l, , seq_len(K)] <- sapply(seq_len(K), \(k) apply(tmp3.2$Tht[1:p, 1:p, k, 1L, ], 3L, function(Th) norm(Th - Theta.y[[k]], "F")))
        MSE[match(.nu, .nu.vec), a, i, "Theta", l, , K + 1L] <- if(K == 1L) MSE[match(.nu, .nu.vec), a, i, "Theta", l, , seq_len(K)] else rowMeans(MSE[match(.nu, .nu.vec), a, i, "Theta", l, , seq_len(K)])
        
        #### PR curve and AUC of B (marginalising on rho)
        for(k in seq_len(K)) {
          out_prcurve <- PrecRecAUC(A = tmp3.2$Tht[1:p, 1:p, k, 1L, ], B = Theta.y[[k]])
          Precision[match(.nu, .nu.vec), a, i, "Theta", l, , k] <- out_prcurve$precision
          Recall[match(.nu, .nu.vec), a, i, "Theta", l, , k] <- out_prcurve$recall
          AUC[match(.nu, .nu.vec), a, i, "Theta", l, k] <- out_prcurve$auc
        }
        
        Precision[match(.nu, .nu.vec), a, i, "Theta", l, , K + 1L] <- if(K == 1L) Precision[match(.nu, .nu.vec), a, i, "Theta", l, , seq_len(K)] else rowMeans(Precision[match(.nu, .nu.vec), a, i, "Theta", l, , seq_len(K)])
        Recall[match(.nu, .nu.vec), a, i, "Theta", l, , K + 1L] <- if(K == 1L) Recall[match(.nu, .nu.vec), a, i, "Theta", l, , seq_len(K)] else rowMeans(Recall[match(.nu, .nu.vec), a, i, "Theta", l, , seq_len(K)])
        AUC[match(.nu, .nu.vec), a, i, "Theta", l, K + 1L] <- if(K == 1L) AUC[match(.nu, .nu.vec), a, i, "Theta", l, seq_len(K)] else mean(AUC[match(.nu, .nu.vec), a, i, "Theta", l, seq_len(K)])
      }
      
      cat("iter n. ", i, "and alpha =", alpha.seq[a])
      cat("\nMSE of B")
      print(MSE[match(.nu, .nu.vec), a, i, "B", , , K + 1L])
      cat("\nMSE of Theta")
      print(MSE[match(.nu, .nu.vec), a, i, "Theta", , , K + 1L])
      cat("\nAUC of B and Theta")
      print(AUC[match(.nu, .nu.vec), a, i, , , K + 1L])
      cat("\n")
    }
    
    # if(i %% 5 == 0) save.image(paste0("~/Downloads/Simul1.1_13032024_n", n, "p", p, "q", q,
    #                                   "NA", perc.NA*100, "Xna", perc.X.na*100, "Yna", perc.Y.na*100,
    #                                   "alpha", alpha.true, "nu", .nu, ".RData"))
  }
}
  
AUC.all <- MSE.all <- NULL
MSE.new <- AUC.new <- vector(mode = "list", length = length(.nu.vec))
for(.nu in .nu.vec) {
  n <- 100L                # sample size
  p <- 200L                # number of response variables
  q <- 50L                 # number of predictors
  perc.NA <- 0.25
  perc.X.na <- 0.25
  perc.Y.na <- 0.25

  # load(paste0("~/Downloads/Simul1.1_13032024_n", n, "p", p, "q", q,
  #             "NA", perc.NA*100, "Xna", perc.X.na*100, "Yna", perc.Y.na*100,
  #             "alpha", alpha.true, "nu", .nu, ".RData"))
  
  AUC.mean <- lapply(c("B", "Theta"), \(stat)
                     lapply(1:length(alpha.seq), \(a)
                            colMeans(AUC[match(.nu, .nu.vec), a, , stat, , K + 1L])))
  names(AUC.mean) <- c("B", "Theta")
  AUC.mean <- lapply(AUC.mean, \(x) {
    names(x) <- alpha.seq
    x
  })
  
  AUC.new[[match(.nu, .nu.vec)]] <- AUC.mean
  
  MSE.mean <- lapply(c("B", "Theta"), \(stat)
                     lapply(1:length(alpha.seq), \(a)
                            sapply(1:length(c3), \(j)
                                   colSums(apply(MSE[match(.nu, .nu.vec), a, , stat , , j, 1:K], 2, colMeans)))))# colMeans(MSE[match(.nu, .nu.vec), a, , stat , , j, K + 1L]))))
                                   
  names(MSE.mean) <- c("B", "Theta")
  MSE.mean <- lapply(MSE.mean, \(x) {
    x <- lapply(x, \(y) {
      colnames(y) <- c3
      y
    })
    names(x) <- alpha.seq
    x
  })
  
  MSE.new[[match(.nu, .nu.vec)]] <- MSE.mean
  
  AUC.df <- data.frame(value = c(unlist(AUC.mean[["B"]]), unlist(AUC.mean[["Theta"]])), 
                       prop = rep(c2, length(alpha.seq) * 2), 
                       alpha = factor(rep(rep(alpha.seq, each = length(c2)), 2), levels = c(alpha.seq)),
                       measure = factor(rep(c("B", "Theta"), each = length(alpha.seq) * length(c2)), 
                                        levels = c("B", "Theta"), 
                                        labels = c("B (rho | lambda_seq)", "Theta (lambda | rho_seq)")),
                       setting = .nu)
  AUC.all <- rbind(AUC.all, AUC.df)
  
  MSE.df <- data.frame(value = c(unlist(MSE.mean[["B"]]), unlist(MSE.mean[["Theta"]])),
                       prop = rep(c2, length(c3) * length(alpha.seq) * 2),
                       prop2 = rep(rep(c3, each = length(c2)), length(alpha.seq) * 2),
                       alpha = factor(rep(rep(alpha.seq, each = length(c2) * length(c3)), 2), levels = c(alpha.seq)),
                       measure = rep(c("B", "Theta"), each = length(alpha.seq) * length(c3) * length(c2)), setting = .nu)
  
  MSE.all <- rbind(MSE.all, MSE.df)
}

names(MSE.new) <- names(AUC.new) <- .nu.vec

MSE.all$prop <- factor(MSE.all$prop, levels = rev(c2), 
                       labels = c("0.25", "0.50", "0.75", "1.00"))

MSE.all$alpha <- factor(MSE.all$alpha, levels = alpha.seq, 
                        labels = c("0.25", "0.50", "0.75"))

AUC.all$alpha <- factor(AUC.all$alpha, levels = alpha.seq, 
                        labels = c("0.25", "0.50", "0.75"))

MSE.all$setting <- factor(MSE.all$setting,
                          levels = .nu.vec,
                          labels = c(
                            expression(paste(nu, " = 0.00")),
                            expression(paste(nu, " = ", nu[BIC])),
                            expression(paste(nu, " = ", nu[max]))))

AUC.all$setting <- factor(AUC.all$setting,
                          levels = .nu.vec,
                          labels = c(
                            expression(paste(nu, " = 0.00")),
                            expression(paste(nu, " = ", nu[BIC])),
                            expression(paste(nu, " = ", nu[max]))))

MSE.allB <- subset(MSE.all, measure == "B")
MSE.allB$prop <- factor(MSE.allB$prop, 
                            levels = levels(MSE.allB$prop),
                            labels = c(
                              expression(paste(rho / rho[max],  " = 0.25")),
                              expression(paste(rho / rho[max],  " = 0.50")),
                              expression(paste(rho / rho[max],  " = 0.75")),
                              expression(paste(rho / rho[max],  " = 1.00"))))

p1 <- ggplot(data = MSE.allB, aes(x = prop2, y = value, color = alpha)) + 
  geom_line(size = .75) + 
  geom_point(size = 2.5) + theme_bw() + 
  scale_color_manual(values = c("red", "blue", "green")) +
  facet_grid(setting ~ prop, labeller = labeller(setting = label_parsed, prop = label_parsed)) + 
  xlab(expression(lambda/lambda[max])) + ylab(bquote(Mean~Squared~Error~of~~'{'~hat(B)~'}')) + 
  theme(axis.text = element_text(size = 14), axis.title = element_text(size = 14), 
        legend.text = element_text(size = 14), legend.title = element_text(size = 14),
        strip.text = element_text(size = 14)) + 
  guides(color = guide_legend(title = "Alpha"))
p1

MSE.allTheta <- subset(MSE.all, measure == "Theta")
MSE.allTheta$prop <- factor(MSE.allTheta$prop, 
                            levels = levels(MSE.allTheta$prop),
                            labels = c(
                              expression(paste(lambda / lambda[max],  " = 0.25")),
                              expression(paste(lambda / lambda[max],  " = 0.50")),
                              expression(paste(lambda / lambda[max],  " = 0.75")),
                              expression(paste(lambda / lambda[max],  " = 1.00"))))

p2 <- ggplot(data = MSE.allTheta, aes(x = prop2, y = value, color = alpha)) + 
  geom_line(size = .75) + 
  geom_point(size = 2.5) + theme_bw() + 
  scale_color_manual(values = c("red", "blue", "green")) +
  facet_grid(setting ~ prop, labeller = labeller(setting = label_parsed, prop = label_parsed)) + 
  xlab(expression(rho/rho[max])) + ylab(bquote(Mean~Squared~Error~of~~'{'~hat(Theta)~'}')) + 
  theme(axis.text = element_text(size = 14), axis.title = element_text(size = 14), 
        legend.text = element_text(size = 14), legend.title = element_text(size = 14),
        strip.text = element_text(size = 14)) + 
  guides(color = guide_legend(title = "Alpha"))
p2

library(gridExtra)
library(lemon)
pall <- grid_arrange_shared_legend(p1, p2, nrow = 2, ncol = 1, position = "right")

# ggsave("~/Downloads/mse_curves_271023.pdf", plot = pall, device = "pdf", units = "in",
#        width = 12, height = 16, scale = .75)
 
AUC.all$measure <- factor(AUC.all$measure, 
                          levels = levels(AUC.all$measure),
                          labels = c(
                            expression(paste("{", B,"}")),
                            expression(paste("{", Theta,"}"))))

p3 <- ggplot(data = AUC.all, aes(x = prop, y = value, color = alpha)) + 
  geom_line(size = .75) + 
  geom_point(size = 2.5) + theme_bw() + 
  scale_color_manual(values = c("red", "blue", "green")) +
  facet_grid(setting ~ measure, labeller = labeller(setting = label_parsed, measure = label_parsed)) + 
  # xlab(expression(lambda/lambda[max])) + 
  ylab("Area under the precision-recall curve") + 
  theme(axis.text = element_text(size = 14), axis.title = element_text(size = 14), 
        legend.text = element_text(size = 14), legend.title = element_text(size = 14),
        strip.text = element_text(size = 14)) + 
  guides(color = guide_legend(title = "Alpha"))
p3

# ggsave("~/Downloads/auc_curves_271023.pdf", plot = p3, device = "pdf", units = "in",
#        width = 12, height = 16, scale = .75)
