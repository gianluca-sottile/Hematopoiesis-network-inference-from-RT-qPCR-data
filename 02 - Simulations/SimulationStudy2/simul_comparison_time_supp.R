rm(list = ls())

# libraries
library(MASS)
library(cglasso)
library(JGL)
library(huge)
library(ggplot2)
library(caret)
library(BDgraph)
library(microbenchmark)

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

n <- c(100, 200, 500)           # sample size
p <- c(25, 50, 100, 200)        # number of response variables
q <- c(25, 50)                  # number of predictors

grid <- expand.grid("n" = n, "p" = p, "q" = q)

ncomp <- 5L
n.sim <- 5L

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

# # Output objects
time_comparison <- matrix(0.0, nrow = nrow(grid), ncol = 5L,
                          dimnames = list(1:nrow(grid), c("JMMLE", "CGLASSO", "JCGL025", "JCGL050", "JCGL075")))

#############################################
# starting simulation study
#############################################

for(j in 1:nrow(grid)) {
  perc <- rep(1E-6, grid[j, "p"])
  
  set.seed(1234)
  
  X <- Ximp <- Theta.x <- Sigma.x <- B <- mu <- E <- 
    Y <- Yimp <- Theta.y <- Sigma.y <- Z <- Zimp <- ZX <- vector(mode = "list", length = K)
  unconnectedTheta <- list(NULL, alpha.true)
  unconnectedOmg <- list(NULL, alpha.true)
  
  indexB <- sapply(seq_len(grid[j, "p"]), \(i) sort(sample(grid[j, "q"], size = 2L, replace = FALSE)))
  unconnectedB <- list(NULL, alpha.true)
  
  for(k in 1:K) {
    Theta.y[[k]] <- simulTheta(p = grid[j, "p"], unconnected = unconnectedTheta[[k]],
                               ncomp = ncomp, mintht = Thtmin, maxtht = Thtmax)
    Sigma.y[[k]] <- solve(Theta.y[[k]])
    
    Theta.x[[k]] <- simulOmega(p = grid[j, "q"], unconnected = unconnectedTheta[[k]],
                               ncomp = ncomp, mintht = Omgmin, maxtht = Omgmax)
    Sigma.x[[k]] <- solve(Theta.x[[k]])
    
    X[[k]] <- mvrnorm(n = grid[j, "n"], mu = rep(0.0, grid[j, "q"]), Sigma = Sigma.x[[k]])
    B[[k]] <- simulB(X[[k]], perc, up, grid[j, "q"], grid[j, "p"], Bmin, Bmax, indexB, unconnectedB[[k]])
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
  
  for(k in 1:K) {
    X[[k]] <- mvrnorm(n = grid[j, "n"], mu = rep(0.0, grid[j, "q"]), Sigma = Sigma.x[[k]])
    B[[k]] <- simulB(X[[k]], perc, up, grid[j, "q"], grid[j, "p"], Bmin, Bmax, indexB, unconnectedB[[k]])
    eta <- cbind(1.0, X[[k]]) %*% B[[k]]
    mu[[k]] <- eta
    
    X[[k]][sample(grid[j, "n"], floor(grid[j, "n"] * perc.NA)), sample(grid[j, "q"], floor(grid[j, "q"] * perc.X.na))] <- NA
    Ximp[[k]] <- impute(X[[k]])
    
    E[[k]] <- mvrnorm(n = grid[j, "n"], mu = rep(0, grid[j, "p"]), Sigma = Sigma.y[[k]])
    Y[[k]] <- mu[[k]] + E[[k]]
    
    Y[[k]][sample(grid[j, "n"], floor(grid[j, "n"] * perc.NA)), sample(grid[j, "p"], floor(grid[j, "p"] * perc.Y.na))] <- NA
    Yimp[[k]] <- impute(Y[[k]])
  }
  ZX <- datajcggm(Y = X) 
  Z <- datajcggm(Y = Y, X = X) 
  Zimp <- datajcggm(Y = Yimp, X = Ximp)
  
  ##### JMMLE SECTION #####
  timeJMMLE <- microbenchmark(tmpJMMLE <- JMMLEcust(Ximp, Theta.x, Yimp, Theta.y, B, 
                                                 lambda.vec = sqrt(log(grid[j, "q"])/grid[j, "n"]) * exp(seq(log(1), log(.1), l = 5)), 
                                                 gamma.vec = sqrt(log(grid[j, "p"])/grid[j, "n"]) * exp(seq(log(1), log(.1), l = 5))),
                              times = n.sim, unit = "s")
  
  ##### CGLASSO SECTION #####
  tmpCGLASSO <- vector("list", length = K)
  timeCGLASSO <- microbenchmark(
    for(k in 1:K){
      tmpZ <- datacggm(Y = Y[[k]], X = Ximp[[k]])
      tmpCGLASSO[[k]] <- cglasso(data = tmpZ, nlambda = 5L, lambda.min.ratio = .01,
                                 nrho = 5L, rho.min.ratio = .01, trace = 0L, thr.bcd = 1E-5, thr.em = 1E-3)
    }, times = n.sim, unit = "s")
  
  ##### JCGLASSO SECTION #####
  # Simulation varying alpha
  timeJCGLASSO_025 <- microbenchmark(tmpJCGLASSO_025 <- jcglasso(data = Z, nlambda = 5L, lambda.min.ratio = .01, nu = select(tmpX, BIC)$rho,
                                                                 nrho = 5L, rho.min.ratio = .01, trace = 0L, thr.bcd = 1E-4, thr.em = 1E-3,
                                                                 alpha1 = alpha.seq[1], alpha2 = alpha.seq[1], alpha3 = alpha.seq[1]), 
                                     times = n.sim, unit = "s")
  
  timeJCGLASSO_050 <- microbenchmark(tmpJCGLASSO_025 <- jcglasso(data = Z, nlambda = 5L, lambda.min.ratio = .01, 
                                                                 nrho = 5L, rho.min.ratio = .01, trace = 0L, thr.bcd = 1E-4, thr.em = 1E-3,
                                                                 alpha1 = alpha.seq[2], alpha2 = alpha.seq[2], alpha3 = alpha.seq[2]),
                                     times = n.sim, unit = "s")
  
  timeJCGLASSO_075 <- microbenchmark(tmpJCGLASSO_025 <- jcglasso(data = Z, nlambda = 5L, lambda.min.ratio = .01, 
                                                                 nrho = 5L, rho.min.ratio = .01, trace = 0L, thr.bcd = 1E-4, thr.em = 1E-3,
                                                                 alpha1 = alpha.seq[3], alpha2 = alpha.seq[3], alpha3 = alpha.seq[3]), 
                                     times = n.sim, unit = "s")
    
  time_comparison[j, ] <- c(summary(timeJMMLE, "s")[1, "median"],
                            summary(timeCGLASSO, "s")[1, "median"],
                            summary(timeJCGLASSO_025, "s")[1, "median"],
                            summary(timeJCGLASSO_050, "s")[1, "median"],
                            summary(timeJCGLASSO_075, "s")[1, "median"])  
  
  ##### PRINTING SECTION #####
  cat("n =", grid[j, "n"], "\tp =", grid[j, "p"], "\tq =", grid[j, "q"])
  cat("\nComputational time\n")
  print(time_comparison[j, ])
  cat("\n")
  
  # if(i %% 5 == 0) save.image("~/Downloads/Simul3_130324_computational_time.RData")
}

ALL <- rbind(cbind(grid, "value" = time_comparison[, 1])[1:24, ],
             cbind(grid, "value" = time_comparison[, 2])[1:24, ],
             cbind(grid, "value" = time_comparison[, 3])[1:24, ],
             cbind(grid, "value" = time_comparison[, 4])[1:24, ],
             cbind(grid, "value" = time_comparison[, 5])[1:24, ])
ALL$method <- factor(rep(c("JMMLE", "CGLASSO", "JCGL025", "JCGL050", "JCGL075"), each = 24),
                     levels = c("JMMLE", "CGLASSO", "JCGL025", "JCGL050", "JCGL075"))
ALL$q <- factor(ALL$q, labels = c("q = 25", "q = 50"))
ALL$p <- factor(ALL$p, labels = c("p = 25", "p = 50", "p = 100", "p = 200"))


library(ggplot2)

p1 <- ggplot(ALL) + geom_line(aes(x = n, y = value, col = method), size = 1.2) + 
  facet_grid(rows = vars(q), cols = vars(p)) + theme_bw() + 
  labs(y = "Computational time (seconds)", x = "Sample size (n)") + 
  theme(legend.title = element_blank(), legend.position = "bottom", legend.text = element_text(size = 14),
        strip.text = element_text(size = 14), axis.title = element_text(size = 14), axis.text = element_text(size = 12)
        )
p1
# ggsave("~/Downloads/computational_time.pdf", device = "pdf", scale = .75, width = 14, height = 10, units = "in")
