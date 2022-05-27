############################################################################################
#
# Author:   Gianluca Sottile and Luigi Augugliaro
# e-mail:   gianluca.sottile@unipa.it
# e-mail:   luigi.augugliaro@unipa.it
# data:     16-09-2021
#
# Description: Comparisons between Joint glasso and Joint cglasso
#
# Journal: RSS-Series C
############################################################################################

index <- expand.grid(i.H = c(1,2), i.p = c(1,2), i.perc = c(1,2))

# for(.id in c(5, 6, 7, 8)) {
  .id <- 8

  # clean the workspace
  rm(list = setdiff(ls(), "index"))
  gc(reset = TRUE)
  
  i.H <- index[.id, 1]
  i.p <- index[.id, 2]
  i.q <- 1L
  i.perc <- index[.id, 3]
  
  model <- "censoring"
    
  alpha <- .25
  # alpha <- .5
  # alpha <- .75
  
  percorso <- paste0("results/01_alpha_", 
                     ifelse(alpha == .25, "025", ifelse(alpha == .5, "050", "075")))
  if(!dir.exists(percorso)) dir.create(percorso)
  
  # libraries
  library(MASS)
  library(cglasso)
  library(JGL)
  library(huge)
  source("../../../../01 - RCode/jcglasso.R")
  
  simulTheta <- function(p, ncomp, unconnected = NULL, mintht = 0.4, maxtht = 0.5){
    Tht <- matrix(0, p, p)
    subTht <- diag(ncomp)
    subTht[1L, ] <- subTht[, 1L] <- sort(runif(ncomp, min = mintht, maxtht), TRUE)
    block <- matrix(seq_len(p), nrow = ncomp)
    if(!is.null(unconnected)) block <- block[, -unconnected]
    for(h in seq_len(ncol(block))) Tht[block[, h], block[, h]] <- subTht
    diag(Tht) <- 1
    Tht
  }
  
  nruns <- 10
  subnsim <- 5
  nsim <- nruns * subnsim
  up <- 40                        # right censoring value
  n <- 100                        # sample size
  p <- c(50, 200)                 # number of response variables
  q <- c(50, 200)                 # number of predictors
  K <- 3                          # number of groups
  pH <- c(0.20, 0.40)
  perc0 <- c(0.20, 0.40)
  unconnected <- list(NULL, 9, 9:10)
  c2 <- c(1, .75, .5, .25, .1) #seq(1, .1, l = 10)
  Bmin <- 0.3
  Bmax <- 0.7
  
  # rho / rho_max = c(1, .75, .5, .25, .1, .01) => rho = rho_max * c(1, .75, .5, .25, .1, .01)
  
  H <- p[i.p] * pH[i.H]
  H2 <- q[i.q] * pH[i.H]
  nrho1 <- length(c2)                              # number of rho1-values
  rho1.min.ratio <- 1
  nlambda1 <- length(c2)                           # number of lambda1-values
  lambda1.min.ratio <- 1
  
  perc <- rep(c(perc0[i.perc], 1e-6), c(H, p[i.p] - H))         # percentage of expected missing for Y
  perc2 <- rep(c(perc0[i.perc], 1e-6), c(H2, q[i.q] - H2))      # percentage of expected missing for X
  
  ##############
  # output
  ##############
  conv <- matrix("", nsim, 1)
  colnames(conv) <- c("jcglasso")
  
  # filenames <- paste0(paste0(percorso, "/run."), seq_len(nruns), ".RData")
  filenames <- paste0(percorso, "/run.", seq_len(nruns), ".RData")
  
  Tht.dim <- c(p[i.p], p[i.p], K, nlambda1, nrho1)
  Thth <- array(0, dim = Tht.dim)
  Tht.array <- array(0, dim = c(prod(Tht.dim), 1, subnsim), 
                     dimnames = list(NULL, c("jcglasso"), NULL) )
  
  Thtx.dim <- c(q[i.q], q[i.q], K, nlambda1, nrho1)
  Ththx <- array(0, dim = Thtx.dim)
  Thtx.array <- array(0, dim = c(prod(Thtx.dim), 1, subnsim), 
                      dimnames = list(NULL, c("jcglasso"), NULL) )
  
  B.dim <- c(q[i.q] + 1L, p[i.p], K, nlambda1, nrho1)
  Bh <- array(0, dim = B.dim)
  B.array <- array(0, dim = c(prod(B.dim), 1, subnsim), 
                   dimnames = list(NULL, c("jcglasso"), NULL))
  
  #### MSE ####
  mse_B <- array(0, dim = c(nsim, K, nlambda1, nrho1, 1), 
                 dimnames = list(nsim = NULL, K = paste0("group", seq_len(K)),
                                 nlambda1 = 1L:nlambda1, nrho1 = 1:nrho1, 
                                 estimator = c("jcglasso")))
  
  mse_tbl_B <- array(0, dim = c(2, 1, nrho1, K), 
                     dimnames = list(type = c("mean", "sd"), 
                                     estimator = c("jcglasso"), 
                                     rho1 = c2,
                                     K = paste0("group", seq_len(K))))
  
  
  mse_Tht <- array(0, dim = c(nsim, K, nlambda1, nrho1, 1), 
                   dimnames = list(nsim = NULL, K = paste0("group", seq_len(K)), 
                                   nlambda1 = 1:nlambda1, nrho1 = 1L:nrho1, 
                                   estimator = c("jcglasso")))
  
  mse_tbl_Tht <- array(0, dim = c(2, 1, nlambda1, K), 
                       dimnames = list(type = c("mean", "sd"), 
                                       estimator = c("jcglasso"), 
                                       lambda1 = c2,
                                       K = paste0("group", seq_len(K))))
  
  
  mse_Thtx <- array(0, dim = c(nsim, K, nlambda1, nrho1, 1), 
                    dimnames = list(nsim = NULL, K = paste0("group", seq_len(K)), 
                                    nlambda1 = 1:nlambda1, nrho1 = 1L:nrho1, 
                                    estimator = c("jcglasso")))
  
  mse_tbl_Thtx <- array(0, dim = c(2, 1, nlambda1, K), 
                        dimnames = list(type = c("mean", "sd"), 
                                        estimator = c("jcglasso"), 
                                        lambda1 = c2,
                                        K = paste0("group", seq_len(K))))
  
  #### PR-CURVE and AUC ####
  
  prcurveTheta <- function(thetah, thetat, nrho){
    U <- upper.tri(thetat, diag = FALSE)
    thetat_v <- thetat[U]
    A <- which(abs(thetat_v) > 0)
    thetah_m <- apply(thetah, 3, function(M) M[U])
    
    recall <- if(length(A) == 1) apply(t(abs(thetah_m[A, ]) > 0), 2, mean) else apply(abs(thetah_m[A, ]) > 0, 2, mean)
    precision <- vector(mode = "numeric", length = nrho)
    
    for(i in 1:nrho){
      id <- which(abs(thetah_m[, i]) > 0)
      precision[i] <- ifelse(length(id) == 0, NA, mean(abs(thetat_v[id]) > 0))
    }
    
    # prcurve <- function(tab) c(precision = tab[1,1] / sum(tab[, 1]), recall = tab[1,1] / sum(tab[1,]))
    # tab <- t(sapply(1:nrho, function(i) prcurve(table(factor(1*(abs(thetat_v) > 0), levels = c(1,0)), 
    #                                          factor(1*(abs(thetah_m[, i]) > 0), levels = c(1,0))))))
    # # tab <- na.omit(tab)
    # min_prec <- length(A) / length(thetat_v)
    # tab <- rbind(c(1, 0), tab, c(min_prec, 1))
    
    # plot(tab[,2], tab[,1], type = "b", xlim = c(0,1), ylim = c(0,1))
    # lines(seq(0, 1, l=10), splinefun(tab[,2], tab[,1], method = "monoH.FC")(seq(0, 1, l=10)), col = 2)
    
    # auc <- integrate(splinefun(tab[,2], tab[,1], method = "monoH.FC"), 0, 1)$value
    
    id <- order(recall)
    precision <- precision[id]
    recall <- recall[id]
    
    min_prec <- length(A) / length(thetat_v)
    recall.auc <- c(0, recall, 1)
    precision.auc <- c(1, precision, min_prec)
    tab <- cbind(recall.auc, precision.auc)
    tab <- na.omit(tab)
    recall.auc <- tab[,1]
    precision.auc <- tab[,2]
    
    dprecision <- c(diff(precision.auc), 0)
    drecall <- c(diff(recall.auc), 0)
    auc <- sum(na.omit(precision.auc * drecall)) + sum(na.omit(dprecision * drecall)) / 2
    out <- list(recall = recall, precision = precision, auc = auc)
    # out <- list(recall = tab[-c(1,7),2], precision = tab[-c(1,7),1], auc = auc)
    out
  }
  
  prcurve_Tht <- array(0, dim = c(nsim, K, nlambda1, nrho1, 2, 1), 
                       dimnames = list(nsim = NULL, K = paste0("group", seq_len(K)), 
                                       lambda1 = 1L:nlambda1, nrho1 = 1L:nrho1, 
                                       type = c("recall", "precision"), 
                                       estimator = c("jcglasso")))
  
  auc_tbl_Tht <- array(0, dim = c(nsim, 1, K, nlambda1), 
                       dimnames = list(nsim = NULL, estimator = c("jcglasso"), 
                                       K = paste0("group", seq_len(K)), lambda1 = 1:nlambda1))
  
  prcurve_data_Tht <- array(0, dim = c(nlambda1, nrho1, K, 2, 1), 
                            dimnames = list(lambda1 = 1L:nlambda1, rho1 = 1L:nrho1, 
                                            K = paste0("group", seq_len(K)), 
                                            type = c("recall", "precision"), 
                                            estimator = c("jcglasso")))
  
  auc_Tht  <- array(0, dim = c(2, 1, nlambda1, K), 
                    dimnames = list(type = c("mean", "sd"), 
                                    estimator = c("jcglasso"), 
                                    lambda1 = NULL, K = paste0("group", seq_len(K))))
  
  prcurve_Thtx <- array(0, dim = c(nsim, K, nlambda1, nrho1, 2, 1), 
                        dimnames = list(nsim = NULL, K = paste0("group", seq_len(K)), 
                                        lambda1 = 1:nlambda1, nrho1 = 1L:nrho1, 
                                        type = c("recall", "precision"), 
                                        estimator = c("jcglasso")))
  
  auc_tbl_Thtx <- array(0, dim = c(nsim, 1, K, nlambda1), 
                        dimnames = list(nsim = NULL, estimator = c("jcglasso"), 
                                        K = paste0("group", seq_len(K)), lambda1 = 1:nlambda1))
  
  prcurve_data_Thtx <- array(0, dim = c(nlambda1, nrho1, K, 2, 1), 
                             dimnames = list(lambda1 = 1:nlambda1, rho1 = 1L:nrho1, 
                                             K = paste0("group", seq_len(K)), 
                                             type = c("recall", "precision"), 
                                             estimator = c("jcglasso")))
  
  auc_Thtx  <- array(0, dim = c(2, 1, nlambda1, K), 
                     dimnames = list(type = c("mean", "sd"), 
                                     estimator = c("jcglasso"), 
                                     lambda1 = NULL, K = paste0("group", seq_len(K))))
  
  prcurveB <- function(thetah, thetat, nrho){
    thetat_v <- c(thetat[-1, ])
    A <- which(abs(thetat_v) > 0)
    thetah_m <- apply(thetah, 3, function(M) c(M[-1,]))
    
    recall <- if(length(A) == 1) apply(t(abs(thetah_m[A, ]) > 0), 2, mean) else apply(abs(thetah_m[A, ]) > 0, 2, mean)
    precision <- vector(mode = "numeric", length = nrho)
    for(i in 1:nrho){
      id <- which(abs(thetah_m[, i]) > 0)
      precision[i] <- ifelse(length(id) == 0, NA, mean(abs(thetat_v[id]) > 0))
    }
    id <- order(recall)
    precision <- precision[id]
    recall <- recall[id]
    
    # prcurve <- function(tab) c(precision = tab[1,1] / sum(tab[, 1]), recall = tab[1,1] / sum(tab[1,]))
    # tab <- t(sapply(1:nrho, function(i) prcurve(table(factor(1*(abs(thetat_v) > 0), levels = c(1,0)), 
    #                                                   factor(1*(abs(thetah_m[, i]) > 0), levels = c(1,0))))))
    # 
    # min_prec <- length(A) / length(thetat_v)
    # tab <- rbind(c(1, 0), tab, c(min_prec, 1))
    # 
    # # plot(tab[,2], tab[,1], type = "b")
    # # lines(seq(0, 1, l=10), splinefun(tab[,2], tab[,1], method = "monoH.FC")(seq(0, 1, l=10)), col = 2)
    # 
    # auc <- integrate(splinefun(tab[,2], tab[,1], method = "monoH.FC"), 0, 1)$value
    # out <- list(recall = tab[-c(1,7),2], precision = tab[-c(1,7),1], auc = auc)
    
    min_prec <- length(A) / length(thetat_v)
    recall.auc <- c(0, recall, 1)
    precision.auc <- c(1, precision, min_prec)
    tab <- cbind(recall.auc, precision.auc)
    tab <- na.omit(tab)
    recall.auc <- tab[,1]
    precision.auc <- tab[,2]
    
    dprecision <- c(diff(precision.auc), 0)
    drecall <- c(diff(recall.auc), 0)
    auc <- sum(na.omit(precision.auc * drecall)) + sum(na.omit(dprecision * drecall)) / 2
    out <- list(recall = recall, precision = precision, auc = auc)
    
    out
  }
  
  prcurve_B <- array(0, dim = c(nsim, K, nlambda1, nrho1, 2, 1), 
                     dimnames = list(nsim = NULL, K = paste0("group", seq_len(K)), 
                                     lambda1 = 1L:nlambda1, nrho1 = 1:nrho1, 
                                     type = c("recall", "precision"), 
                                     estimator = c("jcglasso")))
  
  auc_tbl_B <- array(0, dim = c(nsim, 1, K, nrho1), 
                     dimnames = list(nsim = NULL, estimator = c("jcglasso"), 
                                     K = paste0("group", seq_len(K)), rho1 = 1:nrho1))
  
  prcurve_data_B <- array(0, dim = c(nlambda1, nrho1, K, 2, 1), 
                          dimnames = list(lambda1 = 1L:nlambda1, rho1 = 1:nrho1, 
                                          K = paste0("group", seq_len(K)), 
                                          type = c("recall", "precision"), 
                                          estimator = c("jcglasso")))
  
  auc_B  <- array(0, dim = c(2, 1, nrho1, K), 
                  dimnames = list(type = c("mean", "sd"), 
                                  estimator = c("jcglasso"), 
                                  rho1 = NULL, K = paste0("group", seq_len(K))))
  
  
  #############################################
  # starting simulation study
  #############################################
  set.seed(123)
  
  Xorig <- X <- Theta.x <- Sigma.x <- Theta.y <- Sigma.y <- B <- mu <- E <- Y <- S <- Z <- vector(mode = "list", length = K)
  
  for(k in seq_len(K)) {
    Theta.y[[k]] <- simulTheta(p = p[i.p], unconnected = unconnected[[k]], 
                               ncomp = 5, mintht = 0.3, maxtht = 0.5)
    Sigma.y[[k]] <- solve(Theta.y[[k]])
    
    sim <- huge.generator(n = n, d = q[i.q], graph = "hub", prob = 0.95)
    Sigma.x[[k]] <- round(sim$sigma, 5)
    Theta.x[[k]] <- round(solve(Sigma.x[[k]]), 5)
    
    X[[k]] <- mvrnorm(n = n, mu = rep(0, q[i.q]), Sigma = Sigma.x[[k]])
    
    B[[k]] <- matrix(0, nrow = q[i.q] + 1, ncol = p[i.p])
    B[[k]][-1L, ][1:2, ] <- runif(2*q[i.q], Bmin, Bmax)
    
    eta <- X[[k]] %*% B[[k]][-1L, ]
    
    cutoff <- function(b0, i, up, perc = perc) (mean(pnorm(up - b0 - eta[, i], lower.tail =  FALSE))) - perc
    B[[k]][1L, ] <- sapply(1:p[i.p], function(.i) uniroot(cutoff, interval = c(0, 100), i = .i, up = up, perc = perc[.i])$root)
    
    mu[[k]] <- cbind(1, X[[k]]) %*% B[[k]]
    
    Xorig[[k]] <- X[[k]]
    X[[k]] <- sapply(seq_len(q[i.q]), function(i) {
      x <- X[[k]][, i]
      id <- sample(n, floor(n * perc2[i]))
      x[id] <- NA
      x
    })
  }
  
  pb <- txtProgressBar(min = 0L, max = nsim, style = 3L)
  
  jj <- 0L
  
  U <- outer(1:p[i.p], 1:p[i.p], "<")
  U1 <- outer(1:p[i.p], 1:p[i.p], "<=")
  U2 <- outer(1:q[i.q], 1:q[i.q], "<=")
  
  for(h in 1:(nruns)) {
    for(j in 1:(subnsim)) {
      jj <- jj + 1L
      setTxtProgressBar(pb, jj)
      
      Yorig <- Y
      for(k in seq_len(K)) {
        E[[k]] <- mvrnorm(n = n, mu = rep(0, p[i.p]), Sigma = Sigma.y[[k]])
        Y[[k]] <- mu[[k]] + E[[k]]
        S[[k]] <- cov(Y[[k]]) * (n - 1) / n
        Yorig[[k]] <- Y[[k]]
        if(model == "censoring") {
          Y[[k]][Y[[k]] > up] <- up
          Z[[k]] <- datacjggm(Y[[k]], up = up, X = X[[k]])
        } 
        else {
          Y[[k]] <- sapply(seq_len(p[i.p]), function(i) {
            x <- Y[[k]][, i]
            id <- sample(n, floor(n * perc[i]))
            x[id] <- NA
            x
          })
          Z[[k]] <- datacjggm(Y[[k]], X = X[[k]])
          Y[[k]][is.na(Y[[k]])] <- 40L
        }
      }
      
      ###################
      # section jcglasso #
      ###################
      dim(Thth) <- Tht.dim
      dim(Ththx) <- Thtx.dim
      dim(Bh) <- B.dim
      
      out <- jcglasso(data = Z, nrho = 1L, nlambda = 1L, thr.bcd = 1L, thr.em = 1L, 
                      trace = 0L, penalty = "group", alpha = alpha)
      # rho1 <- max(out$rho) * c2 / out$alpha
      # lambda1 <- max(out$lambda) * c2 / out$alpha

      rho1 <- (max(out$rho) / out$alpha) * c2
      lambda1 <- (max(out$lambda) / out$alpha) * c2
      out <- jcglasso(data = Z, rho = rho1, lambda = lambda1, trace = 0L, 
                      penalty = "group", alpha = alpha, thr.em = 1e-2)
      
      for(o in seq_len(nrho1)){
        for(i in seq_len(nlambda1)) {
          Thth[, , , i, o] <- list2array(coef(out, type = "Theta", rho.id = o, lambda.id = i))
          mse_Tht[jj, , i, o, "jcglasso"] <- sapply(seq_len(K), function(k) sum((Thth[, , k, i, o][U1] - Theta.y[[k]][U1])^2))
          Ththx[, , , i, o] <- list2array(coef(out, type = "Omega", rho.id = o, lambda.id = i))
          mse_Thtx[jj, , i, o, "jcglasso"] <- sapply(seq_len(K), function(k) sum((Ththx[, , k, i, o][U2] - Theta.x[[k]][U2])^2))
          Bh[, , , i, o] <- list2array(coef(out, type = "B", rho.id = o, lambda.id = i))
          mse_B[jj, , i, o, "jcglasso"] <- sapply(seq_len(K), function(k) sum((Bh[, , k, i, o] - B[[k]])^2))
        }
      }
      
      for(i in seq_len(nlambda1)){
        #################################
        # Precision-Recall Curve (THETA)
        
        for(k in seq_len(K)) {
          # jcglasso
          out_prcurve <- prcurveTheta(thetah = Thth[, , k, i, ], thetat = Theta.y[[k]], nrho = nrho1)
          prcurve_Tht[jj, k, i, , "recall", "jcglasso"] <- out_prcurve$recall
          prcurve_Tht[jj, k, i, , "precision", "jcglasso"] <- out_prcurve$precision
          auc_tbl_Tht[jj, "jcglasso", k, i] <- out_prcurve$auc
          
          out_prcurve <- prcurveTheta(thetah = Ththx[, , k, i, ], thetat = Theta.x[[k]], nrho = nrho1)
          prcurve_Thtx[jj, k, i, , "recall", "jcglasso"] <- out_prcurve$recall
          prcurve_Thtx[jj, k, i, , "precision", "jcglasso"] <- out_prcurve$precision
          auc_tbl_Thtx[jj, "jcglasso", k, i] <- out_prcurve$auc
        }
      }
      
      for(o in seq_len(nrho1)){
        #################################
        # Precision-Recall Curve (B)
        
        for(k in seq_len(K)) {
          # jcglasso
          out_prcurve <- prcurveB(thetah = Bh[, , k, , o], thetat = B[[k]], nrho = nlambda1)
          prcurve_B[jj, k, , o, "recall", "jcglasso"] <- out_prcurve$recall
          prcurve_B[jj, k, , o, "precision", "jcglasso"] <- out_prcurve$precision
          auc_tbl_B[jj, "jcglasso", k, o] <- out_prcurve$auc
        }
      }
      
      dim(Thth) <- NULL
      Tht.array[, "jcglasso", j] <- Thth
      dim(Ththx) <- NULL
      Thtx.array[, "jcglasso", j] <- Ththx
      dim(Bh) <- NULL
      B.array[, "jcglasso", j] <- Bh
    }
    
    save(file = filenames[h], list = c("Tht.dim", "Theta.y", "Tht.array", "mse_Tht", "prcurve_Tht", "auc_tbl_Tht",
                                       "Thtx.dim", "Theta.x", "Thtx.array", "mse_Thtx", "prcurve_Thtx", "auc_tbl_Thtx",
                                       "B.dim", "B", "B.array", "mse_B", "prcurve_B", "auc_tbl_B", 
                                       "nruns", "subnsim", "nsim", "nrho1", "nlambda1", "c2", "filenames"))
  }
  
  close(pb)
  
  for(i in seq_len(nlambda1)){
    mse_tbl_Tht["mean", "jcglasso", i, ] <- colMeans(sapply(seq_len(K), function(k) apply(mse_Tht[, k, i, 1:length(c2), "jcglasso"], 1L, min)))
    mse_tbl_Tht["sd", "jcglasso", i, ] <- apply(sapply(seq_len(K), function(k) apply(mse_Tht[, k, i, 1:length(c2), "jcglasso"], 1L, min)), 2L, sd)
    
    mse_tbl_Thtx["mean", "jcglasso", i, ] <- colMeans(sapply(seq_len(K), function(k) apply(mse_Thtx[, k, i, 1:length(c2), "jcglasso"], 1L, min)))
    mse_tbl_Thtx["sd", "jcglasso", i, ] <- apply(sapply(seq_len(K), function(k) apply(mse_Thtx[, k, i, 1:length(c2), "jcglasso"], 1L, min)), 2L, sd)
  }
  
  print(ftab <- ftable(mse_tbl_Tht, row.vars = c("lambda1", "K"), col.vars = c("estimator", "type")))
  print(ftab2 <- ftable(mse_tbl_Thtx, row.vars = c("lambda1", "K"), col.vars = c("estimator", "type")))
  
  save(file = paste0(percorso, "/mse_tbl_Tht.RData"), list = c("mse_tbl_Tht"))
  save(file = paste0(percorso, "/mse_tbl_Thtx.RData"), list = c("mse_tbl_Thtx"))
  
  for(o in seq_len(nrho1)){
    mse_tbl_B["mean", "jcglasso", o, ] <- colMeans(sapply(seq_len(K), function(k) apply(mse_B[, k, 1:length(c2), o, "jcglasso"], 1L, min)))
    mse_tbl_B["sd", "jcglasso", o, ] <- apply(sapply(seq_len(K), function(k) apply(mse_B[, k, 1:length(c2), o, "jcglasso"], 1L, min)), 2, sd)
  }
  
  print(ftab3 <- ftable(mse_tbl_B, row.vars = c("rho1", "K"), col.vars = c("estimator", "type")))
  
  save(file = paste0(percorso, "/mse_tbl_B.RData"), list = c("mse_tbl_B"))
  
  
  for(i in seq_len(nlambda1)){
    prcurve_data_Tht[i, 1:length(c2), , "recall", "jcglasso"] <- sapply(seq_len(K), function(k) apply(prcurve_Tht[, k, i, 1:length(c2), "recall", "jcglasso"], 2, mean, na.rm=TRUE))
    prcurve_data_Tht[i, 1:length(c2), , "precision", "jcglasso"] <- sapply(seq_len(K), function(k) apply(prcurve_Tht[, k, i, 1:length(c2), "precision", "jcglasso"], 2, mean, na.rm=TRUE))
    
    prcurve_data_Thtx[i, 1:length(c2), , "recall", "jcglasso"] <- sapply(seq_len(K), function(k) apply(prcurve_Thtx[, k, i, 1:length(c2), "recall", "jcglasso"], 2, mean, na.rm=TRUE))
    prcurve_data_Thtx[i, 1:length(c2), , "precision", "jcglasso"] <- sapply(seq_len(K), function(k) apply(prcurve_Thtx[, k, i, 1:length(c2), "precision", "jcglasso"], 2, mean, na.rm=TRUE))
  }
  
  auc_Tht["mean", "jcglasso", , ] <- sapply(seq_len(K), function(k) apply(auc_tbl_Tht[, "jcglasso", k, ], 2, mean))
  auc_Tht["sd", "jcglasso", ,] <- sapply(seq_len(K), function(k) apply(auc_tbl_Tht[, "jcglasso", k, ], 2, sd))
  
  auc_Thtx["mean", "jcglasso", , ] <- sapply(seq_len(K), function(k) apply(auc_tbl_Thtx[, "jcglasso", k, ], 2, mean))
  auc_Thtx["sd", "jcglasso", ,] <- sapply(seq_len(K), function(k) apply(auc_tbl_Thtx[, "jcglasso", k, ], 2, sd))
  
  print(ftab.2 <- ftable(auc_Tht, row.vars = c("lambda1", "K"), col.vars = c("estimator", "type")))
  print(ftab2.2 <- ftable(auc_Thtx, row.vars = c("lambda1", "K"), col.vars = c("estimator", "type")))
  
  save(file = paste0(percorso, "/auc_tbl_Tht.RData"), list = c("auc_tbl_Tht"))
  save(file = paste0(percorso, "/auc_tbl_Thtx.RData"), list = c("auc_tbl_Thtx"))
  
  for(o in seq_len(nrho1)){
    prcurve_data_B[1:length(c2), o, , "recall", "jcglasso"] <- sapply(seq_len(K), function(k) apply(prcurve_B[, k, 1:length(c2), o, "recall", "jcglasso"], 2, mean, na.rm=TRUE))
    prcurve_data_B[1:length(c2), o, , "precision", "jcglasso"] <- sapply(seq_len(K), function(k) apply(prcurve_B[, k, 1:length(c2), o, "precision", "jcglasso"], 2, mean, na.rm=TRUE))
  }
  
  auc_B["mean", "jcglasso", , ] <- sapply(seq_len(K), function(k) apply(auc_tbl_B[, "jcglasso", k, ], 2, mean))
  auc_B["sd", "jcglasso", ,] <- sapply(seq_len(K), function(k) apply(auc_tbl_B[, "jcglasso", k, ], 2, sd))
  
  print(ftab3.3 <- ftable(auc_B, row.vars = c("rho1", "K"), col.vars = c("estimator", "type")))
  
  save(file = paste0(percorso, "/auc_tbl_B.RData"), list = c("auc_tbl_B"))
  
  save.image(file = paste0(percorso, "/analisi.RData"))
  
# }