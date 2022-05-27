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

# for(.id in c(1, 2, 5, 6, 3, 4, 7, 8)) {
  .id <- 8
  
  # clean the workspace
  rm(list = setdiff(ls(), "index"))
  
  i.H <- index[.id, 1]
  i.p <- index[.id, 2]
  i.perc <- index[.id, 3]
  
  model <- "censoring"
  
  alpha <- .5
  
  percorso <- paste0("results/")
  if(!dir.exists(percorso)) dir.create(percorso)
  
  # libraries
  library(MASS)
  library(cglasso)
  library(JGL)
  source('../../../../01 - RCode/jcglasso.R')
  
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
  K <- 3                          # number of groups
  pH <- c(0.20, 0.40)
  perc <- c(0.20, 0.40)
  unconnected <- list(NULL, 9, 9:10)
  c2 <- c(1, .75, .5, .25, .1)
  
  # rho / rho_max = c(1, .75, .5, .25, .1, .01) => rho = rho_max * c(1, .75, .5, .25, .1, .01)
  
  H <- p[i.p] * pH[i.H]
  nrho1 <- length(c2)                              # number of rho1-values
  rho1.min.ratio <- 1E-2
  nrho2 <- 1L                                      # number of rho2-values
  rho2.min.ratio <- 1L
  perc <- rep(c(perc[i.perc], 1e-6), c(H, p[i.p] - H))      # percentage of expected missing
  
  ##############
  # output
  ##############
  conv <- matrix("", nsim, 2)
  colnames(conv) <- c("jgl", "jcglasso")
  
  Tht.dim <- c(p[i.p], p[i.p], K, nrho2, nrho1)
  Thth <- array(0, dim = Tht.dim)
  Tht.array <- array(0, dim = c(prod(Tht.dim), 2, subnsim), 
                     dimnames = list(NULL, c("jgl", "jcglasso"), NULL) )
  
  filenames <- paste0(percorso, "run.", seq_len(nruns), ".RData")
  
  mse_Tht <- array(0, dim = c(nsim, K, nrho2, nrho1, 2), 
                   dimnames = list(nsim = NULL, K = paste0("group", seq_len(K)), 
                                   nrho2 = 1:nrho2, nrho1 = 1:nrho1, 
                                   estimator = c("jcglasso", "jgl")))
  
  mse_tbl_Tht <- array(0, dim = c(2, 2, nrho1, K), 
                       dimnames = list(type = c("mean", "sd"), 
                                       estimator = c("jcglasso", "jgl"), 
                                       rho1 = c2, K = paste0("group", seq_len(K))))
  
  #############################################
  # starting simulation study
  #############################################
  set.seed(123)
  
  Theta.y <- Sigma.y <- B <- mu <- E <- Y <- S <- Z <- vector(mode = "list", length = K)
  eta <- matrix(0, n, p[i.p])
  X <- cbind(rep(1, n))
  
  for(k in seq_len(K)) {
    Theta.y[[k]] <- simulTheta(p = p[i.p], unconnected = unconnected[[k]], 
                               ncomp = 5, mintht = 0.3, maxtht = 0.5)
    Sigma.y[[k]] <- solve(Theta.y[[k]])
    
    B[[k]] <- matrix(0, nrow = 1L, ncol = p[i.p])
    
    cutoff <- function(b0, i, up, perc = perc) (mean(pnorm(up - b0 - eta[, i], lower.tail =  FALSE))) - perc
    B[[k]][1L, ] <- if(model == "censoring") sapply(1:p[i.p], function(.i) uniroot(cutoff, interval = c(0, 100), i = .i, up = up, perc = perc[.i])$root)
    else sapply(1:p[i.p], function(.i) uniroot(cutoff, interval = c(0, 100), i = .i, up = up, perc = 1E-6)$root)
    
    mu[[k]] <- X %*% B[[k]]
  }
  
  pb <- txtProgressBar(min = 0L, max = nsim, style = 3L)
  
  jj <- 0L
  
  U <- outer(1:p[i.p], 1:p[i.p], "<")
  U1 <- outer(1:p[i.p], 1:p[i.p], "<=")
  
  for(h in seq_len(nruns)) {
    for(j in seq_len(subnsim)) {
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
          Z[[k]] <- datacjggm(Y[[k]], up = up)
        } 
        else {
          Y[[k]] <- sapply(seq_len(p[i.p]), function(i) {
            x <- Y[[k]][, i]
            id <- sample(n, floor(n * perc[i]))
            x[id] <- NA
            x
          })
          Z[[k]] <- datacjggm(Y[[k]])
          Y[[k]][is.na(Y[[k]])] <- up
        }
      }
      
      # rho1_max <- (1 / K) * max(sapply(S, function(x) max(abs(x[U]))))
      # rho1_min <- rho1.min.ratio * rho1_max
      # rho1 <- exp(seq(from = log(rho1_max), to = log(rho1_min), length = nrho1))
      # 
      # rho2_max <- rho1_max
      # rho2_min <- rho2.min.ratio * rho2_max
      # rho2 <- rho2_max * c2
      # rho2 <- exp(seq(from = log(rho2_max), to = log(rho2_min), length = nrho2))
      
      ###################
      # section jcglasso #
      ###################
      dim(Thth) <- Tht.dim
      # out <- jcglasso(data = Z, nrho1 = nrho1, rho1.min.ratio = rho1.min.ratio, nrho2 = nrho2, rho2.min.ratio = rho2.min.ratio, 
      #                 trace = 1L, penalty = "group")
      out <- jcglasso(data = Z, nrho = 1L, thr.bcd = 1L, thr.em = 1L, trace = 0L, penalty = "group", alpha = alpha)
      rho1 <- max(out$rho) * c2 / out$alpha
      
      out <- jcglasso(data = Z, rho = rho1, trace = 0L, penalty = "group", alpha = alpha, thr.em = 1E-2)
      i <- 1L
      for(o in seq_len(nrho1)) {
        #   # for(i in seq_len(nrho2)) {
        #     # out <- jcglasso(data = Z, penalty = "group", rho1 = rho1[o], rho2 = rho2[i], trace = 0L)
        #     # Thth[, , , i, o] <- list2array(out$Tht)
        #     # mse_Tht[jj, , i, o, "jcglasso"] <- sapply(seq_len(K), function(k) sum((Thth[, , k, i, o][U1] - Theta.y[[k]][U1])^2))
        Thth[, , , i, o] <- list2array(lapply(out$Tht, function(x) x[, , i, o]))
        mse_Tht[jj, , i, o, "jcglasso"] <- sapply(seq_len(K), function(k) sum((Thth[, , k, i, o][U1] - Theta.y[[k]][U1])^2))
        #   # }
      }
      
      dim(Thth) <- NULL
      Tht.array[, "jcglasso", j] <- Thth
      
      rho1 <- out$rho
      rho2 <- out$rho2
      
      ##################
      # section JGL #
      ##################
      dim(Thth) <- Tht.dim
      for(o in seq_len(nrho1)) {
        # for(i in seq_len(nrho2)) {
        out <- JGL(Y = Y, penalty = "group", lambda1 = rho1[o], lambda2 = rho2[o], tol = 1E-2, 
                   return.whole.theta = TRUE, truncate = 1E-6, weights = "sample.size")
        Thth[, , , i, o] <- list2array(out$theta)
        mse_Tht[jj, , i, o, "jgl"] <- sapply(seq_len(K), function(k) sum((Thth[, , k, i, o][U1] - Theta.y[[k]][U1])^2))
        # }
      }

      dim(Thth) <- NULL
      Tht.array[, "jgl", j] <- Thth
    }
    
    save(file = filenames[h], list = c("Tht.dim", "Theta.y", "Tht.array", "nruns", 
                                       "subnsim", "nsim", "nrho1", "nrho2", "filenames"))
  }
  
  close(pb)
  
  mse_tbl_Tht["mean", "jgl", , ] <- sapply(seq_len(K), function(k) apply(mse_Tht[, k, 1L, , "jgl"], 2L, mean))
  mse_tbl_Tht["sd", "jgl", , ] <- sapply(seq_len(K), function(k) apply(mse_Tht[, k, 1L, , "jgl"], 2L, sd))
  mse_tbl_Tht["mean", "jcglasso", , ] <- sapply(seq_len(K), function(k) apply(mse_Tht[, k, 1L, , "jcglasso"], 2L, mean))
  mse_tbl_Tht["sd", "jcglasso", , ] <- sapply(seq_len(K), function(k) apply(mse_Tht[, k, 1L, , "jcglasso"], 2L, sd))
  
  ftab <- ftable(mse_tbl_Tht, row.vars = c("rho1", "K"), col.vars = c("estimator", "type"))
  
  print(ftab)
  
  save(file = paste0(percorso, "mse_tbl.RData"), list = c("mse_tbl_Tht"))
  
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
    out
  }
  
  ###################
  # output
  ###################
  
  prcurve_Tht <- array(0, dim = c(nsim, K, nrho2, nrho1, 2, 2), 
                       dimnames = list(nsim = NULL, K = paste0("group", seq_len(K)), 
                                       nrho2 = 1:nrho2, nrho1 = 1:nrho1, 
                                       type = c("recall", "precision"), 
                                       estimator = c("jcglasso", "jgl")))
  
  auc_tbl_Tht <- array(0, dim = c(nsim, 2, K, nrho2), 
                       dimnames = list(nsim = NULL, estimator = c("jcglasso", "jgl"), 
                                       K = paste0("group", seq_len(K)), rho2 = 1:nrho2))
  
  prcurve_data_Tht <- array(0, dim = c(nrho2, nrho1, K, 2, 2), 
                            dimnames = list(rho2 = 1:nrho2, rho1 = 1:nrho1, 
                                            K = paste0("group", seq_len(K)), 
                                            type = c("recall", "precision"), 
                                            estimator = c("jcglasso", "jgl")))
  
  auc_Tht  <- array(0, dim = c(2, 2, nrho2, K), 
                    dimnames = list(type = c("mean", "sd"), 
                                    estimator = c("jcglasso", "jgl"), 
                                    rho2 = NULL, K = paste0("group", seq_len(K))))
  
  pb <- txtProgressBar(min = 0, max = nsim, style = 3)
  
  iii <- 0
  for(h in seq_len(nruns)) {
    
    load(filenames[h])
    U1 <- upper.tri(Theta.y, diag = TRUE)
    
    for(ii in seq_len(subnsim)) {
      iii <- iii + 1L
      setTxtProgressBar(pb, iii)
      
      ##################################
      # Analisys: Precision matrix
      ##################################
      
      Thth.jgl <- Tht.array[, "jgl", ii]
      Thth.jcglasso <- Tht.array[, "jcglasso", ii]
      dim(Thth.jgl) <- dim(Thth.jcglasso) <- Tht.dim
      # rho.mat <- rho.array[, , ii]
      
      for(j in seq_len(nrho2)){
        #################################
        # Precision-Recall Curve
        
        # jgl
        for(k in seq_len(K)) {
          out_prcurve <- prcurveTheta(thetah = Thth.jgl[, , k, j, ], thetat = Theta.y[[k]], nrho = nrho1)
          prcurve_Tht[iii, k, j, , "recall", "jgl"] <- out_prcurve$recall
          prcurve_Tht[iii, k, j, , "precision", "jgl"] <- out_prcurve$precision
          auc_tbl_Tht[iii, "jgl", k, j] <- out_prcurve$auc
          
          # jcglasso
          out_prcurve <- prcurveTheta(thetah = Thth.jcglasso[, , k, j, ], thetat = Theta.y[[k]], nrho = nrho1)
          prcurve_Tht[iii, k, j, , "recall", "jcglasso"] <- out_prcurve$recall
          prcurve_Tht[iii, k, j, , "precision", "jcglasso"] <- out_prcurve$precision
          auc_tbl_Tht[iii, "jcglasso", k, j] <- out_prcurve$auc
        }
      }
    }
  }
  
  ##################################
  # Analisys: Precision matrix
  ##################################
  
  for(j in seq_len(nrho2)){
    prcurve_data_Tht[j, , , "recall", "jgl"] <- sapply(seq_len(K), function(k) apply(prcurve_Tht[, k, j, , "recall", "jgl"], 2, mean, na.rm=TRUE))
    prcurve_data_Tht[j, , , "precision", "jgl"] <- sapply(seq_len(K), function(k) apply(prcurve_Tht[, k, j, , "precision", "jgl"], 2, mean, na.rm=TRUE))
    
    prcurve_data_Tht[j, , , "recall", "jcglasso"] <- sapply(seq_len(K), function(k) apply(prcurve_Tht[, k, j, , "recall", "jcglasso"], 2, mean, na.rm=TRUE))
    prcurve_data_Tht[j, , , "precision", "jcglasso"] <- sapply(seq_len(K), function(k) apply(prcurve_Tht[, k, j, , "precision", "jcglasso"], 2, mean, na.rm=TRUE))
  }
  
  auc_Tht["mean", "jcglasso", , ] <- sapply(seq_len(K), function(k) apply(auc_tbl_Tht[, "jcglasso", k, , drop=FALSE], 2, mean))
  auc_Tht["sd", "jcglasso", ,] <- sapply(seq_len(K), function(k) apply(auc_tbl_Tht[, "jcglasso", k, , drop=FALSE], 2, sd))
  
  auc_Tht["mean", "jgl", ,] <- sapply(seq_len(K), function(k) apply(auc_tbl_Tht[, "jgl", k, , drop=FALSE], 2, mean))
  auc_Tht["sd", "jgl", ,] <- sapply(seq_len(K), function(k) apply(auc_tbl_Tht[, "jgl", k, , drop=FALSE], 2, sd))
  
  save(file = paste0(percorso, "auc_tbl.RData"), list = c("auc_Tht"))
  
  # table and figure
  
  library(xtable)
  
  print(ftab, digits = 2)
  
  xtable(as.matrix(ftab))
  
  ftab2 <- ftable(auc_Tht, row.vars = c("rho2", "K"), col.vars = c("estimator", "type"))
  
  print(ftab2, digits = 2)
  
  xtable(as.matrix(ftab2))
  
  save.image(file = paste0(percorso, "analisi.RData"))
  
# }
