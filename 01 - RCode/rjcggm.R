#### main function to generate simulated data
rjcggm <- function(n = 300L, p = 25L, q = 10L, K = 2L, ncompPrecision = 5L, ncompB = 2L, 
                   percConnectedCompK = .5, perc.NA = .25, perc.X.na = .25, perc.Y.na = .25, 
                   up = 40.0, perc.cens = .25, perc.Y.cens = .25, Thtmin = .3, Thtmax = .5, 
                   Omgmin = .3, Omgmax = .5, Bmin = .5, Bmax = .7) {
  
  perc.cens <- pmax(perc.cens, 1E-6)
  
  indexB <- sapply(seq_len(p), \(i) sort(sample(q, size = ncompB, replace = FALSE)))
  unconnectedB <- c(NA, percConnectedCompK)
  unconnectedTheta <- c(NA, percConnectedCompK)
  unconnectedOmg <- c(NA, percConnectedCompK)
  
  Theta.x <- Sigma.x <- X <- B <- Theta.y <- Sigma.y <- Y <-  Xorig <- Yorig <- 
    ZX <- Z <- vector(mode = "list", length = K)
  for(k in 1:K) {
    Theta.x[[k]] <- simulOmega(p = q, unconnected = unconnectedOmg[k],
                               ncomp = ncompPrecision, mintht = Omgmin, maxtht = Omgmax)
    Sigma.x[[k]] <- solve(Theta.x[[k]])
    
    X[[k]] <- mvrnorm(n = n, mu = rep(0.0, q), Sigma = Sigma.x[[k]])
    perc <- rep(1E-6, p)
    perc[sample(p, floor(p * perc.Y.cens))] <- perc.cens
    B[[k]] <- simulB(X[[k]], perc, up, q, p, Bmin, Bmax, indexB, unconnectedB[k])
    eta <- cbind(1.0, X[[k]]) %*% B[[k]]
    mu <- eta
    
    Theta.y[[k]] <- simulTheta(p = p, unconnected = unconnectedTheta[k],
                               ncomp = ncompPrecision, mintht = Thtmin, maxtht = Thtmax)
    Sigma.y[[k]] <- solve(Theta.y[[k]])
    
    E <- mvrnorm(n = n, mu = rep(0, p), Sigma = Sigma.y[[k]])
    Y[[k]] <- mu + E
    
    Xorig[[k]] <- X[[k]]
    X[[k]][sample(n, floor(n * perc.NA)), sample(q, floor(q * perc.X.na))] <- NA
    Yorig[[k]] <- Y[[k]]
    Y[[k]][sample(n, floor(n * perc.NA)), sample(p, floor(p * perc.Y.na))] <- NA
  }
  ZX <- datajcggm(Y = X, up = up) 
  Z <- datajcggm(Y = Y, X = X, up = up)
  
  list(ZX = ZX, Z = Z, originals = list(X = Xorig, B = B, Y = Yorig))
}


#### auxiliary functions to simulate block sparse structures
simulOmega <- function(p, ncomp, unconnected = NA, mintht = 0.3, maxtht = 0.5){
  if(p %% ncomp != 0) warning("problem in dimensions p and ncomp")
  frac <- p / ncomp
  block <- matrix(0.0, ncomp, ncomp)
  block[1L, ] <- block[, 1L] <- sort(runif(ncomp, min = mintht, maxtht), TRUE)
  diag(block) <- 1.0
  subTht <- diag(frac)
  Tht <- kronecker(subTht, block)
  
  if(!is.na(unconnected)) {
    tmp <- which(Tht[upper.tri(Tht)] != 0)
    id <- sample(tmp, size = floor(unconnected * length(tmp)))
    Tht[upper.tri(Tht)][id] <- 0.0
    Tht <- t(Tht)
    Tht[upper.tri(Tht)][id] <- 0.0
  }
  
  diag(Tht) <- 1.0
  Tht
}
simulTheta <- function(p, ncomp, unconnected = NA, mintht = 0.3, maxtht = 0.5){
  if(p %% ncomp != 0) warning("problem in dimensions p and ncomp")
  frac <- p / ncomp
  block <- matrix(0.0, ncomp, ncomp)
  block[1L, ] <- block[, 1L] <- sort(runif(ncomp, min = mintht, maxtht), TRUE)
  diag(block) <- 1.0
  subTht <- diag(frac)
  Tht <- kronecker(subTht, block)
  
  if(!is.na(unconnected)) {
    tmp <- which(Tht[upper.tri(Tht)] != 0)
    id <- sample(tmp, size = floor(unconnected * length(tmp)))
    Tht[upper.tri(Tht)][id] <- 0.0
    Tht <- t(Tht)
    Tht[upper.tri(Tht)][id] <- 0.0
  }
  
  diag(Tht) <- 1.0
  Tht
}
simulB <- function(X, perc, up, q, p, Bmin = .5, Bmax = .7, indexB, unconnectedB = NA) {
  B <- matrix(0.0, nrow = q + 1L, ncol = p)
  for(j in seq_len(p)) B[1L + indexB[, j], j] <- runif(2, Bmin, Bmax)
  
  eta <- X %*% B[-1L, ]
  if(up > 0) {
    cutoff <- function(b0, i, up, perc = perc) (mean(pnorm(up - b0 - eta[, i], lower.tail =  FALSE))) - perc
    B[1L, ] <- sapply(1:p, function(.i) uniroot(cutoff, interval = c(0, 100), i = .i, up = up, perc = perc[.i])$root)
  }
  
  if(!is.na(unconnectedB)) {
    tmp <- which(B[-1L, ] != 0)
    id <- sample(tmp, size = floor(unconnectedB * length(tmp)), replace = FALSE)
    B[-1L, ][id] <- 0.0
  }
  
  return(B)
}
