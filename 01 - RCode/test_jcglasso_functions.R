#### a simple R script to learn how to use jcglasso main class and S3 methods

library(cglasso)
source("datajcggm.R")
source("jcglasso.R")
source("jcggm.R")
source("gof.R")
source("to_graph.R")

library(MASS)
library(huge)

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
  cutoff <- function(b0, i, up, perc = perc) (mean(pnorm(up - b0 - eta[, i], lower.tail =  FALSE))) - perc
  B[1L, ] <- sapply(1:p, function(.i) uniroot(cutoff, interval = c(0, 100), i = .i, up = up, perc = perc[.i])$root)
  
  if(!is.na(unconnectedB)) {
    tmp <- which(B[-1L, ] != 0)
    id <- sample(tmp, size = floor(unconnectedB * length(tmp)), replace = FALSE)
    B[-1L, ][id] <- 0.0
  }
  
  return(B)
}
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

n <- 300L                # sample size (integer)
p <- 25L                 # number of response variables (integer)
q <- 10L                 # number of predictors (integer)
K <- 2L                  # number of subpopulation (integer)
ncompPrecision <- 5L     # number of connected components in the precision matrices Omega e Theta for each k (integer)
ncompB <- 2L             # number of connected components in the coefficient regression matrix B for each k (integer)
percConnectedCompK <- .5 # percentages of connected components for Omega, Theta and B across k (numeric vector of length k - 1)
perc.NA <- 0.25          # percentage of missing-at-random data for each column
perc.X.na <- 0.25        # percentage of X columns with missing-at-random data
perc.Y.na <- 0.25        # percentage of Y columns with missing-at-random data
up <- 40.0               # upper censored value for the response matrix
perc.cens <- 0.25        # percentage of censored data for each column
perc.Y.cens <- 0.25      # percentage of Y columns with censored data
Thtmin <- 0.3            # minimum value of the precision matrix Theta
Thtmax <- 0.5            # maximum value of the precision matrix Theta
Omgmin <- 0.3            # minimum value of the precision matrix Omega
Omgmax <- 0.5            # maximum value of the precision matrix Omega
Bmin <- 0.3              # minimum value of coefficient regression matrix B
Bmax <- 0.7              # maximum value of coefficient regression matrix B

sim <- rjcggm(n = n, p = p, q = q, K = K, ncompPrecision = ncompPrecision, ncompB = ncompB, 
              percConnectedCompK = percConnectedCompK, perc.NA = perc.NA, perc.X.na = perc.X.na, perc.Y.na = perc.Y.na, 
              up = up, perc.cens = perc.cens, perc.Y.cens = perc.Y.cens, Thtmin = Thtmin, Thtmax = Thtmax, 
              Omgmin = Omgmin, Omgmax = Omgmax, Bmin = Bmin, Bmax = Bmax)

##### datajcggm object with only the X matrix
sim$ZX

##### datajcggm object with both the X and Y matrices
sim$Z

##### produces summaries of an object of class ‘datajcggm’ for the first 
summary(sim$Z)

##### fits the joint graphical lasso model to X datasets with censored and/or missing values
model1 <- jcglasso(data = sim$ZX, nrho = 25L, rho.min.ratio = .01, alpha1 = .5, trace = 1L)
model1

##### choice of the optimal nu parameter fixed alpha = 0.5
par(mfrow = c(2, 2))
plot(AIC(model1))
plot(BIC(model1))

qfun <- QFun2(model1, mle = TRUE)
plot(AIC(model1, Qfun = qfun))
plot(BIC(model1, Qfun = qfun))

##### best model obtained minimizing the BIC criterion
summary(model1, GoF = BIC(model1, Qfun = qfun), print.all = FALSE)
model1best <- select.jcglasso(model1, GoF = BIC(model1, Qfun = qfun))
model1best

##### simple grayscale plot of the precision matrices Omega
plot(model1best)

##### fits the joint conditional graphical lasso model to the datasets with censored and/or missing values
##### fixing nu to its optimal value and alpha3 = 0.5
model2 <- jcglasso(data = sim$Z, nrho = 10L, rho.min.ratio = .01, alpha1 = .5,
                   nlambda = 10L, lambda.min.ratio = .01, alpha2 = .5,
                   nu = model1best$rho, alpha3 = model1best$alpha1, trace = 1L)

##### check for KKT conditions
gradB(model2, lambda.id = 5, rho.id = 5)$maxEps
gradTht(model2, lambda.id = 5, rho.id = 5)$maxEps

##### choice of the optimal lambda and rho parameters fixed alpha1 = alpha2 = 0.5
##### and fixing nu to its optimal value and alpha3 = 0.5
par(mfrow = c(2, 2))
plot(AIC(model2))
plot(BIC(model2))

qfun2 <- QFun2(model2, mle = TRUE)
plot(AIC(model2, Qfun = qfun2))
plot(BIC(model2, Qfun = qfun2))

##### best model obtained minimizing the BIC criterion
summary(model2, GoF = BIC(model2, Qfun = qfun2), print.all = FALSE)
model2best <- select.jcglasso(model2, GoF = BIC(model2, Qfun = qfun2))
model2best

##### check for KKT conditions
gradB(model2best, lambda.id = 1, rho.id = 1)$maxEps
gradTht(model2best, lambda.id = 1, rho.id = 1)$maxEps

##### simple grayscale plot of the precision matrices Theta
plot(model2best)

#### perform post-hoc maximum likelihood refitting of a selected joint 
#### conditional graphical lasso model with censored and/or missing values.
model2bestMLE <- jcggm(model2best)
model2bestMLE

coef(model2bestMLE, type = "B")
coef(model2bestMLE, type = "Theta")
coef(model2bestMLE, type = "Omega")

##### returns a named list of graphs using the results of an R object of class ‘jcglasso’ or ‘jcggm’.

graphs <- to_graph2(model2bestMLE)
graphs

par(mfrow = c(1, 2))
# plot of an undirected graph representing the conditional dependence structure among the p response variables (i.e., Theta)
plot(graphs, type = "Gyy")
# plot of a directed graph representing the effetcs of the q predictors onto the p response variables (i.e., B)
plot(graphs, type = "Gxy")
# plot of an undirected graph representing the conditional dependence structure among the q covariates (i.e., Omega)
plot(graphs, type = "Gxx")
# plot of both the undirected graph representing the conditional dependence structure among the p response variables (i.e.,Theta) 
# and the directed graph representing the effects of the q predictors on the p response variables (i.e., B).
plot(graphs, type = "conditional")
# overall plot of both the undirected graph representing the conditional dependence structure among the p response variables (i.e., Theta) 
# and the undirected graph representing the conditional dependence structure among the q covariates (i.e., Omega), 
# plus the directed graph representing the effects of the q predictors on the p response variables (i.e., B).
plot(graphs, type = "bipartite")
