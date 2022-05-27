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


# clean the workspace
rm(list = ls())
gc(reset = TRUE)

i.H <- 2L
i.p <- 2L
i.q <- 1L
i.perc <- 2L

model <- "censoring"

alpha <- .5

percorso <- paste0("results/")
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

nruns <- 1
subnsim <- 10
nsim <- nruns * subnsim
up <- 40                        # right censoring value
n <- 100                        # sample size
p <- c(50, 200)                 # number of response variables
q <- c(50, 200)                 # number of predictors
K <- 3                          # number of groups
pH <- c(0.20, 0.40)
perc0 <- c(0.20, 0.40)
unconnected <- list(NULL, 9, 9:10)
c2 <- c(1, .75, .5, .25, .1) #seq(1, .1, l = 10) #c(1, .75, .5, .25, .1)
Bmin <- 0.3
Bmax <- 0.7

# rho / rho_max = c(1, .75, .5, .25, .1, .01) => rho = rho_max * c(1, .75, .5, .25, .1, .01)

H <- p[i.p] * pH[i.H]
H2 <- q[i.q] * pH[i.H]
nrho1 <- length(c2)                              # number of rho1-values
rho1.min.ratio <- 1
nlambda1 <- length(c2)                           # number of lambda1-values
lambda1.min.ratio <- 1
nnu1 <- length(c2)                           # number of nu1-values
nu1.min.ratio <- 1

perc <- rep(c(perc0[i.perc], 1e-6), c(H, p[i.p] - H))         # percentage of expected missing for Y
perc2 <- rep(c(perc0[i.perc], 1e-6), c(H2, q[i.q] - H2))      # percentage of expected missing for X

##############
# output
##############
conv <- matrix("", nsim, 1)
colnames(conv) <- c("jcglasso")

# filenames <- paste0(paste0(percorso, "/run."), seq_len(nruns), ".RData")
filenames <- paste0(percorso, "run.", seq_len(nruns), ".RData")

mean.dim <- c(nsim, p[i.p] + q[i.q], K, nlambda1, nrho1, nnu1)
MeanZ <- array(0, dim = mean.dim)
mseMeanZ <- array(0, dim = c(nsim, p[i.p] + q[i.q], nlambda1, nrho1, nnu1),
                  dimnames = list(1:nsim, c(paste0("Y", 1:p[i.p]), paste0("X", 1:q[i.q])),
                                  paste0("lambda", 1:nlambda1), paste0("rho", 1:nrho1), 
                                  paste0("nu", 1:nnu1)))

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

Z2 <- lapply(X, datacjggm)
out0 <- jcglasso(data = Z2, nrho = 1, thr.em = 1, 
                 trace = 0L, penalty = "group", alpha = alpha)
nu <- (out0$rho + out0$rho2) * c2

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
    true.mean <- sapply(1:K, function(k) colMeans(cbind(Yorig[[k]], Xorig[[k]])))
    
    ####################
    # section jcglasso #
    ####################
    out <- jcglasso(data = Z, nrho = 1L, nlambda = 1L, thr.bcd = 1L, thr.em = 1L, 
                    trace = 0L, penalty = "group", alpha = alpha, rho.x = 0)
    rho1 <- (out$rho + out$rho2) * c2
    lambda1 <- (out$lambda + out$lambda2) * c2
    
    for(o in 1:nnu1){
      out <- jcglasso(data = Z, rho = rho1, lambda = lambda1, trace = 0L, rho.x = nu[o],
                      penalty = "group", alpha = alpha, thr.em = 1e-2)
      for(i in 1:nlambda1){
        for(ii in 1:nrho1){
          MeanZ[jj, , , i, ii, o] <- sapply(1:K, function(k) colMeans(out$Zipt[[k]][, , i, ii]))
          mseMeanZ[jj, , i, ii, o] <- apply(MeanZ[jj, , , i, ii, o] - true.mean, 1, function(x) mean(x^2))
        }
      }
    }
  }
  
  save(file = filenames[h], list = c("MeanZ", "mseMeanZ", "mean.dim",
                                     "nruns", "subnsim", "nsim", 
                                     "nrho1", "nlambda1", "nnu1", "c2", "filenames"))
}

close(pb)

save.image(file = paste0(percorso, "analisi.RData"))

percent <- function(x, digits = 2, format = "f", ...) {
  paste0(formatC(100 * x, format = format, digits = digits, ...), "%")
}

c3 <- c2
c3 <- paste0(formatC(c3, format = "f", digits = 2))

boxplotMeanZ <- data.frame(value = NA, lambda = rep(c3, each = 20*5), rho = rep(rep(c3, each = 20), 5))
for(i in 1:nlambda1){
    boxplotMeanZ$value[((i-1)*100+1):(i*100)] <- c(sapply(1:5, function(o) (rowMeans(sapply(1:nsim, function(.j) rowMeans(sapply(1:K, function(k) apply(MeanZ[.j, 201:220, k, i, o, ], 1, sd))))))))
}
library(ggplot2)
g1 <- ggplot(boxplotMeanZ, aes(x = lambda, y = value, fill = rho)) + 
  geom_boxplot(alpha = 0.75, notch = FALSE, notchwidth = 0.8) + 
  xlab(expression(lambda/lambda[max])) + ylab(expression(paste(max[nu],"SE({", hat(mu),"})"))) + 
  labs(fill = expression(rho/rho[max])) + theme_bw() + ylim(c(0, 0.03)) +
  theme(axis.title.x = element_text(size = 16), 
        axis.title.y = element_text(size = 16, face = "bold"),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 16), 
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16))
g1
ggsave("figs/fig3supp_boxplot_mu.pdf", plot = g1, device = "pdf", units = "in", width = .7*12, height = .7*8)


