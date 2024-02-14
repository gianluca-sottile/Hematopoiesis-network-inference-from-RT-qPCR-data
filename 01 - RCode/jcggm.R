##### jcggm functions, used to fit a post-hoc maximum likelihood refit of a joint conditional graphical lasso estimator with partially observed data #####
jcggm <- function(object, GoF = AIC, rho.id, lambda.id, tp.min = 1.0E-6, trace = 0L, refit = FALSE, ...) {
  if (!inherits(object, "jcglasso")) stop(sQuote("object"), " is not an object of class ", sQuote("jcglasso"))
  # testing 'tp.min'
  if (!is.vector(tp.min)) stop(sQuote("tp.min"), " is not a vector")
  if (length(tp.min) != 1L) stop(sQuote("tp.min"), " is not a vector of length ", sQuote("1"))
  if (tp.min < 0L) stop(sQuote("tp.min"), " is not a positive value")
  # Testing 'trace'
  if (!is.vector(trace)) stop(sQuote("trace"), " is not a vector")
  if (length(trace) != 1) stop(sQuote("trace"), " is not an object of length ", sQuote(1))
  if (is.logical(trace)) stop(sQuote("trace"), " is not an object of type ", dQuote("integer"))
  if (abs(as.integer(trace) - trace) > 0) stop(sQuote("trace"), " is not an object of type ", dQuote("integer"))
  if (!is.element(trace, c(0L, 1L, 2L))) stop("not allowed value in ", sQuote("trace"), ". Please, choice ", sQuote("0"), ", ", sQuote("1"), " or ", sQuote("2"))
  
  tp.min <- max(min(tp.min, min(object$lambda), min(object$rho)), 1e-13)
  
  this.call <- match.call()
  num_class <- length(object$nobs)
  q <- object$npred
  
  if (missing(rho.id) & missing(lambda.id)) {
    if (!is.element(class(GoF), c("function", "GoF2")))
      stop (sQuote("GoF"), " is not either a goodness-of-fit function (AIC or BIC) neither an object of class ", sQuote("GoF"))
    dots <- list(...)
    if (is.function(GoF)) {
      if (is.null(dots$type)) dots$type <- ifelse(q == 0L, "FD", "CC")
      GoF.name <- deparse(substitute(GoF))
      if (!is.element(GoF.name, c("AIC", "BIC")))
        stop(sQuote(GoF.name), " is not a valid function. Please, use ", sQuote("AIC"), " or ", sQuote("BIC"))
      GoF <- switch(GoF.name,
                    AIC = do.call(function(...) AIC(object, ...), dots),
                    BIC = do.call(function(...) BIC(object, ...), dots))
    }
    object <- select.jcglasso(object, GoF = GoF)
  } 
  else {
    # testing 'lambda.id'
    nlambda <- object$nlambda
    if (missing(lambda.id)) {
      if (nlambda == 1L) lambda.id <- 1L
      else stop(sQuote("lambda.id"), " is missing")
    } 
    else {
      if (!is.vector(lambda.id)) stop(sQuote("lambda.id"), " is not a vector")
      if (length(lambda.id) != 1L) stop(sQuote("lambda.id"), " is not a vector of length ", sQuote("1"))
      if (any(abs(as.integer(lambda.id) - lambda.id) > 0)) stop(sQuote("pen2.id"), " is not an object of type ", dQuote("integer"))
      if (lambda.id <= 0) stop(sQuote("pen2.id"), " is not a positive integer")
      if (lambda.id > nlambda) stop("some entry in ", sQuote("lambda.id"), " is larger than ", sQuote(nlambda))
    }
    # testing 'rho.id'
    nrho <- object$nrho
    if (missing(rho.id)) {
      if (nrho == 1L) rho.id <- 1L
      else stop(sQuote("rho.id"), " is missing")
    }
    else {
      if (!is.vector(rho.id)) stop(sQuote("rho.id"), " is not a vector")
      if (length(rho.id) != 1L) stop(sQuote("rho.id"), " is not a vector of length ", sQuote("1"))
      if (any(abs(as.integer(rho.id) - rho.id) > 0)) stop(sQuote("rho.id"), " is not an object of type ", dQuote("integer"))
      if (rho.id <= 0L) stop(sQuote("rho.id"), " is not a positive integer")
      if (rho.id > nrho) stop("some entry in ", sQuote("rho.id"), " is larger than ", sQuote(nrho))
    }
    
    object <- select.jcglasso(object, GoF = GoF, rho.id = rho.id, lambda.id = lambda.id)
  }
  nrho <- 1L
  nlambda <- 1L
  id_X <- object$InfoStructure$id_X
  id_Y <- object$InfoStructure$id_Y
  
  wTht <- object$weights.Tht
  wB <- object$weights.B
  for(i in 1:num_class){
    wTht[id_Y, id_Y, i][object$InfoStructure$Adj[id_Y, id_Y, i, 1L, 1L] == 0] <- .Machine$double.xmax
    diag(wTht[id_Y, id_Y, i]) <- 0.0
    if(q > 0L) {
      wTht[id_X, id_X, i][object$InfoStructure$Adj[id_X, id_X, i, 1L, 1L] == 0] <- .Machine$double.xmax
      if(is.matrix(wTht[id_X, id_X, i])) diag(wTht[id_X, id_X, i]) <- 0.0 else wTht[id_X, id_X, i] <- 0.0
      wB[, , i][object$InfoStructure$Adj[id_X, id_Y, i, 1L, 1L] == 0] <- .Machine$double.xmax
    }
  }
  
  rho.true <- object$rho
  lambda.true <- object$lambda
  alpha1 <- object$alpha1
  alpha2 <- object$alpha2
  alpha3 <- object$alpha3
  model <- paste0("MLE of a ", object$model)
  
  object$weights.B <- wB
  object$weights.Tht <- wTht
  
  if(!refit) {
    out <- jcggm.fit(object, tp.min = tp.min, trace = trace, ...)
  }
  else {
    out <- jcglasso.fit(Z = object$Z, diagonal = object$diagonal, weights.B = wB, weights.Tht = wTht, 
                        penalty = object$penalty, nrho = nrho, rho.min.ratio = object$rho.min.ratio, 
                        rho = object$rho, nlambda = nlambda, lambda.min.ratio = object$lambda.min.ratio, 
                        lambda = object$lambda, nu = object$nu, 
                        alpha1 = object$alpha1, alpha2 = object$alpha2, alpha3 = object$alpha3, 
                        maxit.em = object$maxit.em, 
                        thr.em = object$thr.em, maxit.bcd = object$maxit.bcd, thr.bcd = object$thr.bcd, 
                        trace = trace, offset = object$offset, covar2corr = object$covar2corr, 
                        truncate = object$truncate)
  }
  
  InfoStructure <- list(Adj = out$Adj, ncomp = out$ncomp, Ck = out$Ck, pk = out$pk,
                        id_X = out$id_X, id_Y = out$id_Y)
  for(i in seq_len(num_class)) {
    wTht[, , i][wTht[, , i] == .Machine$double.xmax] <- +Inf
    if(q > 0L) wB[, , i][wB[, , i] == .Machine$double.xmax] <- +Inf
  }
  out.jcggm <- list(call = this.call, Zipt = out$Zipt, B = out$B,
                    mu = out$mu, R = out$R, S = out$S, Sgm = out$Sgm, Tht = out$Tht,
                    Omega = out$Omega, dfB = out$dfB, dfTht = out$dfTht, dfOmg = out$dfOmg,
                    InfoStructure = InfoStructure, nit = out$nit, Z = object$Z, diagonal = object$diagonal, 
                    weights.B = wB, weights.Tht = wTht,
                    nlambda = nlambda, lambda.min.ratio = object$lambda.min.ratio, lambda = lambda.true, 
                    nrho = nrho, rho.min.ratio = object$rho.min.ratio, rho = rho.true, 
                    nu = object$nu, alpha1 = object$alpha1, alpha2 = object$alpha2, alpha3 = object$alpha3, 
                    connected = out$connected, penalty = object$penalty,
                    model = model, maxit.em = object$maxit.em, thr.em = object$thr.em,
                    maxit.bcd = object$maxit.bcd, thr.bcd = object$thr.bcd, conv = out$conv,
                    offset = object$offset, covar2corr = object$covar2corr, truncate = object$truncate,
                    subrout = out$subrout, trace = trace, nobs = object$nobs, nresp = object$nresp, npred = object$npred)
  
  class(out.jcggm) <- c("jcggm", "jcglasso")
  out.jcggm
}

jcggm.fit <- function (object, tp.min = 1.0E-6, trace = FALSE, ...) {
  
  Z <- object$Z
  
  diagonal <- object$diagonal
  storage.mode(diagonal) <- "integer"
  penalty <- object$penalty
  maxit.em <- object$maxit.em
  storage.mode(maxit.em) <- "integer"
  thr.em <- object$thr.em
  storage.mode(thr.em) <- "double"
  maxit.bcd <- object$maxit.bcd
  storage.mode(maxit.bcd) <- "integer"
  thr.bcd <- object$thr.bcd
  storage.mode(thr.bcd) <- "double"
  covar2corr <- object$covar2corr
  truncate <- object$truncate
  
  alpha1 <- object$alpha1
  storage.mode(alpha1) <- "double"
  alpha2 <- object$alpha2
  storage.mode(alpha2) <- "double"
  alpha3 <- object$alpha3
  storage.mode(alpha3) <- "double"
  
  storage.mode(tp.min) <- "double"
  storage.mode(trace) <- "integer"
  
  nrho <- 1L
  nlambda <- 1L
  
  K <- length(Z)
  storage.mode(K) <- "integer"
  
  n <- nobs(Z)
  storage.mode(n) <- "integer"
  weights <- proportions(n)
  storage.mode(weights) <- "double"
  
  p <- object$nresp
  storage.mode(p) <- "integer"
  
  q <- object$npred
  storage.mode(q) <- "integer"
  
  k_lab <- paste0("class", seq_len(K))
  rho_lab <- paste0("rho_", seq_len(nrho))
  lambda_lab <- paste0("lambda_", seq_len(nlambda))
  
  ynames <- colNames2(Z)[[1]]$Y
  xnames <- colNames2(Z)[[1]]$X
  znames <- c(ynames, xnames)
  
  X.null <- q == 0L
  
  seq_k <- seq_len(K)
  id_Y <- object$InfoStructure$id_Y
  id_X <- object$InfoStructure$id_X
  dim_Z <- p + q
  id_Z <- seq_len(dim_Z)
  storage.mode(dim_Z) <- "integer"
  
  weights.Tht <- object$weights.Tht
  if(!X.null) {
    weights.B <- object$weights.B
    dfOmg <- object$dfOmg
  } else weights.B <- object$weights.B
  
  ##### computing output and working objects #####
  Zmat <- array(0.0, dim = c(max(n), dim_Z, K),
                dimnames = list(NULL, c(ynames, xnames), k_lab))
  zm <- object$InfoStructure$zm
  zv <- object$InfoStructure$zv
  loz <- object$InfoStructure$lo
  upz <- object$InfoStructure$up
  Id <- object$InfoStructure$Id
  InfoP <- object$InfoStructure$InfoP
  nP <- object$InfoStructure$np
  Zipt <- object$Zipt
  B <- object$B
  mu <- object$mu
  R <- object$R
  S <- object$S
  Sgm <- object$Sgm
  Tht <- object$Tht
  Thtxx <- object$Omega
  B_n <- array(drop(B), dim = dim(B)[1:3], 
               dimnames = dimnames(B)[1:3])
  mu_n <- array(drop(mu), dim = dim(mu)[1:3], 
                dimnames = dimnames(mu)[1:3])
  R_n <- array(drop(R), dim = dim(R)[1:3], 
               dimnames = dimnames(R)[1:3])
  S_n <- array(drop(S), dim = dim(S)[1:3], 
               dimnames = dimnames(S)[1:3])
  Sgm_n <- array(drop(Sgm), dim = dim(Sgm)[1:3], 
                 dimnames = dimnames(Sgm)[1:3])
  Tht_n <- array(drop(Tht), dim = dim(Tht)[1:3], 
                 dimnames = dimnames(Tht)[1:3])
  Adj <- object$InfoStructure$Adj
  dfB <- object$dfB
  dfTht <- object$dfTht
  ncomp <- object$InfoStructure$ncomp
  Ck <- object$InfoStructure$Ck
  pk <- object$InfoStructure$Ck
  
  T1o <- zm * 0.0
  T2o <- S_n * 0.0
  T1 <- T1o[, 1]
  T2 <- T2o[, , 1]
  
  offset <- object$offset
  m_offset <- double(K)
  
  for(k in seq_k){
    Zmat[seq_len(n[k]), , k] <- cbind(getMatrix2(Z, name = "Y", ordered = TRUE)[[k]], 
                                      getMatrix2(Z, name = "X", ordered = TRUE)[[k]])
    Zmat[, , k][is.na(Zmat[, , k])] <- 0.0
    
    row.order <- object$Z[[k]]$Info$order
    
    Zipt[seq_len(n[k]), , k, , ] <- Zipt[seq_len(n[k]), , k, 1L, 1L][row.order, ] 
    mu[seq_len(n[k]), , k, , ] <- mu[seq_len(n[k]), , k, 1L, 1L][row.order, ] 
    mu_n[seq_len(n[k]), , k] <- mu_n[seq_len(n[k]), , k][row.order, ]
    R[seq_len(n[k]), , k, , ] <- R[seq_len(n[k]), , k, 1L, 1L][row.order, ] 
    R_n[seq_len(n[k]), , k] <- R_n[seq_len(n[k]), , k][row.order, ]
    
    T1o[, k] <- sapply(id_Z, function(j) { 
      id <- Id[[k]][, j] == 0L
      sum(Zmat[seq_len(n[k]), j, k][id]) 
    })
    
    for(i in id_Z) {
      for(j in i:dim_Z) {
        if(any(id <- Id[[k]][, i] == 0L & Id[[k]][, j] == 0L)) {
          T2o[i, j, k] <- sum(Zmat[seq_len(n[k]), i, k][id] * Zmat[seq_len(n[k]), j, k][id])
          T2o[j, i, k] <- T2o[i, j, k]
        }
      }
    }
    
    m_offset[k] <- mean(offset[seq_len(n[k]), k])
  }
  
  Zipt_lo <- zm - 3.0 * sqrt(zv)
  Zipt_up <- zm + 3.0 * sqrt(zv)
  Zipt_n <- array(drop(Zipt), dim = dim(Zipt)[1:3], 
                  dimnames = dimnames(Zipt)[1:3])
  
  nit.tot <- object$nit * 0L
  storage.mode(nit.tot) <- "integer"
  nit <- integer(2L)
  cnnctd <- object$connected[, 1L, 1L]
  conv <- integer(1)
  subrout <- integer(1)
  
  X <- array(0.0, c(max(n), q, K))
  # Sxx <- array(0.0, c(q, q, K))
  # Sxy <- array(0.0, c(q, p, K))
  
  ##### starting EM algorithm #####
  if(trace == 2) {
    if(!X.null){
      cat("\n*************************************************************************\n",
          "Fitting jcglasso model number ", 
          formatC(nrho, digits = 0, width = 5, format = "d"),
          "\n\t\t      rho = ",
          formatC(object$rho, digits = 6, width = 9, format = "f"),
          "\n\t\t   alpha1 = ",
          formatC(alpha1, digits = 3, width = 6, format = "f"),
          "\n\t\t   lambda = ",
          formatC(object$lambda, digits = 6, width = 9, format = "f"),
          "\n\t\t   alpha2 = ",
          formatC(alpha2, digits = 3, width = 6, format = "f"),
          "\n\t\t       nu = ", 
          formatC(object$nu, digits = 6, width = 9, format = "f"),
          "\n\t\t    alpha3 = ",
          formatC(alpha3, digits = 3, width = 6, format = "f"),
          "\n")
    } 
    else {
      cat("\n*************************************************************************\n",
          "Fitting jcglasso model number ", 
          formatC(nrho, digits = 0, width = 5, format = "d"),
          "\n\t\t      rho = ",
          formatC(object$rho, digits = 6, width = 9, format = "f"),
          "\n\t\t   alpha1 = ",
          formatC(alpha1, digits = 3, width = 6, format = "f"),
          "\n")
    }
  }
  
  dOmg <- 0.0
  for(ii in 1:maxit.em) {
    B_o <- B_n
    Tht_o <- Tht_n
    if(!X.null) Thtxx_o <- array(Thtxx[, , , 1L, 1L], dim = dim(Thtxx)[1:3], 
                                 dimnames = dimnames(Thtxx)[1:3])
    
    ##### computing E step #####
    for(k in seq_k) {
      temp <- .Fortran(cglasso:::C_e_step_v2, n = n[k], p = dim_Z, Y = Zmat[seq_len(n[k]), , k], 
                       lo = loz[, k], up = upz[, k], nP = nP[k], 
                       InfoP = InfoP[[k]], T1o = T1o[, k], T2o = T2o[, , k], 
                       mu_n = mu_n[seq_len(n[k]), , k], Sgm = Sgm_n[, , k], Tht = Tht_n[, , k], 
                       Yipt_lo = Zipt_lo[, k], Yipt_up = Zipt_up[, k], Yipt_n = Zipt_n[seq_len(n[k]), , k], 
                       T1 = T1, T2 = T2, conv = conv)
      if(temp$conv != 0) { 
        subrout <- 1
        stop("error in E step!") 
      }
      zm[, k] <- temp$T1 / n[k]
      Zipt_n[seq_len(n[k]), , k] <- temp$Yipt_n
      S_n[, , k] <- temp$T2
      B_n[1L, , k] <- zm[id_Y, k] - m_offset[k]
      mu_n[seq_len(n[k]), id_Y, k] <- outer(offset[seq_len(n[k]), k], B_n[1L, , k], FUN = "+")
      R_n[seq_len(n[k]), , k] <- Zipt_n[seq_len(n[k]), id_Y, k] - mu_n[seq_len(n[k]), id_Y, k]
      # S_n[[k]][id_Y, id_Y] <- crossprod(R_n[[k]]) / n[k]
      
      if(!X.null) {
        mu_n[seq_len(n[k]), id_X, k] <- rep(zm[id_X, k], each = n[k])
        X[seq_len(n[k]), , k] <- sweep(matrix(Zipt_n[seq_len(n[k]), id_X, k], n[k], q), 2, zm[id_X, k], "-")
      }
    }
    if(trace == 2) cat("\n\tE-step completed!")
    
    ##### fitting multilasso model #####
    if(!X.null) {
      if(trace == 2) cat("\n\tM-step:\n\t       fitting multivariate sparse group lasso model with ",
                         "\n\t       lambda = ", formatC(object$lambda, digits = 6, width = 9, format = "f"), 
                         "and alpha2 = ", formatC(alpha2, digits = 3, width = 4, format = "f"), 
                         "\n\n")
      
      # if(trace == 2) {
      #   cat("\t\tADMM step\t||B_new(k) - B_old(k)||_1/||B_old(k)||_1\n")
      #   # cat("\t\t   ", formatC(tmpB$nit[2], digits = 0, width = 5, format = "d"),
      #   #     "\t\t\t", formatC(tempmul$diff, digits = 8, width = 10, format = "f"),"\n")
      # }
      
      # for(k in seq_k) {
      # X[seq_len(n[k]), , k] <- sweep(Zipt_n[seq_len(n[k]), id_X, k], 2, zm[id_X, k], "-")
      # Sxx[, , k] <- crossprod(X[seq_len(n[k]), , k]) / n[k]
      # Sxy[, , k] <- crossprod(X[seq_len(n[k]), , k], R_n[seq_len(n[k]), , k]) / n[k]
      # }
      
      tmpB <- .Fortran(cglasso:::C_apg, p = p*K, q = q*K, n = max(n)*K, k = K,
                       nk = as.double(rep(n, each = q)), fk = as.double(rep(weights, each = q)),
                       A = do.call(blockdiag, array2list(X)), b = do.call(blockdiag, array2list(R_n)),
                       Tht = do.call(blockdiag, array2list(Tht_n[id_Y, id_Y, , drop = FALSE])),
                       weights = do.call(blockdiag, array2list(weights.B)),
                       lambda = tp.min, alpha = alpha2, xm = zm[id_X, ], ym = zm[id_Y, ],
                       x = do.call(blockdiag, array2list(B_n[-1L, , , drop = FALSE])),
                       beta = B_n,
                       maxit = as.integer(maxit.bcd), thr = thr.bcd, trace = trace,
                       df = matrix(0L, p + 1L, K), nit = integer(1), conv = integer(1))
      dfB[, , 1L, 1L] <- tmpB$df
      nit[2L] <- nit[2L] + tmpB$nit[1L]
      B_n <- tmpB$beta
      
      # tmpB <- .Fortran(cglasso:::C_admm_b_joint, p = p, q = q, qp = as.integer(p*q), K = K, 
      #                  fk = weights, xtx = Sxx, xtr = Sxy, B = B_n[-1, , , drop = FALSE], 
      #                  Tht = Tht_n[id_Y, id_Y, , drop = FALSE], wB = weights.B, 
      #                  lmb1 = tp.min, alpha = alpha, maxit = maxit.bcd, thr = thr.bcd, rho = 1.0, trace = trace, 
      #                  conv = integer(1), subrout = integer(1), nit = integer(2), df = matrix(0L, p + 1, K))
      # if(tmpB$conv != 0) {
      #   conv <- tmpB$conv
      #   subrout <- 2
      #   stop("error in multilasso step!")
      # }
      # dfB[, , 1L, 1L] <- tmpB$df
      # nit[2] <- nit[2] + tmpB$nit[2]
      for(k in seq_k){
        # B_n[1L, , k] <- zm[id_Y, k] - zm[id_X, k] %*% tmpB$B[, , k]
        # B_n[-1L, , k] <- tmpB$B[, , k]
        mu_n[seq_len(n[k]), id_Y, k] <- cbind(1, Zipt_n[seq_len(n[k]), id_X, k]) %*% B_n[, , k]
        R_n[seq_len(n[k]), , k] <- Zipt_n[seq_len(n[k]), id_Y, k] - mu_n[seq_len(n[k]), id_Y, k] - offset[seq_len(n[k]), k]
      }
    }
    for(k in seq_k){
      YM <- crossprod(Zipt_n[seq_len(n[k]), , k], mu_n[seq_len(n[k]), , k])
      S_n[, , k] <- (S_n[, , k] + crossprod(mu_n[seq_len(n[k]), , k]) - YM - t(YM)) / n[k]
      if(covar2corr) S_n[[k]] <- cov2cor(S_n[[k]])
    }
    
    ##### computing theta and S values for connected components #####
    if(trace == 2) {
      if(X.null) {
        cat("\n")
        cat("\tM-step:\n\t       fitting joint glasso model with",
            "rho = ", formatC(object$rho, digits = 6, width = 9, format = "f"),
            "and alpha1 = ", formatC(alpha1, digits = 3, width = 4, format = "f"),
            "\n\n")
      } 
      else {
        cat("\tM-step:\n\t       fitting joint glasso model with",
            "\n\t      rho = ", formatC(object$rho, digits = 6, width = 9, format = "f"),
            "and alpha1 = ", formatC(alpha1, digits = 3, width = 4, format = "f"),
            "\n\t       nu = ", formatC(object$nu, digits = 6, width = 9, format = "f"), 
            "and alpha3 = ", formatC(alpha3, digits = 3, width = 4, format = "f"),
            "\n\n")
      }
    }
    tmpTht <- .Fortran(cglasso:::C_admm_tht_sub, p = p, N = K, fk = weights, S = S_n[id_Y, id_Y, , drop = FALSE], 
                       wTht = weights.Tht[id_Y, id_Y, , drop = FALSE], pendiag = diagonal, rho = tp.min, alpha = alpha1, 
                       maxit = maxit.bcd, thr = thr.bcd, Tht = Tht_n[id_Y, id_Y, , drop = FALSE], k = integer(1),
                       Ck = integer(p), pk = integer(p), nit = integer(1), df = integer(K), 
                       conv = integer(1), trace = trace)
    nit[2] <- nit[2] + tmpTht$nit
    Tht_n[id_Y, id_Y, ] <- tmpTht$Tht
    
    if(!X.null) { 
      tmpOmg <- .Fortran(cglasso:::C_admm_tht_sub, p = q, N = K, fk = weights, S = S_n[id_X, id_X, , drop = FALSE], 
                         wTht = weights.Tht[id_X, id_X, , drop = FALSE], pendiag = diagonal, rho = tp.min, alpha = alpha3, 
                         maxit = maxit.bcd, thr = thr.bcd, Tht = Tht_n[id_X, id_X, , drop = FALSE], k = integer(1),
                         Ck = integer(q), pk = integer(q), nit = integer(1), df = integer(K), 
                         conv = integer(1), trace = trace)
      nit[2] <- nit[2] + tmpOmg$nit
      Thtxx[, , , 1L, 1L] <- tmpOmg$Tht
    }
    
    nit[1] <- ii
    
    for(k in seq_k) {
      Adj[id_Y, id_Y, k, 1L, 1L] <- 1*(Tht_n[id_Y, id_Y, k] != 0)
      if(!X.null) { 
        Adj[id_X, id_X, k, 1L, 1L] <- 1*(Thtxx[, , k, 1L, 1L] != 0)
        Tht_n[id_X, id_X, k] <- Thtxx[, , k, 1L, 1L] + B_n[-1L, , k] %*% Tht_n[id_Y, id_Y, k] %*% t(matrix(B_n[-1L, , k], q, p))
        Tht_n[id_X, id_Y, k] <- -(B_n[-1L, , k] %*% Tht_n[id_Y, id_Y, k])
        Tht_n[id_Y, id_X, k] <- t(matrix(Tht_n[id_X, id_Y, k], q, p))
      }
      Sgm_n[, , k] <- solve(Tht_n[, , k])
    } 
    
    dB <- if(!X.null)
      sum(weights * sapply(seq_len(K), function(k) norm(B_n[, , k] - B_o[, , k], type = "F") / (p + dfB[p + 1, k, 1L, 1L])))
    else sum(weights * sapply(seq_len(K), function(k) norm(B_n[, , k] - B_o[, , k], type = "2") / (p + dfB[p + 1, k, 1L, 1L])))
    dTht <- sum(weights * sapply(seq_len(K), function(k) norm(Tht_n[id_Y, id_Y, k] - Tht_o[id_Y, id_Y, k], type = "F") / (p + tmpTht$df[k])))
    if(!X.null) dOmg <- sum(weights * sapply(seq_len(K), function(k) norm(matrix(Thtxx[, , k, 1L, 1L], q, q) - matrix(Thtxx_o[, , k], q, q), "F") / (q + tmpOmg$df[k])))
    
    # dB <- if(!X.null) sqrt(sum(B_n - B_o)^2) / sqrt(sum(B_n^2))
    # dTht <- sqrt(sum(Tht_n[id_Y, id_Y, ] - Tht_o[id_Y, id_Y, ])^2) / sqrt(sum(Tht_n[id_Y, id_Y, ]^2))
    # if(!X.null) dOmg <- sqrt(sum(Thtxx[, , , 1L, 1L] - Thtxx_o)^2) / sqrt(sum(Thtxx[, , , 1L, 1L]^2))
    
    if(trace == 2) {
      if(!X.null) {
        cat("\n\tChecking convergence criterion (threshold = ", 
            formatC(thr.em, digits = 6, width = 8, format = "f"),
            ")\n\t||B_old - B_new|_F / dfB  = ",
            formatC(dB, digits = 7, width = 9, format = "f"),
            "\n\t||Tht_old - Tht_new||_F / dfTht = ",
            formatC(dTht, digits = 7, width = 9, format = "f"), 
            "\n\t||Omg_old - Omg_new||_F / dfOmg = ", 
            formatC(dOmg, digits = 7, width = 9, format = "f"), 
            "\n")
        if(max(dB, dTht, dOmg) <= thr.em) cat("\tConvergence criterion is met!\n\n")
      }
      else {
        cat("\n\tChecking convergence criterion (threshold = ", 
            formatC(thr.em, digits = 6, width = 8, format = "f"),
            ")\n\t||B_old - B_new|_F / dfB  = ",
            formatC(dB, digits = 7, width = 9, format = "f"),
            "\n\t||Tht_old - Tht_new||_F / dfTht = ",
            formatC(dTht, digits = 7, width = 9, format = "f"), 
            "\n")
        if(max(dB, dTht) <= thr.em) cat("\tConvergence criterion is met!\n\n")
      }
    } 
    if(max(dB, dTht, dOmg) <= thr.em) break
  }
  if(ii >= maxit.em) {
    conv <- 1
    subrout <- 3
    stop("error in M step (unconnected components)!")
  }
  if(trace == 1) {
    if(!X.null){
      cat("\njcglasso model number ", 
          formatC(nrho, digits = 0, width = 5, format = "d"),
          "\n\t        rho = ",
          formatC(object$rho, digits = 6, width = 9, format = "f"),
          "\n\t     alpha1 = ",
          formatC(alpha1, digits = 3, width = 6, format = "f"),
          "\n\t     lambda = ",
          formatC(object$lambda, digits = 6, width = 9, format = "f"),
          "\n\t     alpha2 = ",
          formatC(alpha2, digits = 3, width = 6, format = "f"),
          "\n\t         nu = ",
          formatC(object$nu, digits = 6, width = 9, format = "f"),
          "\n\t     alpha3 = ",
          formatC(alpha3, digits = 3, width = 6, format = "f"),
          "\n\t   EM-steps = ",
          formatC(nit[1], digits = 0, width = 5, format = "d"),
          "\n\t      steps = ",
          formatC(nit[2], digits = 0, width = 5, format = "d"),
          "\n\n"
      )
    } 
    else {
      cat("\njcglasso model number ", 
          formatC(nrho, digits = 0, width = 5, format = "d"),
          "\n\t        rho = ",
          formatC(object$rho, digits = 6, width = 9, format = "f"),
          "\n\t     alpha1 = ",
          formatC(alpha1, digits = 3, width = 6, format = "f"),
          "\n\t   EM-steps = ",
          formatC(nit[1], digits = 0, width = 5, format = "d"),
          "\n\t      steps = ",
          formatC(nit[2], digits = 0, width = 5, format = "d"),
          "\n\n"
      )
    }
  }
  
  ##### buiding output components #####
  for(k in seq_len(K)) {
    Zipt[, , k, 1L, 1L] <- Zipt_n[, , k]
    B[, , k, 1L, 1L] <- B_n[, , k]
    mu[, , k, 1L, 1L] <- mu_n[, , k]
    R[, , k, 1L, 1L] <- R_n[, , k]
    S[, , k, 1L, 1L] <- S_n[, , k]
    Sgm[, , k, 1L, 1L] <- Sgm_n[, , k]
    Tht[, , k, 1L, 1L] <- Tht_n[, , k]
    
    diag(Adj[, , k, 1L, 1L]) <- 0
    
    dfTht[k, 1L, 1L] <- tmpTht$df[k]
    if(!X.null) {
      Adj[id_X, id_Y, k, 1L, 1L] <- 1*(B[-1, , k, 1L, 1L] != 0)
      Adj[id_Y, id_X, k, 1L, 1L] <- t(Adj[id_X, id_Y, k, 1L, 1L])
      dfOmg[k, 1L, 1L] <- tmpOmg$df[k]
    }
    
    row.order <- order(Z[[k]]$Info$order)
    Zipt[seq_len(n[k]), , k, 1L, 1L] <- Zipt[seq_len(n[k]), , k, 1L, 1L][row.order, , drop = FALSE]
    mu[seq_len(n[k]), , k, 1L, 1L] <- mu[seq_len(n[k]), , k, 1L, 1L][row.order, , drop = FALSE]
    R[seq_len(n[k]), , k, 1L, 1L] <- R[seq_len(n[k]), , k, 1L, 1L][row.order, , drop = FALSE]
  }
  nit.tot[, 1L, 1L] <- nit
  
  names(Z) <- k_lab
  
  out <- list(n = n, p = p, q = q, id_X = id_X, id_Y = id_Y, Z = Zmat, 
              Id = Id, nP = nP, InfoP = InfoP, ncomp = ncomp, Ck = Ck, pk = pk, 
              lo = loz, up = upz, zm = zm, zv = zv, wTht = weights.Tht, 
              nrho = nrho, rhoratio = object$rho.min.ratio, rho = tp.min,
              nlambda = nlambda, lambdaratio = object$lambda.min.ratio, lambda = tp.min, 
              nu = tp.min, alpha1 = alpha1, alpha2 = alpha2, alpha3 = alpha3, pendiag = diagonal, connected = cnnctd, 
              maxit_em = maxit.em, thr_em = thr.em, maxit_bcd = maxit.bcd, thr_bcd = thr.bcd, 
              Zipt = Zipt, B = B, mu = mu, R = R, S = S, Sgm = Sgm, Tht = Tht, 
              Adj = Adj, Omega = Thtxx, dfB = dfB, dfTht = dfTht, dfOmg = NULL,
              nit = nit.tot, conv = conv, offset = offset,
              subrout = subrout, trace = trace)
  if(!X.null) out$dfOmg <- dfOmg
  
  out$conv <- switch(as.character(out$conv),
                     "-1" = "memory allocation error",
                     "0" = "Ok",
                     "1" = "maximum number of iterations has been exceeded",
                     "2" = "error in E-step",
                     "3" = "matrix inversion failed")
  return(out)
}
