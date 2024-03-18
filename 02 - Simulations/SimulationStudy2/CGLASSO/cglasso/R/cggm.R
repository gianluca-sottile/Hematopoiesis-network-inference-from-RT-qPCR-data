cggm <- function(object, GoF = AIC, lambda.id, rho.id, tp.min = 1.0E-6, ntp = 100L, maxit.em = 1.0E+4, thr.em = 1.0E-3, maxit.bcd = 1.0E+5,
                 thr.bcd = 1.0E-4, trace = 0L, algorithm = c("glasso", "admm"), ...){
    this.call <- match.call()
    algorithm <- match.arg(algorithm)
    # testing 'object'
    if (!inherits(object, "cglasso")) stop(sQuote("object"), " is not an object of class ", sQuote("cglasso"))
    Z <- object$Z
    n <- nobs(Z)
    p <- nresp(Z)
    q <- npred(Z)
    if (missing(lambda.id) & missing(rho.id)) {
        if (!is.element(class(GoF), c("function", "GoF")))
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
        object <- select.cglasso(object, GoF = GoF)
        nlambda <- object$nlambda
        lambda.id <- 1L
        nrho <- object$nrho
        rho.id <- 1L
    } else {
        # testing 'lambda.id'
        nlambda <- object$nlambda
        if (missing(lambda.id)) {
            if (q == 0 | nlambda == 1L) lambda.id <- 1L
            else stop(sQuote("lambda.id"), " is missing")
        } else {
            if (!is.vector(lambda.id)) stop(sQuote("lambda.id"), " is not a vector")
            if (length(lambda.id) != 1L) stop(sQuote("lambda.id"), " is not a vector of length ", sQuote("1"))
            if (any(abs(as.integer(lambda.id) - lambda.id) > 0)) stop(sQuote("lambda.id"), " is not an object of type ", dQuote("integer"))
            if (lambda.id <= 0) stop(sQuote("lambda.id"), " is not a positive integer")
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
    }
    # testing 'tp.min'
    if (!is.vector(tp.min)) stop(sQuote("tp.min"), " is not a vector")
    if (length(tp.min) != 1L) stop(sQuote("tp.min"), " is not a vector of length ", sQuote("1"))
    if (tp.min < 0L) stop(sQuote("tp.min"), " is not a positive value")
    # testing 'ntp'
    if (!is.vector(ntp)) stop(sQuote("ntp"), " is not a vector")
    if (length(ntp) != 1L) stop(sQuote("ntp"), " is not a vector of length ", sQuote("1"))
    if (any(abs(as.integer(ntp) - ntp) > 0)) stop(sQuote("ntp"), " is not an object of type ", dQuote("integer"))
    if (ntp <= 0L) stop(sQuote("ntp"), " is not a positive integer")
    # Testing 'maxit.em'
    if (!is.vector(maxit.em)) stop(sQuote("maxit.em"), " is not a vector")
    if (length(maxit.em) != 1) stop(sQuote("maxit.em"), " is not an object of length ", sQuote(1))
    if (abs(as.integer(maxit.em) - maxit.em) > 0) stop(sQuote("maxit.em"), " is not an object of type ", dQuote("integer"))
    if (maxit.em <= 0) stop(sQuote("maxit.em"), " is not a positive integer")
    # Testing 'thr.em'
    if (!is.vector(thr.em)) stop(sQuote("thr.em"), " is not a vector")
    if (length(thr.em) != 1) stop(sQuote("thr.em"), " is not an object of length ", sQuote(1))
    if (thr.em <= 0 ) stop(sQuote("thr.em"), " is not a positive value")
    # Testing 'maxit.bcd'
    if (!is.vector(maxit.bcd)) stop(sQuote("maxit.bcd"), " is not a vector")
    if (length(maxit.bcd) != 1) stop(sQuote("maxit.bcd"), " is not an object of length ", sQuote(1))
    if (abs(as.integer(maxit.bcd) - maxit.bcd) > 0) stop(sQuote("maxit.bcd"), " is not an object of type ", dQuote("integer"))
    if (maxit.bcd <= 0) stop(sQuote("maxit.bcd"), " is not a positive integer")
    # Testing 'thr.bcd'
    if (!is.vector(thr.bcd)) stop(sQuote("thr.bcd"), " is not a vector")
    if (length(thr.bcd) != 1) stop(sQuote("thr.bcd"), " is not an object of length ", sQuote(1))
    if (thr.bcd <= 0 ) stop(sQuote("thr.bcd"), " is not a positive value")
    # Testing 'trace'
    if (!is.vector(trace)) stop(sQuote("trace"), " is not a vector")
    if (length(trace) != 1) stop(sQuote("trace"), " is not an object of length ", sQuote(1))
    if (is.logical(trace)) stop(sQuote("trace"), " is not an object of type ", dQuote("integer"))
    if (abs(as.integer(trace) - trace) > 0) stop(sQuote("trace"), " is not an object of type ", dQuote("integer"))
    if (!is.element(trace, c(0L, 1L, 2L))) stop("not allowed value in ", sQuote("trace"), ". Please, choice ", sQuote("0"), ", ", sQuote("1"), " or ", sQuote("2"))
    pendiag <- object$diagonal
    lambda.max <- object$lambda[lambda.id]
    lambda.min <- object$lambda[nlambda]
    rho.max <- object$rho[rho.id]
    rho.min <- object$rho[nrho]
    B.ini <- object$B[, , lambda.id, rho.id]
    if (is.vector(B.ini)) B.ini <- t(B.ini)
    if (q > 0L) {
        mask.B <- B.ini[-1L, ]
        mask.B[abs(mask.B) > 0] <- 1
    } else mask.B <- NULL
    Tht.ini <- object$Tht[, , lambda.id, rho.id]
    mask.Tht <- Tht.ini
    mask.Tht[abs(mask.Tht) > 0] <- 1
    row.order <- Z$Info$order
    Yipt.ini <- object$Yipt[row.order, , lambda.id, rho.id]
    mu.ini <- object$mu[row.order, , lambda.id, rho.id]
    R.ini <- object$R[row.order, , lambda.id, rho.id]
    S.ini <- object$S[, , lambda.id, rho.id]
    Sgm.ini <- object$Sgm[, , lambda.id, rho.id]
    tp.min <- max(min(tp.min, lambda.min, rho.min), 1e-13)
    out <- cggm.fit(Z = Z, pendiag = pendiag, lambda.max = lambda.max, rho.max = rho.max, tp.min = tp.min,
                    ntp = ntp, mask.B = mask.B, mask.Tht = mask.Tht, Yipt.ini = Yipt.ini, B.ini = B.ini,
                    mu.ini = mu.ini, R.ini = R.ini, S.ini = S.ini, Sgm.ini = Sgm.ini, Tht.ini = Tht.ini,
                    maxit.em = maxit.em, thr.em = thr.em, maxit.bcd = maxit.bcd, thr.bcd = thr.bcd,
                    trace = trace, algorithm = algorithm)
    
    InfoStructure <- list(Adj_yy = object$InfoStructure$Adj_yy[, , lambda.id, rho.id, drop = FALSE], ncomp = object$InfoStructure$ncomp[lambda.id, rho.id, drop = FALSE],
                          Ck = object$InfoStructure$Ck[, lambda.id, rho.id, drop = FALSE], pk = object$InfoStructure$pk[, lambda.id, rho.id, drop = FALSE],
                          Adj_xy = object$InfoStructure$Adj_xy[, , lambda.id, rho.id, drop = FALSE])
                          
    out.cggm <- list(call = this.call, Yipt = out$Yipt, B = out$B, mu = out$mu, R = out$R, S = out$S, Sgm = out$Sgm, Tht = out$Tht,
                     dfB = object$dfB[, lambda.id, rho.id, drop = FALSE], dfTht = object$dfTht[lambda.id, rho.id, drop = FALSE],
                     InfoStructure = InfoStructure, nit = out$nit, Z = Z, nlambda = 1L, lambda = lambda.max, nrho = 1L,
                     rho = rho.max, maxit.em = maxit.em, thr.em = thr.em, maxit.bcd = maxit.bcd, thr.bcd = thr.bcd,
                     conv = out$conv, subrout = out$subrout, trace = trace, nobs = object$n, nresp = object$nresp, npred = object$npred)
    
    class(out.cggm) <- c("cggm", "cglasso")
    out.cggm
}

cggm.fit <- function(Z, pendiag, lambda.max, rho.max, tp.min, ntp, mask.B, mask.Tht, Yipt.ini, B.ini,
                 mu.ini, R.ini, S.ini, Sgm.ini, Tht.ini, maxit.em, thr.em, maxit.bcd, thr.bcd, trace, algorithm) {
    n <- nobs(Z)
    p <- nresp(Z)
    q <- npred(Z)
    Y <- getMatrix(Z, "Y", ordered = TRUE)
    Y[is.na(Y)] <- 0
    X <- getMatrix(Z, "X", ordered = TRUE)
    if (!is.null(X)) X <- as.matrix(X)
    ynames <- colnames(Y)
    xnames <- colnames(X)
    yrnames <- rownames(Y)
    ym <- ColMeans(Z)$Y
    yv <- ColVars(Z)$Y
    names(ym) <- NULL
    names(yv) <- NULL
    lo <- Z$Info$lo
    up <- Z$Info$up
    Id <- event(Z, ordered = TRUE)
    InfoP <- Z$Info$Pattern
    nP <- dim(InfoP)[1L]
    tp_lab <- paste0("tp", seq_len(ntp))
    lambda <- seq(from = lambda.max, to = tp.min, length = ntp)
    rho <- seq(from = rho.max, to = tp.min, length = ntp)
    B <- array(0, dim = c(q + 1L, p, 1L, 1L), dimnames = list(coef = c("Int.", xnames), response = ynames, lambda = "lmb", rho = "rho"))
    B[, , 1L, 1L] <- B.ini
    Yipt <- array(0, dim = c(n, p, 1L, 1L), dimnames = list(yrnames, response = ynames, lambda = "lmb", rho = "rho"))
    Yipt[, , 1L, 1L] <- Yipt.ini
    mu <- array(0, dim = c(n, p, 1L, 1L), dimnames = list(yrnames, response = ynames, lambda = "lmb", rho = "rho"))
    mu[, , 1L, 1L] <- mu.ini
    R <- array(0, dim = c(n, p, 1L, 1L), dimnames = list(yrnames, response = ynames, lambda = "lmb", rho = "rho"))
    R[, , 1L, 1L] <- R.ini
    S <- array(0, dim = c(p, p, 1L, 1L), dimnames = list(response = ynames, response = ynames, lambda = "lmb", rho = "rho"))
    S[, , 1L, 1L] <- S.ini
    Sgm <- array(0, dim = c(p, p, 1L, 1L), dimnames = list(response = ynames, response = ynames, lambda = "lmb", rho = "rho"))
    Sgm[, , 1L, 1L] <- Sgm.ini
    Tht <- array(0, dim = c(p, p, 1L, 1L), dimnames = list(response = ynames, response = ynames, lambda = "lmb", rho = "rho"))
    Tht[, , 1L, 1L] <- Tht.ini
    nit <- array(0L, dim = c(2L, ntp), dimnames = list(steps = c("EM", "nit"), tp = tp_lab))
    # setting storage mode
    storage.mode(n) <- "integer"
    storage.mode(p) <- "integer"
    storage.mode(Y) <- "double"
    storage.mode(Id) <- "integer"
    storage.mode(nP) <- "integer"
    storage.mode(InfoP) <- "integer"
    storage.mode(lo) <- "double"
    storage.mode(up) <- "double"
    storage.mode(ym) <- "double"
    storage.mode(yv) <- "double"
    storage.mode(pendiag) <- "integer"
    storage.mode(mask.Tht) <- "double"
    storage.mode(ntp) <- "integer"
    storage.mode(rho) <- "double"
    storage.mode(maxit.em) <- "integer"
    storage.mode(thr.em) <- "double"
    storage.mode(maxit.bcd) <- "integer"
    storage.mode(thr.bcd) <- "double"
    storage.mode(Yipt) <- "double"
    storage.mode(B) <- "double"
    storage.mode(mu) <- "double"
    storage.mode(R) <- "double"
    storage.mode(S) <- "double"
    storage.mode(Sgm) <- "double"
    storage.mode(Tht) <- "double"
    storage.mode(nit) <- "integer"
    conv <- integer(1)
    subrout <- integer(1)
    storage.mode(trace) <- "integer"
    if (q == 0L) {
        if(algorithm == "glasso"){
            out <- .Fortran(C_cggm_v1, n = n, p = p, Y = Y, Id = Id, nP = nP, InfoP = InfoP, lo = lo, up = up,
                            ym = ym, yv = yv, pendiag = pendiag, wTht = mask.Tht, ntp = ntp, rho = rho,
                            maxit_em = maxit.em, thr_em = thr.em, maxit_bcd = maxit.bcd, thr_bcd = thr.bcd,
                            Yipt = Yipt, B = B, mu = mu, R = R, S = S, Sgm = Sgm, Tht = Tht, nit = nit,
                            conv = conv, subrout = subrout, trace = trace)
        } else {
            out <- .Fortran(C_cggm_v3, n = n, p = p, Y = Y, Id = Id, nP = nP, InfoP = InfoP, lo = lo, up = up,
                            ym = ym, yv = yv, pendiag = pendiag, wTht = mask.Tht, ntp = ntp, rho = rho,
                            maxit_em = maxit.em, thr_em = thr.em, maxit_bcd = maxit.bcd, thr_bcd = thr.bcd,
                            Yipt = Yipt, B = B, mu = mu, R = R, S = S, Sgm = Sgm, Tht = Tht, nit = nit,
                            conv = conv, subrout = subrout, trace = trace)
        }
        out$subrout <- switch(as.character(out$subrout),
                              "0" = "",
                              "1" = "e_step",
                              "2" = "glassosub",
                              "3" = "cggm_v1")
    }
    else {
        storage.mode(q) <- "integer"
        storage.mode(X) <- "double"
        storage.mode(mask.B) <- "double"
        storage.mode(lambda) <- "double"
        if(algorithm == "glasso"){
            out <- .Fortran(C_cggm_v2, n = n, q = q, X = X, p = p, Y = Y, Id = Id, nP = nP, InfoP = InfoP,
                            lo = lo, up = up, ym = ym, yv = yv, wB = mask.B, pendiag = pendiag, wTht = mask.Tht ,
                            ntp = ntp, lambda = lambda ,rho = rho, maxit_em = maxit.em, thr_em = thr.em,
                            maxit_bcd = maxit.bcd, thr_bcd = thr.bcd, Yipt = Yipt, B = B, mu = mu, R = R, S = S,
                            Sgm = Sgm, Tht = Tht, nit = nit, conv = conv, subrout = subrout, trace = trace)
        } else {
            out <- .Fortran(C_cggm_v4, n = n, q = q, X = X, p = p, Y = Y, Id = Id, nP = nP, InfoP = InfoP,
                            lo = lo, up = up, ym = ym, yv = yv, wB = mask.B, pendiag = pendiag, wTht = mask.Tht ,
                            ntp = ntp, lambda = lambda ,rho = rho, maxit_em = maxit.em, thr_em = thr.em,
                            maxit_bcd = maxit.bcd, thr_bcd = thr.bcd, Yipt = Yipt, B = B, mu = mu, R = R, S = S,
                            Sgm = Sgm, Tht = Tht, nit = nit, conv = conv, subrout = subrout, trace = trace)
        }
        out$subrout <- switch(as.character(out$subrout),
                              "0" = "",
                              "1" = "e_step",
                              "2" = "multilasso",
                              "3" = "glassosub",
                              "4" = "cggm_v2")
    }
    row.order <- order(Z$Info$order)
    out$Yipt <- out$Yipt[row.order, , , , drop = FALSE]
    out$mu <- out$mu[row.order, , , , drop = FALSE]
    out$R <- out$R[row.order, , , , drop = FALSE]
    out$conv <- switch(as.character(out$conv),
                       "-1" = "memory allocation error",
                       "0" = "Ok",
                       "1" = "maximum number of iterations has been exceeded",
                       "2" = "error in E-step",
                       "3" = "matrix inversion failed")
    out
}
