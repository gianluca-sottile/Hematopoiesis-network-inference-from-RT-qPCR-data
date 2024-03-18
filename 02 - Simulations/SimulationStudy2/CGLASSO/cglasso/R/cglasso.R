cglasso <- function(formula, data, subset, contrasts = NULL, diagonal = FALSE, weights.B = NULL,
                    weights.Tht = NULL, nlambda, lambda.min.ratio, lambda, nrho, rho.min.ratio, rho,
                    maxit.em = 1.0E+4, thr.em = 1.0E-3, maxit.bcd = 1.0E+5, thr.bcd = 1.0E-4, trace = 0L,
                    algorithm = c("glasso", "admm")) {
    this.call <- match.call()
    algorithm <- match.arg(algorithm)
    zero <- 1.0E-6
    # testing 'formula'
    if (missing(formula)) formula <- . ~ . # stop(sQuote("formula"), " is missing")
    # testing 'data'
    if (!is.datacggm(data)) stop(sQuote("data"), " is not an object of class ", sQuote("datacggm"))
    ####################
    # formula2datacggm #
    ####################
    # testing LHS 'formula'
    fmlTerms <- function(fml, pos) paste(deparse(fml[[pos]], width.cutoff = 500L), collapse = " ")
    if (fmlTerms(formula, 2L) == ".")
#    if (deparse1(formula[[2L]]) == ".")
        formula <- formula(paste0(paste0("cbind(", paste(colNames(data)$Y, collapse = ", "), ")"), " ~ ", fmlTerms(formula, 3L)))
#        formula <- formula(paste0(paste0("cbind(", paste(colNames(data)$Y, collapse = ", "), ")"), " ~ ", deparse1(formula[[3L]])))
    if (as.character(formula[[2L]])[1L] != "cbind")
        stop("Please use ", sQuote("cbind"), " to specify LHS in ", sQuote("formula"), " object")
    Y.LHS <- unlist(lapply(formula[[2L]][-1L], deparse))
    # Y.LHS <- all.vars(formula[[2L]], unique = FALSE)
    if (any(table(Y.LHS) > 1L))
        stop("repeated response variables are not permitted in ", sQuote("formula"), " object")
    noVars <- !is.element(Y.LHS, colNames(data)$Y)
    if (any(noVars)) stop("Following variables are not stored as rensponse variables: ", Y.LHS[noVars])
    # testing RHS 'formula'
    if (is.null(getMatrix(data, name = "X"))) {
        if (fmlTerms(formula, 3L) == ".") formula <- update(formula, . ~ 1)
#        if (deparse1(formula[[3L]]) == ".") formula <- update(formula, . ~ 1)
        mt <- terms(formula)
        if (length(attr(mt, which = "term.labels")) != 0L)
            stop("Predictors are not stored in ", sQuote("data"))
        if (!as.logical(attr(mt, "intercept"))) {
            warning("Current version does not fit models without intercept term, thus it is added to the current formula.")
            formula <- update(formula, . ~ 1)
            mt <- terms(formula)
        }
        data.df <- data.frame(getMatrix(data, name = "Y")[, Y.LHS, drop = FALSE])
        nResp <- length(Y.LHS)
        nPred <- 0L
    } else {
        FML.RHS <- attr(terms(formula, data = getMatrix(data, name = "X")), "term.labels")
        if (length(FML.RHS) > 0L) {
            FML.RHS <- paste(FML.RHS, collapse = " + ")
            formula <- formula(paste0(fmlTerms(formula, 2L), " ~ ", FML.RHS))
#            formula <- formula(paste0(deparse1(formula[[2L]]), " ~ ", FML.RHS))
            X.RHS <- all.vars(formula[[3L]])
            noVars <- !is.element(X.RHS, colNames(data)$X)
            if (any(noVars)) stop("Following variables are not stored as predictors: ", X.RHS[noVars])
            data.df <- data.frame(getMatrix(data, name = "Y")[, Y.LHS, drop = FALSE], getMatrix(data, name = "X")[, X.RHS, drop = FALSE])
            nPred <- length(X.RHS)
        } else {
            data.df <- data.frame(getMatrix(data, name = "Y")[, Y.LHS, drop = FALSE])
            nPred <- 0L
        }
        nResp <- length(Y.LHS)
    }
    mt <- terms(formula)
    if (!as.logical(attr(mt, "intercept"))) {
        warning("Current version does not fit models without intercept term, thus it has been added to the current formula.")
        formula <- update(formula, . ~ . + 1)
        mt <- terms(formula)
    }
    # creating model.frame
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "subset"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf$na.action <- na.pass
    mf$drop.unused.levels <- TRUE
    mf$formula <- formula
    mf$data <- data.df
    mf[[1L]] <- quote(stats::model.frame)
    mf <- eval(mf, parent.frame())
    # updating 'datacggm' obejct
    Y <- model.response(mf, type = "numeric")
    lo <- lower(data)[colnames(Y)]
    up <- upper(data)[colnames(Y)]
    if (length(attr(mt, which = "term.labels")) == 0L) {
        Z <- datacggm(Y = Y, lo = lo, up = up)
    } else {
        X <- model.matrix(mt, mf, contrasts)[, -1L, drop = FALSE]
        Z <- datacggm(Y = Y, X = X, lo = lo, up = up)
    }
    # starting main code cglasso
    X.null <- is.null(getMatrix(Z, "X"))
    n <- nobs(Z)
    p <- nresp(Z)
    q <- npred(Z)
    if (p == 1L) stop("number of response variables is equal to ", sQuote(1))
    xnames <- colNames(Z)$X
    ynames <- colNames(Z)$Y
    # Testing 'diagonal'
    if (!is.vector(diagonal)) stop(sQuote("diagonal"), " is not a vector")
    if (length(diagonal) > 1L) stop(sQuote("diagonal"), " is not an object of length ", sQuote(1))
    if (!is.logical(diagonal)) stop(sQuote("diagonal"), " is not a logical object")
    # Testing 'weights.B'
    if (!is.null(weights.B)) {
        if (X.null) stop("Argument ", sQuote("weights.B"), " can not be used because ", sQuote("X"), " is missing")
        if (!is.matrix(weights.B)) stop(sQuote("weights.B"), " is not a matrix")
        weights.B.dim <- c(q, p)
        if (any(dim(weights.B) != weights.B.dim)) stop(sQuote("weights.B"), " is not a matrix with 'dim' attribute equal to ", weights.B.dim[1L], " ", weights.B.dim[2L])
        if (any(weights.B < 0)) stop("negative weights in ", sQuote("weights.B "), " are not allowed")
        if (all(weights.B == 0)) stop("all entries in ", sQuote("weights.B "), " are equal to zero")
        weights.B[weights.B == +Inf] <- .Machine$double.xmax
        rownames(weights.B) <- xnames
        colnames(weights.B) <- ynames
    } else {
        if (!X.null) weights.B <- matrix(1, nrow = q, ncol = p, dimnames = list(xnames, ynames))
    }
    # Testing 'weights.Tht'
    if (!is.null(weights.Tht)) {
        if (!is.matrix(weights.Tht)) stop(sQuote("weights.Tht"), " is not a matrix")
        weights.Tht.dim <- c(p, p)
        if (any(dim(weights.Tht) != weights.Tht.dim)) stop(sQuote("weights.Tht"), " is not a matrix with dim attribute equal to ", weights.Tht.dim[1L], " ", weights.Tht.dim[2L])
        if(!isSymmetric(weights.Tht)) stop(sQuote("weights.Tht"), " is not a symmetric matrix")
        if (any(weights.Tht < 0)) stop("negative weights in ", sQuote("weights.Tht "), " are not allowed")
        if (all(weights.Tht == 0)) stop("all entries in ", sQuote("weights.Tht "), " are equal to zero")
        weights.Tht[weights.Tht == +Inf] <- .Machine$double.xmax
        if (!diagonal & any(diag(weights.Tht) != 0)) stop(sQuote("diagonal = FALSE"), " but some diagonal entry in ", sQuote("weights.Tht"), " is not zero")
        if (diagonal & all(diag(weights.Tht) == 0)) stop(sQuote("diagonal = TRUE"), " but all diagonal entries in ", sQuote("weights.Tht"), " are zero")
    } else {
        weights.Tht <- matrix(1, nrow = p, ncol = p)
        diag(weights.Tht) <- ifelse(diagonal, 1, 0)
    }
    rownames(weights.Tht) <- ynames
    colnames(weights.Tht) <- ynames
    # Testing 'nlambda' and 'lambda'
    if (!missing(nlambda) & !missing(lambda))
        warning("Argumnet ", sQuote("nlambda"), " is overwritten using length(lambda)")
    # Testing 'lambda.min.ratio' and 'rho'
    if (!missing(lambda.min.ratio) & !missing(lambda))
        warning("Argumnet ", sQuote("lambda.min.ratio"), " is overwritten using lambda")
    # Testing 'nlambda'
    if (!missing(nlambda)) {
        if (X.null) stop("Argument ", sQuote("nlambda"), " can not be used because ", sQuote("X"), " is missing")
        if (!is.vector(nlambda)) stop(sQuote("nlambda"), " is not a vector")
        if (length(nlambda) != 1) stop(sQuote("nlambda"), " is not an object of length ", sQuote(1))
        if (abs(as.integer(nlambda) - nlambda) > 0) stop(sQuote("nlambda"), " is not an integer value")
        if (nlambda <= 0) stop(sQuote("nlambda"), " is not a strictly positive integer value")
    } else {
        nlambda <- ifelse(X.null, 1L, 10L)
    }
    # Testing 'lambda.min.ratio'
    if (!missing(lambda.min.ratio)) {
        if (X.null) stop("Argument ", sQuote("lambda.min.ratio"), " can not be used because ", sQuote("X"), " is missing")
        if (!is.vector(lambda.min.ratio)) stop(sQuote("lambda.min.ratio"), " is not a vector")
        if (length(lambda.min.ratio) != 1) stop(sQuote("lambda.min.ratio"), " is not an object of length ", sQuote(1))
        if (lambda.min.ratio < 0 | lambda.min.ratio > 1) stop(sQuote("lambda.min.ratio"), " does not belong to the closed interval [0, 1]")
    } else {
        if (X.null) lambda.min.ratio <- 0
        else lambda.min.ratio <- ifelse(q < n, zero, 1.0e-2)
    }
    # Testing 'lambda'
    if (!missing(lambda)) {
        if (X.null) stop("Argument ", sQuote("lambda"), " can not be used because ", sQuote("X"), " is missing")
        if(!is.vector(lambda)) stop(sQuote("lambda"), " is not a vector")
        if(any(lambda < 0)) stop("some entry in ", sQuote("lambda"), " is negative")
        lambda <- sort(lambda, decreasing = TRUE)
        id <- lambda <= zero
        if(any(id)) lambda <- c(lambda[!id], zero)
        if (!is.null(nlambda)) nlambda <- length(lambda)
        if (!is.null(lambda.min.ratio)) lambda.min.ratio <- lambda[nlambda] / lambda[1L]
    } else {
        if (X.null) lambda <- 0
        else lambda <- double(nlambda)
    }
    # Testing 'nrho' and 'rho'
    if (!missing(nrho) & !missing(rho))
        warning("Argumnet ", sQuote("nrho"), " is overwritten using length(rho)")
    # Testing 'rho.min.ratio' and 'rho'
    if (!missing(rho.min.ratio) & !missing(rho))
        warning("Argumnet ", sQuote("rho.min.ratio"), " is overwritten using rho")
    # Testing 'nrho'
    if (!missing(nrho)) {
        if (!is.vector(nrho)) stop(sQuote("nrho"), " is not a vector")
        if (length(nrho) != 1) stop(sQuote("nrho"), " is not an object of length ", sQuote(1))
        if (abs(as.integer(nrho) - nrho) > 0) stop(sQuote("nrho"), " is not an integer value")
        if (nrho <= 0) stop(sQuote("nrho"), " is not a strictly positive integer value")
    } else nrho <- 10L
    # Testing 'rho.min.ratio'
    if (!missing(rho.min.ratio)) {
        if (!is.vector(rho.min.ratio)) stop(sQuote("rho.min.ratio"), " is not a vector")
        if (length(rho.min.ratio) != 1) stop(sQuote("rho.min.ratio"), " is not an object of length ", sQuote(1))
        if (rho.min.ratio < 0 | rho.min.ratio > 1) stop(sQuote("rho.min.ratio"), " does not belong to the closed interval [0, 1]")
    } else rho.min.ratio <- ifelse(p < n, zero, 1.0e-2)
    # Testing 'rho'
    if (!missing(rho)) {
        if(!is.vector(rho)) stop(sQuote("rho"), " is not a vector")
        if(any(rho < 0)) stop("some entry in ", sQuote("rho"), " is negative")
        rho <- sort(rho, decreasing = TRUE)
        id <- rho <= zero
        if(any(id)) rho <- c(rho[!id], zero)
        nrho <- length(rho)
        rho.min.ratio <- rho[nrho] / rho[1L]
    } else rho <- double(nrho)
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
    # fitting cglasso model
    out <- cglasso.fit(Z = Z, diagonal = diagonal, weights.B = weights.B, weights.Tht = weights.Tht, nlambda = nlambda,
                       lambda.min.ratio = lambda.min.ratio, lambda = lambda, nrho = nrho, rho.min.ratio = rho.min.ratio,
                       rho = rho, maxit.em = maxit.em, thr.em = thr.em, maxit.bcd = maxit.bcd, thr.bcd = thr.bcd, trace = trace,
                       algorithm = algorithm)
    if (out$conv != "Ok") {
        msg <- paste("cglasso does not converge. Subroutine", sQuote(out$subrout), "returns message:\n",sQuote(out$conv))
        warning(msg)
    }
    weights.B[weights.B == .Machine$double.xmax] <- +Inf
    weights.Tht[weights.Tht == .Machine$double.xmax] <- +Inf
    InfoStructure <- list(Adj_yy = out$Adj_yy, ncomp = out$ncomp, Ck = out$Ck, pk = out$pk, Adj_xy = out$Adj_xy)
    if (Z$Info$Pattern[1L, "i"] > n) model <- "glasso"
    else {
        id.mar <- any(Z$Info$Pattern[, "nmar"] > 0)
        id.cens <- any(Z$Info$Pattern[, c("nlc", "nrc")] > 0)
        if (id.mar & !id.cens) model <- "missglasso"
        if (!id.mar & id.cens) model <- "censored glasso"
        if (id.mar & id.cens) model <- "hybrid glasso"
    }
    if (q > 0L) model <- paste("conditional", model)
    out.cglasso <- list(call = this.call, Yipt = out$Yipt, B = out$B, mu = out$mu, R = out$R, S = out$S, Sgm = out$Sgm,
                        Tht = out$Tht, dfB = out$dfB, dfTht = out$dfTht, InfoStructure = InfoStructure,
                        nit = out$nit, Z = Z, diagonal = diagonal, weights.B = weights.B, weights.Tht = weights.Tht,
                        nlambda = out$nlambda, lambda.min.ratio = out$lambdaratio, lambda = out$lambda,
                        nrho = out$nrho, rho.min.ratio = out$rhoratio, rho = out$rho, model = model, maxit.em = maxit.em,
                        thr.em = thr.em, maxit.bcd = maxit.bcd, thr.bcd = thr.bcd, conv = out$conv,
                        subrout = out$subrout, trace = trace, nobs = n, nresp = nResp, npred = nPred)
    class(out.cglasso) <- "cglasso"
    out.cglasso
}

# cglasso fitting function
cglasso.fit <- function(Z, diagonal, weights.B, weights.Tht, nlambda, lambda.min.ratio, lambda, nrho, rho.min.ratio,
                        rho, maxit.em, thr.em, maxit.bcd, thr.bcd, trace, algorithm) {
    n <- nobs(Z)
    p <- nresp(Z)
    q <- npred(Z)
    Y <- getMatrix(Z, "Y", ordered = TRUE)
    Y[is.na(Y)] <- 0
    X <- getMatrix(Z, "X", ordered = TRUE)
    if (!is.null(X)) X <- as.matrix(X)
    ynames <- colNames(Z)$Y
    xnames <- colNames(Z)$X
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
    # setting output
    lambda_lab <- paste0("lmb", seq_len(nlambda))
    rho_lab <- paste0("rho", seq_len(nrho))
    Yipt <- array(0, dim = c(n, p, nlambda, nrho), dimnames = list(yrnames, response = ynames, lambda = lambda_lab, rho = rho_lab))
    B <- array(0, dim = c(q + 1L, p, nlambda, nrho), dimnames = list(coef = c("Int.", xnames), response = ynames, lambda = lambda_lab, rho = rho_lab))
    mu <- array(0, dim = c(n, p, nlambda, nrho), dimnames = list(yrnames, response = ynames, lambda = lambda_lab, rho = rho_lab))
    R <- array(0, dim = c(n, p, nlambda, nrho), dimnames = list(yrnames, response = ynames, lambda = lambda_lab, rho = rho_lab))
    S <- array(0, dim = c(p, p, nlambda, nrho), dimnames = list(response = ynames, response = ynames, lambda = lambda_lab, rho = rho_lab))
    Sgm <- array(0, dim = c(p, p, nlambda, nrho), dimnames = list(response = ynames, response = ynames, lambda = lambda_lab, rho = rho_lab))
    Tht <- array(0, dim = c(p, p, nlambda, nrho), dimnames = list(response = ynames, response = ynames, lambda = lambda_lab, rho = rho_lab))
    Adj_yy <- array(0L, dim = c(p, p, nlambda, nrho), dimnames = list(response = ynames, response = ynames, lambda = lambda_lab, rho = rho_lab))
    if (q > 0)
        Adj_xy <- array(0L, dim = c(q, p, nlambda, nrho), dimnames = list(predictor = xnames, response = ynames, lambda = lambda_lab, rho = rho_lab))
    dfB <- array(0L, dim = c(p + 1L, nlambda, nrho), dimnames = list(df = c(ynames, "Tot."), lambda = lambda_lab, rho = rho_lab))
    dfTht <- matrix(0L, nrow = nlambda, ncol = nrho, dimnames = list(lambda = lambda_lab, rho = rho_lab))
    ncomp <- matrix(0L, nrow = nlambda, ncol = nrho, dimnames = list(lambda = lambda_lab, rho = rho_lab))
    Ck <- array(0L, dim = c(p, nlambda, nrho), dimnames = list(NULL, lambda = lambda_lab, rho = rho_lab))
    pk <- array(0L, dim = c(p, nlambda, nrho), dimnames = list(NULL, lambda = lambda_lab, rho = rho_lab))
    nit <- array(0L, dim = c(2L, nlambda, nrho), dimnames = list(steps = c("EM", "nit"), lambda = lambda_lab, rho = rho_lab))
    ########################
    # setting storage.mode #
    ########################
    storage.mode(n) <- "integer"
    storage.mode(q) <- "integer"
    if (q != 0L) storage.mode(X) <- "double"
    storage.mode(p) <- "integer"
    storage.mode(Y) <- "double"
    storage.mode(Id) <- "integer"
    storage.mode(nP) <- "integer"
    storage.mode(InfoP) <- "integer"
    storage.mode(lo) <- "double"
    storage.mode(up) <- "double"
    if (!is.null(weights.B)) storage.mode(weights.B) <- "double"
    storage.mode(diagonal) <- "integer"
    storage.mode(weights.Tht) <- "double"
    storage.mode(ym) <- "double"
    storage.mode(yv) <- "double"
    storage.mode(nlambda) <- "integer"
    storage.mode(lambda.min.ratio) <- "double"
    storage.mode(lambda) <- "double"
    storage.mode(nrho) <- "integer"
    storage.mode(rho.min.ratio) <- "double"
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
    storage.mode(Adj_yy) <- "integer"
    if (q > 0) storage.mode(Adj_xy) <- "integer"
    storage.mode(dfB) <- "integer"
    storage.mode(dfTht) <- "integer"
    storage.mode(ncomp) <- "integer"
    storage.mode(Ck) <- "integer"
    storage.mode(pk) <- "integer"
    storage.mode(nit) <- "integer"
    conv <- integer(1)
    subrout <- integer(1)
    storage.mode(trace) <- "integer"
    # fitting models
    if (q == 0L) {
        # fitting standard cglasso model
        if(algorithm == "glasso"){
            out <- .Fortran(C_cglasso_v1, n = n, p = p, Y = Y, Id = Id, nP = nP, InfoP = InfoP, lo = lo, up = up,
                            pendiag = diagonal, wTht = weights.Tht, ym = ym, yv = yv, nrho = nrho, rhoratio = rho.min.ratio,
                            rho = rho, maxit_em = maxit.em, thr_em = thr.em, maxit_bcd = maxit.bcd, thr_bcd = thr.bcd,
                            Yipt = Yipt, B = B, mu = mu, R = R, S = S, Sgm = Sgm, Tht = Tht, Adj_yy = Adj_yy, dfB = dfB, dfTht = dfTht,
                            ncomp = ncomp, Ck = Ck, pk = pk, nit = nit, conv = conv, subrout = subrout, trace = trace)
        } else {
            out <- .Fortran(C_cglasso_v3, n = n, p = p, Y = Y, Id = Id, nP = nP, InfoP = InfoP, lo = lo, up = up,
                            pendiag = diagonal, wTht = weights.Tht, ym = ym, yv = yv, nrho = nrho, rhoratio = rho.min.ratio,
                            rho = rho, maxit_em = maxit.em, thr_em = thr.em, maxit_bcd = maxit.bcd, thr_bcd = thr.bcd,
                            Yipt = Yipt, B = B, mu = mu, R = R, S = S, Sgm = Sgm, Tht = Tht, Adj_yy = Adj_yy, dfB = dfB, dfTht = dfTht,
                            ncomp = ncomp, Ck = Ck, pk = pk, nit = nit, conv = conv, subrout = subrout, trace = trace)
        }
        out$nlambda <- nlambda
        out$lambdaratio <- lambda.min.ratio
        out$lambda <- lambda
        out$subrout <- switch(as.character(out$subrout),
                              "0" = "",
                              "1" = "e_step",
                              "2" = "glassosub",
                              "3" = "cglasso_v1")
    }
    else {
        # fitting conditional cglasso model
        if(algorithm == "glasso"){
            out <- .Fortran(C_cglasso_v2, n = n, q = q, X = X, p = p, Y = Y, Id = Id, nP = nP, InfoP = InfoP,
                            lo = lo, up = up, wB = weights.B, pendiag = diagonal, wTht = weights.Tht, ym = ym,
                            yv = yv, nlambda = nlambda, lambdaratio = lambda.min.ratio, lambda = lambda,
                            nrho = nrho, rhoratio = rho.min.ratio, rho = rho, maxit_em = maxit.em, thr_em = thr.em,
                            maxit_bcd = maxit.bcd, thr_bcd = thr.bcd, Yipt = Yipt, B = B, mu = mu, R = R, S = S,
                            Sgm = Sgm, Tht = Tht, Adj_yy = Adj_yy, dfB = dfB, dfTht = dfTht, ncomp = ncomp, Ck = Ck,
                            pk = pk, Adj_xy = Adj_xy, nit = nit, conv = conv, subrout = subrout, trace = trace)
        } else {
            out <- .Fortran(C_cglasso_v4, n = n, q = q, X = X, p = p, Y = Y, Id = Id, nP = nP, InfoP = InfoP,
                            lo = lo, up = up, wB = weights.B, pendiag = diagonal, wTht = weights.Tht, ym = ym,
                            yv = yv, nlambda = nlambda, lambdaratio = lambda.min.ratio, lambda = lambda,
                            nrho = nrho, rhoratio = rho.min.ratio, rho = rho, maxit_em = maxit.em, thr_em = thr.em,
                            maxit_bcd = maxit.bcd, thr_bcd = thr.bcd, Yipt = Yipt, B = B, mu = mu, R = R, S = S,
                            Sgm = Sgm, Tht = Tht, Adj_yy = Adj_yy, dfB = dfB, dfTht = dfTht, ncomp = ncomp, Ck = Ck,
                            pk = pk, Adj_xy = Adj_xy, nit = nit, conv = conv, subrout = subrout, trace = trace)
        }
        out$subrout <- switch(as.character(out$subrout),
                              "0" = "",
                              "1" = "e_step",
                              "2" = "multilasso",
                              "3" = "glassosub",
                              "4" = "cglasso_v2")
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
