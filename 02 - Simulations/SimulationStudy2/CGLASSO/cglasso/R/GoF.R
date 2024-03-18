QFun <- function(object, mle, verbose = FALSE, ...) {
    # testing class attribute
    if (!inherits(object, c("cggm", "cglasso"))) stop(sQuote("QFun"), "function is not available for an object of class ", sQuote(class(object)))
    # Testing 'mle'
    if (missing(mle)) mle <- FALSE
    else {
        if (class(object)[1L] == "cggm" & !mle) stop("For an object of class ", sQuote("cggm"), " argument ", sQuote("mle"), " must be equal to TRUE")
    }
    if (!is.vector(mle)) stop(sQuote("mle"), " is not a vector")
    if (length(mle) > 1L) stop(sQuote("mle"), " is not an object of length ", sQuote(1))
    if (!is.logical(mle)) stop(sQuote("mle"), " is not a logical object")
    if (class(object)[1L] == "cggm" & mle) mle <- FALSE
    n <- nobs(object$Z)
    p <- nresp(object$Z)
    q <- npred(object$Z)
    Adj_yy <- object$InfoStructure$Adj_yy
    Adj_xy <- object$InfoStructure$Adj_xy
    lambda <- object$lambda
    nlambda <- object$nlambda
    rho <- object$rho
    nrho <- object$nrho
    lambda_lab <- paste0("lmb", seq_len(nlambda))
    rho_lab <- paste0("rho", seq_len(nrho))
    value <- matrix(0, nrow = nlambda, ncol = nrho, dimnames = list(lambda_lab, rho_lab))
    df <- matrix(0, nrow = nlambda, ncol = nrho, dimnames = list(lambda_lab, rho_lab))
    dfB <- matrix(0, nrow = nlambda, ncol = nrho, dimnames = list(lambda_lab, rho_lab))
    dfTht <- matrix(0, nrow = nlambda, ncol = nrho, dimnames = list(lambda_lab, rho_lab))
    const1 <- - 0.5 * n * p * log(2 * pi)
    const2 <- 0.5 * n
    ntp <- nlambda * nrho
    ii <- cbind(i = rep(seq_len(nlambda), each = nrho), j = rep(seq_len(nrho), times = nlambda))
    if(verbose) {
        cat("\nComputing Q-values of the fitted", sQuote(object$model[1L]), ifelse(nrho == 1L & nlambda == 1L, "model", "models"), "\n")
        pb <- txtProgressBar(min = 0, max = ntp, style = 3L)
    }
    for (h in seq_len(ntp)) {
        if(verbose) setTxtProgressBar(pb, h)
        NEW <- TRUE
        i <- ii[h, 1L]
        j <- ii[h, 2L]
        if (mle) {
            for (hh in rev(seq_len(h - 1L))) {
                i1 <- ii[hh, 1L]
                j1 <- ii[hh, 2L]
                test1 <- all(Adj_yy[, , i, j] == Adj_yy[, , i1, j1])
                if (is.null(Adj_xy)) test2 <- TRUE
                else test2 <- all(Adj_xy[, , i, j] == Adj_xy[, , i1, j1])
                if (test1 & test2) {
                    NEW <- FALSE
                    break
                }
            }
            if (NEW) {
                out.mle <- cggm(object, lambda.id = i, rho.id = j, ...)
                Tht <- drop(out.mle$Tht)
                S <- drop(out.mle$S)
                value[i, j] <- const1 + const2 * (determinant(Tht)$modulus - sum(S * Tht))
                dfB[i, j] <- out.mle$dfB[p + 1L, 1L, 1L]
                dfTht[i, j] <- drop(out.mle$dfTht)
                df[i, j] <- dfB[i, j] + dfTht[i, j] + 2 * p
            } else {
                value[i, j] <- value[i1, j1]
                dfB[i, j] <- dfB[i1, j1]
                dfTht[i, j] <- dfTht[i1, j1]
                df[i, j] <- df[i1, j1]
            }
        } else {
            Tht <- object$Tht[, , i, j]
            S <- object$S[, , i, j]
            value[i, j] <- const1 + const2 * (determinant(Tht)$modulus - sum(S * Tht))
            dfB[i, j] <- object$dfB[p + 1L, i, j]
            dfTht[i, j] <- object$dfTht[i, j]
            df[i, j] <- dfB[i, j] + dfTht[i, j] + 2 * p
        }
    }
    if(verbose) close(pb)
    out <- list(value = value, df = df, dfB = dfB, dfTht = dfTht, n = n, p = p, q = q,
                lambda = lambda, nlambda = nlambda, rho = rho, nrho = nrho, model = object$model)
    class(out) <- "QFun"
    out
}

AIC.cglasso <- function(object, k = 2, mle, ...){
    if(k <= 0) stop(sQuote("k"), " is not a positive value")
    type <- ifelse(k == 2, "AIC", "GoF")
    out_QFun <- QFun(object, mle, ...)
    value <-  out_QFun$value
    df <- out_QFun$df
    val <- -2 * value + k * df
    out <- list(value_gof = val, df = df, dfB = out_QFun$dfB, dfTht = out_QFun$dfTht, value = value,
                n = out_QFun$n, p = out_QFun$p, q = out_QFun$q, lambda = out_QFun$lambda, nlambda = out_QFun$nlambda,
                rho = out_QFun$rho, nrho = out_QFun$nrho, type = type, model = out_QFun$model)
    class(out) <- "GoF"
    out
}

BIC.cglasso <- function(object, g = 0, type, mle, ...) {
    # testing 'g'
    if (!is.vector(g)) stop(sQuote("g"), " is not a vector")
    if (length(g) != 1) stop(sQuote("g"), " is not a vector of length ", sQuote(1))
    if ((g < 0) | (g > 1)) stop(sQuote("g"), " does not belong to the closed interval [0, 1]")
    if (g == 0) {
        # classical BIC measure
        out <- AIC(object, k = log(nobs(object$Z)), mle, ...)
        out$type <- "BIC"
    } else {
        out_QFun <- QFun(object, mle, ...)
        n <- out_QFun$n
        p <- out_QFun$p
        q <- out_QFun$q
        value <-  out_QFun$value
        df <- out_QFun$df
#        dfTht <- out_QFun$dfTht
        if (missing(type)) type <- ifelse(q == 0L, "FD", "CC")
        else {
            if (!is.vector(type)) stop(sQuote("type"), "is not a vector")
            if (length(type) != 1L) stop(sQuote("type"), "is not a vector of length ", sQuote("1"))
            if (!is.character(type)) stop(sQuote("type"), "is not an object of type ", sQuote("character"))
            if (!is.element(type, c("FD", "CC")))
                stop(sQuote(type), " is not available. Please, set ", sQuote("type"), " argument equal to ", sQuote("FD"), " or ", sQuote("CC"))
            if (type == "CC" & q == 0L) stop("measure proposed in Chen and Chen (2008, 2012) can not be computed because q = 0")
        }
        if (type == "FD")
            # ebic measure as proposed in Foygel and Drton (2010)
#            val <- -2 * value + dfTht * (log(n) + 4 * g * log(p))
            val <- -2 * value + df * (log(n) + 4 * g * log(p))
        else
            # generalization based on the ebic measure proposed in Chen and Chen (2008)
            # See also Chen and Chen (2012) Statistica Sinica 22, pg. 555-574
            val <- -2 * value + df * (log(n) + 2 * g * log(q))
        type <- paste0("eBIC_", type)
        out <- list(value_gof = val, df = df, dfB = out_QFun$dfB, dfTht = out_QFun$dfTht, value = value,
                    n = out_QFun$n, p = out_QFun$p, q = out_QFun$q, lambda = out_QFun$lambda, nlambda = out_QFun$nlambda,
                    rho = out_QFun$rho, nrho = out_QFun$nrho, type = type, model = out_QFun$model)
    }
    class(out) <- "GoF"
    out
}

select.cglasso <- function(object, GoF = AIC, ...){
    if (!is.element(class(GoF), c("function", "GoF")))
        stop (sQuote("GoF"), " is not either a goodness-of-fit function (AIC or BIC) neither an object of class ", sQuote("GoF"))
    dots <- list(...)
    n <- nobs(object$Z)
    p <- nresp(object$Z)
    q <- npred(object$Z)
    nrho <- object$nrho
    nlambda <- object$nlambda
    if (is.function(GoF)) {
        if (is.null(dots$type)) dots$type <- ifelse(q == 0L, "FD", "CC")
        GoF.name <- deparse(substitute(GoF))
        if (!is.element(GoF.name, c("AIC", "BIC")))
            stop(sQuote(GoF.name), " is not a valid function. Please, use ", sQuote("AIC"), " or ", sQuote("BIC"))
        GoF <- switch(GoF.name,
                        AIC = do.call(function(...) AIC(object, ...), dots),
#                        BIC = do.call(function(...) BIC(object, ...), dots),
                        BIC = do.call(function(...) BIC(object, ...), dots))
    }
    val <- which(GoF$value_gof == min(GoF$value_gof), arr.ind = TRUE)
    if (is.matrix(val)) val <- val[dim(val)[1L], ]
    lambda.id <- val[1L]
    rho.id <- val[2L]
    # reshape cglasso object
    object$Yipt <- object$Yipt[ , , lambda.id, rho.id, drop = FALSE]
    object$B <- object$B[ , , lambda.id, rho.id, drop = FALSE]
    object$mu <- object$mu[ , , lambda.id, rho.id, drop = FALSE]
    object$R <- object$R[ , , lambda.id, rho.id, drop = FALSE]
    object$S <- object$S[ , , lambda.id, rho.id, drop = FALSE]
    object$Sgm <- object$Sgm[ , , lambda.id, rho.id, drop = FALSE]
    object$Tht <- object$Tht[ , , lambda.id, rho.id, drop = FALSE]
    object$dfB <- object$dfB[ , lambda.id, rho.id, drop = FALSE]
    object$dfTht <- object$dfTht[lambda.id, rho.id, drop = FALSE]
    object$InfoStructure$Adj_yy <- object$InfoStructure$Adj_yy[ , , lambda.id, rho.id, drop = FALSE]
    object$InfoStructure$ncomp <- object$InfoStructure$ncomp[lambda.id, rho.id, drop = FALSE]
    object$InfoStructure$Ck <- object$InfoStructure$Ck[ , lambda.id, rho.id, drop = FALSE]
    object$InfoStructure$pk <- object$InfoStructure$pk[ , lambda.id, rho.id, drop = FALSE]
    object$InfoStructure$Adj_xy <- object$InfoStructure$Adj_xy[ , , lambda.id, rho.id, drop = FALSE]
    if (class(object)[1L] == "cglasso") object$nit <- object$nit[ , lambda.id, rho.id, drop = FALSE]
    object$nlambda <- 1L
    object$lambda.min.ratio <- 1
    object$lambda <- object$lambda[lambda.id]
    object$nrho <- 1L
    object$rho.min.ratio <- 1
    object$rho <- object$rho[rho.id]
    object
}
