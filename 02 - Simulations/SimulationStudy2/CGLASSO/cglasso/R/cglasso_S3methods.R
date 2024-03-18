print.cglasso <- function (x, digits = 3L, ...){
    p <- nresp(x$Z)
    q <- npred(x$Z)
    nrho <- x$nrho
    nlambda <- x$nlambda
    rho <- x$rho
    lambda <- x$lambda
    dots <- list(...)
    if (is.null(dots$print.gap)) dots$print.gap <- 2L
    if (is.null(dots$quote)) dots$quote <- FALSE
    if (is.null(dots$row.names)) dots$row.names <- FALSE
    if (nrho == 1L | nlambda == 1L) {
        df.B <- x$dfB[p + 1L, , ] +  p
        df.Tht <- drop(x$dfTht) + p
        df <- df.B + df.Tht
        df.max <- p * (p + 1) / 2 + (q + 1L) * p
        df.per <- formatC(round(df / df.max * 100, digits = digits), format = 'f', digits = digits)
        df.per <- paste("(", df.per, "%)", sep = "")
        ncomp <- drop((x$InfoStructure$ncomp))
        tbl <- data.frame(lambda, rho, df, df.per, ncomp)
        names(tbl) <- c("lambda", "rho", "df", "(df%)", "N. Comp.")
        if (q == 0L) tbl <- tbl[, -1L, drop = FALSE]
        cat("\nCall:  ", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
        do.call(function(...) print.data.frame(tbl, digits = digits, ...), dots)
        cat("\n---\n")
    } else {
        lambda <- rep(lambda, each = nrho)
        rho <- rep(rho, times = nlambda)
        df.B <- as.vector(t(x$dfB[p + 1L, , ])) +  p
        df.Tht <- as.vector(t(x$dfTht)) + p
        df <- df.B + df.Tht
        df.max <- p * (p + 1) / 2 + (q + 1L) * p
        df.per <- formatC(round(df / df.max * 100, digits = digits), format = 'f', digits = digits)
        df.per <- paste("(", df.per, "%)", sep = "")
        ncomp <- as.vector(t(x$InfoStructure$ncomp))
        tbl <- data.frame(lambda, rho, df, df.per, ncomp)
        names(tbl) <- c("lambda", "rho", "df", "", "N. Comp.")
        cat("\nCall:  ", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
        #tbl.list <- split(tbl, f = rep(seq_len(nlambda), each = nrho), drop = FALSE)
        if (nlambda <= nrho) f <- rep(seq_len(nlambda), each = nrho)
        else f <- rep(seq_len(nrho), times = nlambda)
        tbl.list <- split(tbl, f = f, drop = FALSE)
        do.call(function(...) print.listof(tbl.list, digits = digits, ...), dots)
        cat("---\n")
    }
    cat("\nmodel:", sQuote(x$model))
    cat("\n    nObs:", x$nobs)
    cat("\n   nResp:", x$nresp)
    cat("\n   nPred:", x$npred, "\n\n")
    invisible(tbl)
}

coef.cglasso <- function(object, type = c("all", "B", "Sigma", "Theta"), lambda.id, rho.id, drop = TRUE, ...) {
    type <- match.arg(type)
    nrho <- object$nrho
    nlambda <- object$nlambda
    if (!is.logical(drop)) stop(sQuote("drop"), " is not an object of type ", sQuote("integer"))
    if (missing(rho.id)) rho.id <- seq_len(nrho)
    else {
        if(!is.vector(rho.id)) stop(sQuote("rho.id"), " is not a vector")
        if(any(abs(as.integer(rho.id) - rho.id) > 0)) stop(sQuote("rho.id"), " is not an object of type ", dQuote("integer"))
        if(min(rho.id) <= 0) stop("some entry in ", sQuote("rho.id"), " is not a positive integer")
        if(max(rho.id) > nrho) stop("some entry in ", sQuote("rho.id"), " is larger than ", sQuote(nrho))
    }
    if (missing(lambda.id)) lambda.id <- seq_len(nlambda)
    else {
        if(!is.vector(lambda.id)) stop(sQuote("lambda.id"), " is not a vector")
        if(any(abs(as.integer(lambda.id) - lambda.id) > 0)) stop(sQuote("lambda.id"), " is not an object of type ", dQuote("integer"))
        if(min(lambda.id) <= 0) stop("some entry in ", sQuote("lambda.id"), " is not a positive integer")
        if(max(lambda.id) > nlambda) stop("some entry in ", sQuote("lambda.id"), " is larger than ", sQuote(nlambda))
    }
    
    if (type == "all")
        out.coef <- list(B = object$B[, , lambda.id, rho.id, drop = drop],
                         Sigma = object$Sgm[, , lambda.id, rho.id, drop = drop],
                         Theta = object$Tht[, , lambda.id, rho.id, drop = drop])
    else
        out.coef <- switch(type,
                        B = object$B[, , lambda.id, rho.id, drop = drop],
                        Sigma = object$Sgm[, , lambda.id, rho.id, drop = drop],
                        Theta = object$Tht[, , lambda.id, rho.id, drop = drop])
    out.coef
}

fitted.cglasso <- function(object, lambda.id, rho.id, drop = TRUE, ...) {
    nrho <- object$nrho
    nlambda <- object$nlambda
    if (!is.logical(drop)) stop(sQuote("drop"), " is not an object of type ", sQuote("integer"))
    if (missing(rho.id)) rho.id <- seq_len(nrho)
    else {
        if(!is.vector(rho.id)) stop(sQuote("rho.id"), " is not a vector")
        if(any(abs(as.integer(rho.id) - rho.id) > 0)) stop(sQuote("rho.id"), " is not an object of type ", dQuote("integer"))
        if(min(rho.id) <= 0) stop("some entry in ", sQuote("rho.id"), " is not a positive integer")
        if(max(rho.id) > nrho) stop("some entry in ", sQuote("rho.id"), " is larger than ", sQuote(nrho))
    }
    if (missing(lambda.id)) lambda.id <- seq_len(nlambda)
    else {
        if(!is.vector(lambda.id)) stop(sQuote("lambda.id"), " is not a vector")
        if(any(abs(as.integer(lambda.id) - lambda.id) > 0)) stop(sQuote("lambda.id"), " is not an object of type ", dQuote("integer"))
        if(min(lambda.id) <= 0) stop("some entry in ", sQuote("lambda.id"), " is not a positive integer")
        if(max(lambda.id) > nlambda) stop("some entry in ", sQuote("lambda.id"), " is larger than ", sQuote(nlambda))
    }
    out.fitted <- object$mu[, , lambda.id, rho.id, drop = drop]
    out.fitted
}

residuals.cglasso <- function(object, type = c("observed", "working"), lambda.id, rho.id, drop = TRUE, ...) {
    type <- match.arg(type)
    nrho <- object$nrho
    nlambda <- object$nlambda
    if (!is.logical(drop)) stop(sQuote("drop"), " is not an object of type ", sQuote("integer"))
    if (missing(rho.id)) rho.id <- seq_len(nrho)
    else {
        if (!is.vector(rho.id)) stop(sQuote("rho.id"), " is not a vector")
        if (any(abs(as.integer(rho.id) - rho.id) > 0)) stop(sQuote("rho.id"), " is not an object of type ", dQuote("integer"))
        if (min(rho.id) <= 0) stop("some entry in ", sQuote("rho.id"), " is not a positive integer")
        if (max(rho.id) > nrho) stop("some entry in ", sQuote("rho.id"), " is larger than ", sQuote(nrho))
    }
    if (missing(lambda.id)) lambda.id <- seq_len(nlambda)
    else {
        if (!is.vector(lambda.id)) stop(sQuote("lambda.id"), " is not a vector")
        if (any(abs(as.integer(lambda.id) - lambda.id) > 0)) stop(sQuote("lambda.id"), " is not an object of type ", dQuote("integer"))
        if (min(lambda.id) <= 0) stop("some entry in ", sQuote("lambda.id"), " is not a positive integer")
        if (max(lambda.id) > nlambda) stop("some entry in ", sQuote("lambda.id"), " is larger than ", sQuote(nlambda))
    }
    out.residuals <- object$R[, , lambda.id, rho.id, drop = FALSE]
    if (type == "observed") {
        R <- event(object$Z)
        residuals.dim <- dim(out.residuals)
        for (i in 1:residuals.dim[3L]) {
            for (j in 1:residuals.dim[4L]) out.residuals[, , i, j][R != 0] <- NA
        }
    }
    if (drop) out.residuals <- drop(out.residuals)
    out.residuals
}

predict.cglasso <- function(object, type = c("B", "mu", "Sigma", "Theta"), X.new, lambda.new, rho.new, ...) {
    type <- match.arg(type)
    nrho <- object$nrho
    rho <- object$rho
    nlambda <- object$nlambda
    lambda <- object$lambda
    q <- npred(object$Z)
    if (!missing(X.new)) {
        if (q == 0) stop(sQuote("X.new"), " can not be used because no predictors are available in ", sQuote("object"))
        if (type != "mu") stop("type = ", sQuote(type), " can not be used together with ", sQuote("X.new"))
        if (is.vector(X.new)) X.new <- as.matrix(X.new)
        if (!is.matrix(X.new)) stop(sQuote("X.new"), " is not a matrix")
        if (dim(X.new)[2L] != q) stop("wrong dimension in ", sQuote("X.new"),". The number of columns in not equal to ", sQuote(q))
    }
    if (missing(lambda.new)) {
        if (q > 0) stop(sQuote("lambda.new"), " is missing")
        else lambda.new <- lambda
    } else {
        if (q == 0L) stop(sQuote("X.new"), " can not be used because no predictors are available in ", sQuote("object"))
        if (!is.vector(lambda.new)) stop(sQuote("lambda.new"), " is not a vector")
        if (length(lambda.new) != 1L) stop(sQuote("lambda.new"), " is not an object of length ", sQuote(1))
        if (lambda.new < lambda[nlambda]) stop(sQuote("lambda.new"), "is smaller than ", sQuote(lambda[nlambda]))
    }
    if (missing(rho.new)) stop(sQuote("rho.new"), " is missing")
    else {
        if (!is.vector(rho.new)) stop(sQuote("rho.new"), " is not a vector")
        if (length(rho.new) != 1L) stop(sQuote("rho.new"), " is not an object of length ", sQuote(1))
        if (rho.new < rho[nrho]) stop(sQuote("rho.new"), "is smaller than ", sQuote(rho[nrho]))
    }
    if (type == "mu" & !missing(X.new)) type <- "B"
    Min <- switch(type,
                  B = object$B,
                  mu = object$mu,
                  Sigma = object$Sgm,
                  Theta = object$Tht)
    d1 <- dim(Min)[1L]
    d2 <- dim(Min)[2L]
    Mout <- matrix(0, nrow = d1, ncol = d2, dimnames = list(dimnames(Min)[[1L]], dimnames(Min)[[2L]]))
    ########################
    # setting storage.mode #
    ########################
    storage.mode(rho.new) <- "double"
    storage.mode(lambda.new) <- "double"
    storage.mode(nrho) <- "integer"
    storage.mode(rho) <- "double"
    storage.mode(nlambda) <- "integer"
    storage.mode(lambda) <- "double"
    storage.mode(d1) <- "integer"
    storage.mode(d2) <- "integer"
    storage.mode(Min) <- "double"
    storage.mode(Mout) <- "double"
    out.f <- .Fortran(C_predict, newrho = rho.new, newlambda = lambda.new, nrho = nrho,
                      rho = rho, nlambda = nlambda, lambda = lambda, d1 = d1, d2 = d2,
                      Min = Min, Mout = Mout)
    if (type == "B" & !missing(X.new)) out <- cbind(1, X.new) %*% out.f$Mout
    else out <- out.f$Mout
    out
}

summary.cglasso <- function(object, GoF = AIC, print.all = TRUE, digits = 3L,  ...){
    if (!is.element(class(GoF), c("function", "GoF")))
        stop (sQuote("GoF"), " is not a goodness-of-fit function (AIC or BIC) neither an object of class ", sQuote("GoF"))
    dots <- list(...)
    n <- nobs(object$Z)
    p <- nresp(object$Z)
    q <- npred(object$Z)
    if (is.function(GoF)) {
        if (is.null(dots$type)) dots$type <- ifelse(q == 0L, "FD", "CC")
        GoF.name <- deparse(substitute(GoF))
        if (!is.element(GoF.name, c("AIC", "BIC")))
        stop(sQuote(GoF.name), " is not a valid function. Please, use ", sQuote("AIC"), " or ", sQuote("BIC"))
        GoF <- switch(GoF.name,
                      AIC = do.call(function(...) AIC(object, ...), dots),
                      BIC = do.call(function(...) BIC(object, ...), dots))
    }
    if (!is.vector(print.all)) stop(sQuote("print.all"), " is not a vector")
    if (length(print.all) != 1L) stop(sQuote("print.all"), " is not a vector of length ", sQuote("1"))
    if (!is.logical(print.all)) stop(sQuote("print.all"), " is not an object of type ", sQuote("logical"))
    nrho <- object$nrho
    nlambda <- object$nlambda
    rho <- object$rho
    lambda <- object$lambda
    if (is.null(dots$print.gap)) dots$print.gap <- 2L
    if (is.null(dots$quote)) dots$quote <- FALSE
    if (is.null(dots$row.names)) dots$row.names <- FALSE
    if (nrho == 1L | nlambda == 1L) {
        df.B <- object$dfB[p + 1L, , ] +  p
        df.Tht <- drop(object$dfTht) + p
        df <- df.B + df.Tht
        df.max <- p * (p + 1) / 2 + (q + 1L) * p
        df.per <- formatC(round(df / df.max * 100, digits = digits), format = 'f', digits = digits)
        df.per <- paste("(", df.per, "%)", sep = "")
        val <- drop(GoF$value_gof)
        rnk <- rank(val)
        rnk_min <- which.min(rnk)
        rnk <- as.character(rnk)
        rnk[-rnk_min] <- paste(rnk[-rnk_min], "  ")
        rnk[rnk_min] <- paste(rnk[rnk_min], "<-")
        lambda.id <- ifelse (nlambda == 1L, 1L, rnk_min)
        rho.id <- ifelse (nrho == 1L, 1L, rnk_min)
        tbl <- data.frame(lambda, rho, df.B, df.Tht, df, df.per, val, Rank = rnk)
        names(tbl)[6L] <- "(df%)"
        names(tbl)[7L] <- GoF$type
        if (q == 0L) tbl <- tbl[, -1L, drop = FALSE]
        cat("\nCall:  ", paste(deparse(object$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
        if (print.all) do.call(function(...) print.data.frame(tbl, digits = digits, ...), dots)
    } else {
        lambda <- rep(lambda, each = nrho)
        lambda.id <- rep(seq_len(nlambda), each = nrho)
        rho <- rep(rho, times = nlambda)
        rho.id <- rep(seq_len(nrho), times = nlambda)
        df.B <- as.vector(t(object$dfB[p + 1L, , ])) +  p
        df.Tht <- as.vector(t(object$dfTht)) + p
        df <- df.B + df.Tht
        df.max <- p * (p + 1) / 2 + (q + 1L) * p
        df.per <- formatC(round(df / df.max * 100, digits = digits), format = 'f', digits = digits)
        df.per <- paste("(", df.per, "%)", sep = "")
        val <- as.vector(t(GoF$value_gof))
        rnk <- rank(val)
        rnk_min <- which.min(rnk)
        rnk <- as.character(rnk)
        rnk[-rnk_min] <- paste(rnk[-rnk_min], "  ")
        rnk[rnk_min] <- paste(rnk[rnk_min], "<-")
        lambda.id <- lambda.id[rnk_min]
        rho.id <- rho.id[rnk_min]
        tbl <- data.frame(lambda, rho, df.B, df.Tht, df, df.per, val, Rank = rnk)
        names(tbl)[6L] <- "(df%)"
        names(tbl)[7L] <- GoF$type
        cat("\nCall:  ", paste(deparse(object$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
        tbl.list <- split(tbl, f = rep(seq_len(nlambda), each = nrho))
        if (print.all) do.call(function(...) print.listof(tbl.list, digits = digits, ...), dots)
    }
    lbl <- c("model", "nObs", "nResp", "nPred", "lambda", "lambda.id", "rho", "rho.id", GoF$type, "df.B", "df.Tht", "df")
    lbl <- paste("\n", format(lbl, justify = "right"), ":", sep = "")
    cat("\n===============================================================")
    cat("\n\nSummary of the Selected Model\n")
    if (!is.null(object$model)) cat(lbl[1L], sQuote(object$model))
    cat(lbl[2L], n)
    cat(lbl[3L], object$nresp)
    cat(lbl[4L], object$npred)
    if (q > 0L) {
        cat(lbl[5L], tbl[rnk_min, "lambda"])
        cat(lbl[6L], lambda.id)
    }
    cat(lbl[7L], tbl[rnk_min, "rho"])
    cat(lbl[8L], rho.id)
    cat(lbl[9L], tbl[rnk_min, GoF$type])
    cat(lbl[10L], tbl[rnk_min, "df.B"])
    cat(lbl[11L], tbl[rnk_min, "df.Tht"])
    cat(lbl[12L], tbl[rnk_min, "df"])
    cat("\n\n===============================================================\n\n")
    invisible(list(table = tbl, lambda.id = lambda.id, rho.id = rho.id))
}

plot.cglasso <- function(x, what = c("Theta", "diag(Theta)", "b0", "B"), penalty = ifelse(x$nrho >= x$nlambda, "rho", "lambda"),
                         given = NULL, GoF = AIC, add.labels, matplot.arg1, matplot.arg2, labels.arg, abline.arg, mtext.arg, save.plot,
                         grdev = pdf, grdev.arg, digits = 4L, ...) {
    p <- nresp(x$Z)
    q <- npred(x$Z)
    nlambda <- x$nlambda
    nrho <- x$nrho
    if (nrho == 1L & nlambda == 1L) stop("plot method is not available because ", sQuote("nlambda = 1"), " and ", sQuote("nrho = 1"))
    if (inherits(what, "formula")) {
        xy <- function2xy(what)
        penalty <- xy[["penalty"]]
        what <- xy[["what"]]
        given <- xy[["given"]]
    }
    # testing arguments
    penalty <- match.arg(penalty, c("rho", "lambda"))
    if (penalty == "lambda" & nlambda == 1L) stop("plot method is not available because ", sQuote("nlambda = 1"))
    if (penalty == "rho" & nrho == 1L) stop("plot method is not available because ", sQuote("nrho = 1"))
    what <- match.arg(what)
    if (what == "B" & q == 0L) stop("plot method is not available because the number of predictors is zero")
    if (is.null(given)) given <- seq_len(ifelse(penalty == "rho", nlambda, nrho))
    else {
        given.max <- ifelse(penalty == "rho", nlambda, nrho)
        if (!is.vector(given)) stop (sQuote("given"), "is not a vector")
        if (any(abs(as.integer(given) - given) > 0)) stop(sQuote("given"), " is not an object of type ", dQuote("integer"))
        if (min(given) <= 0) stop("some entry in ", sQuote("given"), " is not a positive integer")
        if (any(given > given.max)) stop("some entry in ", sQuote("given"), " is larger than ", sQuote(given.max))
    }
    ngiven <- length(given)
    # testing 'matplot.arg1'
    if (missing(matplot.arg1)) matplot.arg1 <- vector(mode = "list")
    if (is.null(matplot.arg1$type)) matplot.arg1$type <- "l"
    if (is.null(matplot.arg1$col)) matplot.arg1$col <- "black"
    if (is.null(matplot.arg1$lty)) matplot.arg1$lty <- 1L
    if (!is.list(matplot.arg1)) stop(sQuote("matplot.arg1"), " is not an object of type ", dQuote("list"))
    if (is.null(names(matplot.arg1))) stop(sQuote("matplot.arg1"), " is not a naned list")
    # testing 'matplot.arg2'
    if (missing(matplot.arg2)) matplot.arg2 <- vector(mode = "list")
    if (is.null(matplot.arg2$type)) matplot.arg2$type <- "l"
    if (is.null(matplot.arg2$col)) matplot.arg2$col <- "gray70"
    if (is.null(matplot.arg2$lty)) matplot.arg2$lty <- 2L
    if (!is.list(matplot.arg2)) stop(sQuote("matplot.arg2"), " is not an object of type ", dQuote("list"))
    if (is.null(names(matplot.arg2))) stop(sQuote("matplot.arg2"), " is not a naned list")
    # testing 'labels.arg'
    if (missing(labels.arg)) labels.arg <- vector(mode = "list")
    if (is.null(labels.arg$pos)) labels.arg$pos <- 4L
    if (!is.list(labels.arg)) stop(sQuote("labels.arg"), " is not an object of type ", dQuote("list"))
    if (is.null(names(labels.arg))) stop(sQuote("labels.arg"), " is not a naned list")
    # testing 'abline.arg'
    if (missing(abline.arg)) abline.arg <- vector(mode = "list")
    if (is.null(abline.arg$lwd)) abline.arg$lwd <- 2L
    if (is.null(abline.arg$lty)) abline.arg$lty <- 2L
    if (is.null(abline.arg$col)) abline.arg$col <- 2L
    if (!is.list(abline.arg)) stop(sQuote("abline.arg"), " is not an object of type ", dQuote("list"))
    if (is.null(names(abline.arg))) stop(sQuote("abline.arg"), " is not a naned list")
    # testing 'mtext.arg'
    if (missing(mtext.arg)) mtext.arg <- vector(mode = "list")
    else {
        if (!is.list(mtext.arg)) stop(sQuote("mtext.arg"), " is not an object of type ", dQuote("list"))
        if (is.null(names(mtext.arg))) stop(sQuote("mtext.arg"), " is not a naned list")
    }
    # testing 'save.plot'
    if (missing(save.plot)) save.plot <- FALSE
    if (length(save.plot) != 1L) stop(sQuote(save.plot), " is not an object of length ", sQuote("1"))
    if (!inherits(save.plot, c("logical", "character"))) stop(sQuote(save.plot), " is not an object of type ", sQuote("logical"), "or ", , sQuote("character"))
    if(inherits(save.plot, "character")) {
        oldPath <- getwd()
        newPath <- save.plot
        setwd(newPath)
        on.exit(setwd(oldPath), add = TRUE, after = TRUE)
        save.plot <- TRUE
    }
    # testing 'grdev'
    if (!is.function(grdev)) stop(sQuote("grdev"), " is not a function")
    else {
        grdev.name <- deparse(substitute(grdev))
        if (!is.element(grdev.name, c("pdf", "postscript", "svg", "bmp", "jpeg", "png", "tiff")))
            stop(sQuote(grdev.name), " is not a valid graphics device for exporting plots. Please, use one of the following functions:\n",
                 sQuote("pdf"), ", ", sQuote("postscript"), ", ", sQuote("svg"), ", ", sQuote("bmp"), ", ",
                 sQuote("jpeg"), ", ", sQuote("png"), " or ", sQuote("tiff"), ", ")
        if (grdev.name == "postscript") grdev.name <- "ps"
    }
    if (missing(grdev.arg)) grdev.arg <- NULL
    else {
        if (!is.list(grdev.arg)) stop(sQuote("grdev.arg"), " is not an object of type ", sQuote("list"))
    }
    # testing 'digits'
    if (!is.vector(digits)) stop(sQuote("digits"), "is not a vector")
    if (length(digits) > 1L) stop(sQuote("digits"), "is not a vector of length ", sQuote(1))
    if (abs(as.integer(digits) - digits) > 0) stop(sQuote("digits"), " is not an object of type ", dQuote("integer"))
    if (digits < 0L) stop(sQuote("digits"), " is not a positive integer")
    # section: goodness-of-fit
    dots <- list(...)
    if (is.function(GoF)) {
        if (is.null(dots$type)) dots$type <- ifelse(q == 0L, "FD", "CC")
        GoF.name <- deparse(substitute(GoF))
        if (!is.element(GoF.name, c("AIC", "BIC")))
            stop(sQuote(GoF.name), " is not a valid function. Please, use ", sQuote("AIC"), " or ", sQuote("BIC"))
        GoF <- switch(GoF.name,
                      AIC = do.call(function(...) AIC(x, ...), dots),
                      BIC = do.call(function(...) BIC(x, ...), dots))
    }
    # testing 'add.labels'
    if (missing(add.labels)) add.labels <- !is.null(GoF)
    if (!is.logical(add.labels)) stop(sQuote(add.labels), " is not an object of type ", sQuote("logical"))
    if (length(add.labels) != 1L) stop(sQuote(add.labels), " is not an object of length ", sQuote("1"))
    if (add.labels & is.null(GoF)) warning("labels are not added to the plot because GoF is NULL" )
    # starting main code
    tp <- x[[penalty]]
    tp.given <- x[[ifelse(penalty == "rho", "lambda", "rho")]][given]
    tp.given <- formatC(round(tp.given, digits = digits), format = 'f', digits = digits)
    coef.name <- what
    if (coef.name == "diag(Theta)") coef.name <- "Theta"
    if (coef.name == "b0") coef.name <- "B"
    if(save.plot) {
        nplot <- ifelse(what != "B", ngiven, ngiven * p)
        file.names <- paste0(what, "_path_", formatC(seq_len(nplot), width = nchar(nplot), format = "d", flag = "0"), ".", grdev.name)
        cat("Exporting plots\n")
        pb <- txtProgressBar(min = 0L, max = nplot, style = 3L)
        on.exit(close(pb), add = TRUE, after = TRUE)
    } else {
        if (ngiven > 1L | what == "B") {
            op <- par(no.readonly = TRUE)
            on.exit(par(op), add = TRUE)
            devAskNewPage(TRUE)
        }
    }
    if (what == "Theta") U <- .row(c(p, p)) < .col(c(p, p))
    if (add.labels) {
        if (what == "b0") lbls <- colNames(x$Z)$Y
        if (what == "B") lbls <- colNames(x$Z)$X
        if (what == "diag(Theta)") {
            lbls <- seq_len(p)
            lbls <- paste(lbls, lbls, sep = ",")
            lbls <- mapply(function(x) bquote(hat(theta)[.(x)]), lbls)
            names(lbls) <- NULL
            lbls <- sapply(lbls, as.expression)
        }
        if (what == "Theta") {
            lbls <- paste(unlist(mapply(seq, from = 1, to = 1:(p-1))), unlist(mapply(rep, x = 2:p, each = 1:(p-1))), sep = ",")
            lbls <- mapply(function(x) bquote(hat(theta)[.(x)]), lbls)
            names(lbls) <- NULL
            lbls <- sapply(lbls, as.expression)
        }
    }
    hh <- 0
    for (m in seq_len(ngiven)) {
        if (penalty == "lambda") {
            tp.lab <- bquote(lambda ~ " | {" ~ rho %~~% .(tp.given[m]) ~ "}")
            Y <- coef.cglasso(x, type = coef.name, rho.id = given[m], drop = TRUE)
            if (!is.null(GoF))  minGoF.id <- which.min(GoF$value_gof[, given[m]])
        }
        else {
            if (q != 0) tp.lab <- bquote(rho ~ " | {" ~ lambda %~~% .(tp.given[m]) ~ "}")
            else tp.lab <- bquote(rho)
            Y <- coef.cglasso(x, type = coef.name, lambda.id = given[m], drop = TRUE)
            if (!is.null(GoF))   minGoF.id <- which.min(GoF$value_gof[given[m], ])
        }
        ##################
        # plotting 'b0'
        if (what == "b0") {
            hh <- hh + 1L
            if (save.plot) {
                if (is.null(grdev.arg)) grdev(file = file.names[hh])
                else do.call(function(...) grdev(file = file.names[hh], ...), grdev.arg)
            } else dev.hold()
            if (q == 0L) yy <- t(Y)
            else yy <- t(Y[1L, , ])
            xlab <- if (is.null(matplot.arg1$xlab)) tp.lab else matplot.arg1$xlab
            ylab <- if (is.null(matplot.arg1$ylab))  "" else matplot.arg1$ylab
            main <- if (is.null(matplot.arg1$main)) ifelse(q == 0L, "Expected Values", "Intercepts") else matplot.arg1$main
            arg.name <- setdiff(names(matplot.arg1), c("x", "y", "xlab", "ylab", "main"))
            do.call(what = function(...) matplot(x = tp, y = yy, xlab = xlab, ylab = ylab, main = main, ...), args = matplot.arg1[arg.name])
            if (!is.null(GoF)) {
                if(add.labels) {
                    labels <- if(is.null(labels.arg$labels)) lbls else labels.arg$labels
                    arg.name <- setdiff(names(labels.arg), "labels")
                    do.call(what = function(...) text(x = tp[minGoF.id], y = yy[minGoF.id, ], labels = labels, ...), args = labels.arg[arg.name])
                }
                #                abline.arg$v <- tp[minGoF.id]
                #                do.call(what = abline, args = abline.arg)
                abline.arg$v <- tp[minGoF.id]
                do.call(what = abline, args = abline.arg)
                mtext.arg$text <- GoF$type
                mtext.arg$at <- tp[minGoF.id]
                do.call(what = mtext, args = mtext.arg)
                
            }
            if (save.plot) setTxtProgressBar(pb, hh) else dev.flush()
            if(save.plot) dev.off()
        }
        ##################
        # plotting 'B'
        if (what == "B") {
            varname <- colNames(x$Z)$Y
            for (k in seq_len(p)) {
                hh <- hh + 1L
                if (save.plot) {
                    if (is.null(grdev.arg)) grdev(file = file.names[hh])
                    else do.call(function(...) grdev(file = file.names[hh], ...), grdev.arg)
                } else dev.hold()
                yy <- Y[-1L, k, ]
                if (!is.vector(yy)) yy <- t(yy)
                xlab <- if (is.null(matplot.arg1$xlab)) tp.lab else matplot.arg1$xlab
                ylab <- if (is.null(matplot.arg1$ylab)) "Regression Coefficients" else matplot.arg1$ylab
                main <- if (is.null(matplot.arg1$main)) paste0("Response Variable ", varname[k]) else matplot.arg1$main
                arg.name <- setdiff(names(matplot.arg1), c("x", "y", "xlab", "ylab", "main"))
                if (is.null(GoF)) do.call(what = function(...) matplot(x = tp, y = yy, xlab = xlab, ylab = ylab, main = main, ...), args = matplot.arg1[arg.name])
                else {
                    matplot(x = tp, y = yy, xlab = xlab, ylab = ylab, main = main, type = "n")
                    A <- abs(yy[minGoF.id, ]) > 0
                    if(any(!A)) do.call(what = function(...) matpoints(x = tp, y = yy[, !A], ...), args = matplot.arg2[arg.name])
                    if(any(A)) {
                        do.call(what = function(...) matpoints(x = tp, y = yy[, A], ...), args = matplot.arg1[arg.name])
                        if(add.labels) {
                            labels <- if (is.null(labels.arg$labels)) lbls[A] else labels.arg$labels
                            arg.name <- setdiff(names(labels.arg), "labels")
                            do.call(what = function(...) text(x = tp[minGoF.id], y = yy[minGoF.id, A], labels = labels, ...), args = labels.arg[arg.name])
                        }
                    }
                    abline.arg$v <- tp[minGoF.id]
                    do.call(what = abline, args = abline.arg)
                    mtext.arg$text <- GoF$type
                    mtext.arg$at <- tp[minGoF.id]
                    do.call(what = mtext, args = mtext.arg)
                }
                if (save.plot) setTxtProgressBar(pb, hh) else dev.flush()
                if(save.plot) dev.off()
            }
        }
        ##################
        # plotting 'diag(Theta)'
        if (what == "diag(Theta)") {
            hh <- hh + 1L
            if (save.plot) {
                if (is.null(grdev.arg)) grdev(file = file.names[hh])
                else do.call(function(...) grdev(file = file.names[hh], ...), grdev.arg)
            } else dev.hold()
            yy <- t(apply(Y, 3L, diag))
            xlab <- if (is.null(matplot.arg1$xlab)) tp.lab else matplot.arg1$xlab
            ylab <- if (is.null(matplot.arg1$ylab)) "Diagonal Elements" else matplot.arg1$ylab
            main <- if (is.null(matplot.arg1$main)) "Precision Matrix" else matplot.arg1$main
            arg.name <- setdiff(names(matplot.arg1), c("x", "y", "xlab", "ylab", "main"))
            do.call(what = function(...) matplot(x = tp, y = yy, xlab = xlab, ylab = ylab, main = main, ...), args = matplot.arg1[arg.name])
            if (!is.null(GoF)) {
                if(add.labels) {
                    labels <- if(is.null(labels.arg$labels)) lbls else labels.arg$labels
                    arg.name <- setdiff(names(labels.arg), "labels")
                    do.call(what = function(...) text(x = tp[minGoF.id], y = yy[minGoF.id, ], labels = labels, ...), args = labels.arg[arg.name])
                }
                #                abline.arg$v <- tp[minGoF.id]
                #                do.call(what = abline, args = abline.arg)
                abline.arg$v <- tp[minGoF.id]
                do.call(what = abline, args = abline.arg)
                mtext.arg$text <- GoF$type
                mtext.arg$at <- tp[minGoF.id]
                do.call(what = mtext, args = mtext.arg)
            }
            if (save.plot) setTxtProgressBar(pb, hh) else dev.flush()
            if(save.plot) dev.off()
        }
        ##################
        # plotting 'Theta'
        if (what == "Theta") {
            hh <- hh + 1L
            if (save.plot) {
                if (is.null(grdev.arg)) grdev(file = file.names[hh])
                else do.call(function(...) grdev(file = file.names[hh], ...), grdev.arg)
            } else dev.hold()
            yy <- t(apply(Y, 3L, function(M) M[U]))
            xlab <- if (is.null(matplot.arg1$xlab)) tp.lab else matplot.arg1$xlab
            ylab <- if (is.null(matplot.arg1$ylab)) "Off-Diagonal Elements" else matplot.arg1$ylab
            main <- if (is.null(matplot.arg1$main)) "Precision Matrix" else matplot.arg1$main
            arg.name <- setdiff(names(matplot.arg1), c("x", "y", "xlab", "ylab", "main"))
            if (is.null(GoF)) do.call(what = function(...) matplot(x = tp, y = yy, xlab = xlab, ylab = ylab, main = main, ...), args = matplot.arg1[arg.name])
            else {
                matplot(x = tp, y = yy, xlab = xlab, ylab = ylab, main = main, type = "n")
                A <- abs(yy[minGoF.id, ]) > 0
                if(any(!A)) do.call(what = function(...) matpoints(x = tp, y = yy[, !A], ...), args = matplot.arg2[arg.name])
                if(any(A)) {
                    do.call(what = function(...) matpoints(x = tp, y = yy[, A], ...), args = matplot.arg1[arg.name])
                    if(add.labels) {
                        labels <- if(is.null(labels.arg$labels)) lbls[A] else labels.arg$labels
                        arg.name <- setdiff(names(labels.arg), "labels")
                        do.call(what = function(...) text(x = tp[minGoF.id], y = yy[minGoF.id, A], labels = labels, ...), args = labels.arg[arg.name])
                    }
                }
                abline.arg$v <- tp[minGoF.id]
                do.call(what = abline, args = abline.arg)
                mtext.arg$text <- GoF$type
                mtext.arg$at <- tp[minGoF.id]
                do.call(what = mtext, args = mtext.arg)
            }
            if (save.plot) setTxtProgressBar(pb, hh) else dev.flush()
            if(save.plot) dev.off()
        }
    }
    invisible(NULL)
}
