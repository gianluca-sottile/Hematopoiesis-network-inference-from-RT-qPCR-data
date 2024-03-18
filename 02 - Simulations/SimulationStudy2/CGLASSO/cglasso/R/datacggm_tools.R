fit_matrix <- function(x, n, width) {
    x <- head(x, n = n)
    cumwidth <- ifelse(is.null(rownames(x)), nchar(n) + 4L, max(nchar(rownames(x))) + 1L)
    col.width <- apply(x, 2L, function(x) max(nchar(x)))
    if (!is.null(colnames(x))) col.width <- pmax(nchar(colnames(x)), col.width)
    col.width <- col.width + 1L
    cumwidth <- cumwidth + cumsum(col.width)
    x <- x[, cumwidth <= width, drop = FALSE]
    x
}

getMatrix <- function (x, name = c("Y", "X", "both"), ordered = FALSE) {
    if (!is.datacggm(x)) stop(sQuote(x), " is not an object of class ", sQuote("datacggm"))
    if (!is.logical(ordered)) stop(sQuote(ordered), " id not an object of type ", sQuote("logical"))
    name <- match.arg(name)
    if (ordered) {
        row.order <- x$Info$order
        Y <- x$Y[row.order, , drop = FALSE]
        X <- x$X[row.order, , drop = FALSE]
    } else {
        Y <- x$Y
        X <- x$X
    }
    out <- switch(name,
                  Y = Y,
                  X = X,
                  both = list(Y = Y, X = X))
    out
}

event <- function(x, ordered = FALSE) {
    if(!is.datacggm(x)) stop("'x' is not an object of class 'datacggm'")
    if (!is.logical(ordered)) stop(sQuote(ordered), " id not an object of type ", sQuote("logical"))
    R <- x$Info$R
    if (ordered) R <- R[x$Info$order, , drop = FALSE]
    R
}

rowNames <- function(x) {
    if (!is.datacggm(x)) stop(sQuote(x), " is not an object of class ", sQuote("datacggm"))
    Y <- getMatrix(x, "Y")
    X <- getMatrix(x, "X")
    if(is.null(x$X)) out <- list(Y = rownames(Y), X = NULL)
    else out <- list(Y = rownames(Y), X = rownames(X))
    out
}

`rowNames<-` <- function(x, value) {
    if (!is.datacggm(x)) stop(sQuote(x), " is not an object of class ", sQuote("datacggm"))
    nm <- names(value)
    if (length(value) > 2L) stop(sQuote("value"), " can not be a list with length greater than 2")
    if (any(is.null(nm))) stop(sQuote("value"), " is not a named list")
    if (any(!is.element(nm, c("X", "Y")))) stop("invalid 'rowname' given for 'datacggm' object")
    for (i in seq_len(length(nm))) {
        if (nm[i] == "Y") rownames(x$Info$R) <- rownames(x$Y) <- value$Y
        if (nm[i] == "X") rownames(x$X) <- value$X
    }
    x
}

colNames <- function(x) {
    if (!is.datacggm(x)) stop(sQuote(x), " is not an object of class ", sQuote("datacggm"))
    Y <- getMatrix(x, "Y")
    X <- getMatrix(x, "X")
    if(is.null(x$X)) out <- list(Y = colnames(Y), X = NULL)
    else out <- list(Y = colnames(Y), X = colnames(X))
    out
}

`colNames<-` <- function(x, value) {
    if (!is.datacggm(x)) stop(sQuote(x), " is not an object of class ", sQuote("datacggm"))
    nm <- names(value)
    if (length(value) > 2L) stop(sQuote("value"), " can not be a list with length greater than 2")
    if (any(is.null(nm))) stop(sQuote("value"), " is not a named list")
    if (any(!is.element(nm, c("X", "Y")))) stop("invalid 'colname' given for 'datacggm' object")
    for (i in seq_len(length(nm))) {
        if (nm[i] == "Y") {
            colnames(x$Y) <- value$Y
            colnames(x$Info$R) <- value$Y
            names(x$Info$lo) <- value$Y
            names(x$Info$up) <- value$Y
            names(x$Info$ym) <- value$Y
            names(x$Info$yv) <- value$Y
        }
        if (nm[i] == "X") colnames(x$X) <- value$X
    }
    x
}

upper <- function(x) {
    if (!is.datacggm(x)) stop(sQuote(x), " is not an object of class ", sQuote("datacggm"))
    big <- .Machine$double.xmax
    out <- x$Info$up
    out[out == big] <- +Inf
    out
}

lower <- function(x) {
    if (!is.datacggm(x)) stop(sQuote(x), " is not an object of class ", sQuote("datacggm"))
    small <- -.Machine$double.xmax
    out <- x$Info$lo
    out[out == small] <- -Inf
    out
}

ColMeans <- function(x) {
    if (!is.datacggm(x)) stop(sQuote("x"), " is not a ", sQuote("datacggm")," object")
#    out <- if(x$Info$q != 0L) list(Y = x$Info$ym, X = colMeans(x$X)) else list(Y = x$Info$ym, X = NULL)
    if (npred(x) != 0L) {
        X <- getMatrix(x, "X")
        id.numeric <- sapply(X, is.numeric)
        mx.num <- if (any(id.numeric)) {
            colMeans(X[, id.numeric, drop = FALSE])
        } else NULL
        mx.cat <- if (any(!id.numeric)) {
            # statistical mode
            apply(X[, !id.numeric, drop = FALSE], 2L, FUN = function(x) names(which.max(table(x))))
        } else NULL
        out <- list(Y = x$Info$ym, X.numeric = mx.num, X.categorical = mx.cat)
#        out <- list(Y = x$Info$ym, X = sapply(x$X, FUN = function(x) if (!is.numeric(x)) NA else mean(x)))
    } else {
        out <- list(Y = x$Info$ym, X.numeric = NULL, X.categorical = NULL)
    }
    out
}

ColVars <- function(x) {
    if (!is.datacggm(x)) stop(sQuote("x"), " is not a ", sQuote("datacggm")," object")
#    out <- if(x$Info$q != 0L) list(Y = x$Info$yv, X = apply(x$X, 2L, var)) else list(Y = x$Info$yv, X = NULL)
    if (npred(x) != 0L) {
        X <- getMatrix(x, "X")
        id.numeric <- sapply(X, is.numeric)
        vx.num <- if (any(id.numeric)) {
            apply(X[, id.numeric, drop = FALSE], 2L, var)
        } else NULL
        vx.cat <- if (any(!id.numeric)) {
            # Gini-Simpson Index (GSI)
            apply(X[, !id.numeric, drop = FALSE], 2L, FUN = function(x) 1 - sum(prop.table(table(x))^2))
        } else NULL
        out <- list(Y = x$Info$yv, X.numeric = vx.num, X.categorical = vx.cat)
    } else {
        out <- list(Y = x$Info$yv, X.numeric = NULL, X.categorical = NULL)
    }
    out
}

nobs.datacggm <- function(object, ...) {
    object$Info$n
}

nresp <- function(object, ...) {
    UseMethod("nresp")
}

nresp.datacggm <- function(object, ...) {
    object$Info$p
}

npred <- function(object, ...) {
    UseMethod("npred")
}

npred.datacggm <- function(object, ...) {
    object$Info$q
}

rcggm <- function(n, p, b0, X, B, Sigma, probl, probr, probna, ...) {
    ismissing <- c(n = missing(n), p = missing(p), b0 = missing(b0), X = missing(X), B = missing(B), Sigma = missing(Sigma),
                    probl = missing(probl), probr = missing(probr), probna = missing(probna))
    if (all(ismissing[c("n", "X")])) stop("Please, specify at least one of the following arguments: 'n' or 'X'")
    if (ismissing["X"] & !ismissing["B"]) stop(sQuote("X"), " is missing")
    if (!ismissing["X"] & ismissing["B"]) stop(sQuote("B"), " is missing")
    if (all(ismissing[c("p", "b0", "B", "Sigma")])) stop("Please, specify at least one of the following arguments: 'p', 'b0', or 'Sigma'")
    #########################
    # testing arguments
    if (!ismissing["n"]) {
        if (abs(as.integer(n) - n) > 0) stop(sQuote("n"), " is not a positive integer value")
        if (n <= 0) stop(sQuote("n"), " is not a strictly positive integer value")
    }
    if (ismissing["p"]) p <- NULL
    else {
        if (abs(as.integer(p) - p) > 0) stop(sQuote("p"), " is not a positive integer value")
        if (p <= 1L) stop(sQuote("p"), " is not an integer greater than ", sQuote("1"))
    }
    if (!ismissing["b0"]) {
        if (!is.vector(b0)) stop(sQuote("b0"), " is not a vector")
        if (is.null(p)) p <- length(b0)
        else {
            if (length(b0) != p) stop("The length of the vector ", sQuote("b0"), " is not equal to the argument ", sQuote("p"))
        }
    }
    if (ismissing["X"]) X <- NULL
    else {
        if (!is.numeric(X)) stop(sQuote("X"), " is not numeric")
        if (is.vector(X)) X <- matrix(X, nrow = 1L)
        if (!is.matrix(X)) X <- as.matrix(X)
        if (ismissing["n"]) n <- dim(X)[1L]
        else {
            if (dim(X)[1L] != n) stop("The number of rows of the matrix 'X' is not equal to the argument 'n'")
        }
    }
    if (ismissing["B"]) B <- NULL
    else {
        if (!is.numeric(B)) stop(sQuote("B"), " is not numeric")
        if (is.vector(B)) B <- matrix(B, nrow = 1L)
        if (!is.matrix(B)) stop(sQuote("B"), " is not a matrix")
        if (dim(B)[1L] != dim(X)[2L]) stop("The number of rows of the matrix 'B' is not equal to the number of columns of the matrix 'X'")
        if (is.null(p)) p <- dim(B)[2L]
        else {
            if (dim(B)[2L] != p) stop("The number of columns of the matrix 'B' is not equal to the length of the vector 'b0'")
        }
    }
    if (ismissing["Sigma"]) Sigma <- diag(p)
    else {
        if (!is.matrix(Sigma)) stop(sQuote("Sigma"), " is not a matrix")
        if (!isSymmetric(Sigma)) stop(sQuote("Sigma"), " is not a symmetric matrix")
        if (dim(Sigma)[1L] == 1L) stop(sQuote("Sigma"), " is a matrix with dimension 1 x 1")
        if (is.null(p)) p <- dim(Sigma)[1L]
        if (!ismissing["B"]) {
            if (dim(B)[2L] != dim(Sigma)[1L]) stop("the number of columns of the matrix 'B' is not coherent with the dimensions of the matrix 'Sigma'")
        }
    }
    if (ismissing["b0"]) b0 <- rep(0, p)
    if (length(b0) != dim(Sigma)[1L]) stop("the length of the vector 'b0' is not coherent with the dimensions of the matrix 'Sigma'")
    if (!ismissing["probl"]) {
        if (!is.vector(probl)) stop("'probl' is not a vector")
        if (length(probl) == 1L) probl <- rep(probl, p)
        if (!all(0 <= probl & probl <= 1)) stop("some element in 'probl' does not belong to the interval [0, 1]")
    } else probl <- rep(0, p)
    if (!ismissing["probr"]) {
        if (!is.vector(probr)) stop("'probr' is not a vector")
        if (length(probr) == 1L) probr <- rep(probr, p)
        if (!all(0 <= probr & probr <= 1)) stop("some element in 'probr' does not belong to the interval [0, 1]")
    } else probr <- rep(0, p)
    if (any(probl + probr >= 1)) stop("'probl' plus 'probr' can not be greater than or equal to 1")
    if (!ismissing["probna"]) {
        if (!is.vector(probna)) stop("'probna' is not a vector")
        if (length(probna) != 1L) stop("'probna' is not a vector of length ", sQuote("1"))
        if (probna < 0 | probna >= 1) stop("'probna' does not belong to the interval [0, 1)")
    } else probna <- 0
    if (any(probl + probr + probna >= 1)) stop("'probl' plus 'probr' plus 'probna' can not be greater than or equal to 1")
    ########################
    # computing matrix mu
    mu <- matrix(b0, nrow = n, ncol = p, byrow = TRUE)
    if (!ismissing["X"]) mu <- mu + X %*% B
    ########################
    # computing matrix Y
    Y <- mu + mvrnorm(n = n, mu = rep(0, p), Sigma = Sigma, ...)
    if(n == 1L) Y <- matrix(Y, nrow = 1L)
    ########################
    # Computing lo and up
    lo <- double(p)
    up <- double(p)
    lo.find <- function(lo, mu, sd, probl) mean(pnorm(lo, mean = mu, sd = sd, lower.tail = TRUE)) - probl
    up.find <- function(up, mu, sd, probr) mean(pnorm(up, mean = mu, sd = sd, lower.tail = FALSE)) - probr
    for (m in seq_len(p)) {
        mu.m <- mu[, m]
        sd.m <- sqrt(Sigma[m, m])
        int.m <- c(min(mu.m) - 10 * sd.m, max(mu.m) + 10 * sd.m)
        if (probl[m] == 0)
            lo[m] <- -Inf
        else {
            lo[m] <- uniroot(lo.find, interval = int.m, mu = mu.m, sd = sd.m, probl = probl[m])$root
            Y[Y[, m] <= lo[m], m] <- lo[m]
        }
        if (probr[m] == 0)
            up[m] <- +Inf
        else {
            up[m] <- uniroot(up.find, interval = int.m, mu = mu.m, sd = sd.m, probr = probr[m])$root
            Y[Y[, m] >= up[m], m] <- up[m]
        }
        if (probna > 0) {
            id.obs <- lo[m] < Y[, m] & Y[, m] < up[m]
            NOBS <- sum(id.obs)
            if (NOBS > 0) {
                id.na <- rbinom(NOBS, size = 1L, prob = probna / (1 - probl[m] - probr[m]))
                Y[id.obs, m][id.na == 1L] <- NA
            }
        }
    }
    out <- datacggm(Y = Y, lo = lo, up = up, X = X)
    out
}

qqcnorm <- function(x, which, max.plot = 1L, save.plot = FALSE, grdev = pdf, grdev.arg,
                    main = "Censored Normal Q-Q Plot", xlab = "Theoretical Quantiles", ylab = "Sample Quantiles",
                    plot.it = TRUE, plot.pch = c(2L, 20L), plot.col = c(2L, 1L), plot.cex = c(2L, 1L), abline = FALSE,
                    line.col = "gray50", line.lwd = 1L, line.lty = 2L, ...) {
    # testing 'which'
    p <- nresp(x)
    if (missing(which)) which <- seq_len(p)
    else {
        if (!is.vector(which)) stop (sQuote("which"), "is not a vector")
        if (any(abs(as.integer(which) - which) > 0)) stop(sQuote("which"), " is not an object of type ", dQuote("integer"))
        if (min(which) <= 0) stop("some entry in ", sQuote("which"), " is not a positive integer")
        if (any(which > p)) stop("some entry in ", sQuote("which"), " is larger than ", sQuote(p))
    }
    # testing 'max.plot'
    if (!is.vector(max.plot)) stop(sQuote("which"), "is not a vector")
    if (length(max.plot) != 1L) stop(sQuote("max.plot"), "is not a vector of length ", sQuote("1"))
    if (abs(as.integer(max.plot) - max.plot) > 0) stop(sQuote("max.plot"), " is not an object of type ", dQuote("integer"))
    if (max.plot <= 0L) stop(sQuote("max.plot"), " is not a positive integer")
    # testing 'save.plot'
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
    # testing 'grdev.arg'
    if (missing(grdev.arg)) {
        grdev.arg <- NULL
    } else {
        if (!is.list(grdev.arg)) stop(sQuote("grdev.arg"), " is not an object of type ", sQuote("list"))
    }
    n <- nobs(x)
    R <- event(x)
    Y <- getMatrix(x, "Y")
    lo <- lower(x)
    up <- upper(x)
    ymu <- ColMeans(x)$Y
    ysd <- sqrt(ColVars(x)$Y)
    nwhich <- length(which)
    ngrp <- ceiling(nwhich / max.plot)
    out <- vector(mode = "list", length = nwhich)
    names(out) <- colNames(x)$Y[which]
    op <- par(no.readonly = TRUE)
    op$new <- FALSE
    on.exit(par(op), add = TRUE)
    if(save.plot) {
        cat("Exporting plots\n")
        pb <- txtProgressBar(min = 0L, max = nwhich, style = 3L)
        on.exit(close(pb), add = TRUE, after = TRUE)
    }
    hh <- 0
    if (ngrp > 1L) devAskNewPage(TRUE)
    grp <- rep(seq_len(ngrp), each = max.plot, length = nwhich)
    if (save.plot)
        file.names <- paste0("qqcnorm_group_", formatC(seq_len(ngrp), width = nchar(ngrp), format = "d", flag = "0"), ".", grdev.name)
    for(h in seq_len(ngrp)) {
        if (save.plot) {
            if (is.null(grdev.arg)) grdev(file = file.names[h])
            else do.call(function(...) grdev(file = file.names[h], ...), grdev.arg)
        }
        if(ngrp > 1L) dev.hold()
        nplot <- length(which[grp == h])
        if(nwhich > 1) par(mfrow = rev(n2mfrow(nplot)))
        for (i in which[grp == h]) {
            y <- Y[R[, i] != +9L, i]
            status <- R[R[, i] != +9L, i]
            ny <- length(y)
            pp <- ppoints(ny)
            y.q <- qnorm(pp, mean = ymu[i], sd = ysd[i])[order(order(y))]
            y.q[ y.q <= lo[i] ] <- lo[i]
            status[y.q <= lo[i] ] <- -1L
            y.q[ y.q >= up[i] ] <- up[i]
            status[ y.q >= up[i] ] <- +1L
            if (plot.it) {
                pch <- ifelse (status != 0, plot.pch[1L], plot.pch[2L])
                col <- ifelse (status != 0, plot.col[1L], plot.col[2L])
                cex <- ifelse (status != 0, plot.cex[1L], plot.cex[2L])
                plot(y.q, y, main = main, xlab = xlab, ylab = ylab, pch = pch,
                     col = col, cex = cex, ...)
                if (abline)
                    abline(0, 1, col = line.col, lwd = line.lwd, lty = line.lty)
            }
            
            hh <- hh + 1
            out[[hh]] <- cbind(x = y.q, y = y)
            if (save.plot)  setTxtProgressBar(pb, hh)
        }
        if(ngrp > 1L) dev.flush()
        if(save.plot) dev.off()
    }
    invisible(out)
}
