is.datacggm <- function(x) inherits(x, 'datacggm')

dim.datacggm <- function(x) {
    out <- list(Y = dim(getMatrix(x, "Y")), X = dim(getMatrix(x, "X")))
    out
}

dimnames.datacggm <- function(x) {
    Y <- getMatrix(x, "Y")
    out <- list(Y = dimnames(Y), X = NULL)
    if (npred(x) != 0) {
        X <- getMatrix(x, "X")
        out$X <- dimnames(X)
    }
    out
}

`dimnames<-.datacggm` <- function(x, value) {
    if (!is.list(value)) stop(sQuote("value"), " is not a list")
    if (length(value) > 2L) stop(sQuote("value"), " can not be a list with length greater than 2")
    nm <- names(value)
    if (any(is.null(nm))) stop(sQuote("value"), " is not a named list")
    if (any(!is.element(nm, c("X", "Y")))) stop("Wrong names in ", sQuote("value"))
    for (i in seq_len(length(nm))) {
        if (nm[i] == "Y") {
            dimnames(x$Y) <- value$Y
            dimnames(x$Info$R) <- value$Y
            names(x$Info$lo) <- colnames(x$Y)
            names(x$Info$up) <- colnames(x$Y)
            names(x$Info$ym) <- colnames(x$Y)
            names(x$Info$yv) <- colnames(x$Y)
        } else {
            rownames(x$X) <- value$X[[1L]]
            colnames(x$X) <- value$X[[2L]]
        }
    }
    x
}

Math.datacggm <- function(...)  stop("Invalid operation on a 'datacggm' object")

Ops.datacggm  <- function(...)  stop("Invalid operation on a 'datacggm' object")

as.character.datacggm <- function (x, digits, ...) {
    Y <- format(getMatrix(x, "Y"), digits = digits, ...)
    marker <- event(x)
    Y[marker == 1] <- paste0(Y[marker == 1], "+")
    Y[marker == -1] <- paste0(Y[marker == -1], "-")
    Y
}

print.datacggm <- function(x, digits = 3L, n = 10L, width = getOption("width"), ...) {
    # testing arguments
    if (!is.infinite(n)) {
        if (abs(as.integer(n) - n) > 0) stop(sQuote("n"), " is not a positive integer value")
        if (n <= 0) stop(sQuote("n"), " is not a strictly positive integer value")
    }
    if (!is.infinite(width)) {
        if (abs(as.integer(width) - width) > 0) stop(sQuote("width"), " is not a positive integer value")
        if (width <= 0) stop(sQuote("width"), " is not a strictly positive integer value")
    }
    out <- as.character.datacggm(x, digits = digits, ...)
    out <- fit_matrix(out, n = n, width = width)
    cat("Printing", sQuote("datacggm"), "object\n\n")
    cat("Y:", nobs(x), "x", nresp(x), "matrix\n\n")
    print(out, quote = FALSE)
    dn <- nobs(x) - n
    dp <- nresp(x) - dim(out)[2L]
    if (dn > 0 | dp > 0) {
        msg <- paste0("# with ", ifelse(dn > 0, paste0(dn, " more ", ifelse(dn > 1, "rows", "row")), ""),
                    ifelse(dn > 0 & dp > 0, ", and ", ""),
                    ifelse(dp > 0, paste0(dp, " more ", ifelse(dp > 1, "variables", "variable")), ""))
        cat(msg)
    }
    if (npred(x) == 0L) cat("\n\nX: NULL")
    else {
        out <- format(getMatrix(x, "X"), digits = digits, ...)
        out <- fit_matrix(out, n = n, width = width)
        cat("\n\nX:", nobs(x), "x", npred(x), "matrix\n\n")
        print(out, quote = FALSE)
        dq <- npred(x) - dim(out)[2L]
        if (dn > 0 | dq > 0) {
            msg <- paste0("# with ", ifelse(dn > 0, paste0(dn, " more ", ifelse(dn > 1, "rows", "row")), ""),
                        ifelse(dn > 0 & dq > 0, ", and ", ""),
                        ifelse(dq > 0, paste0(dq, " more ", ifelse(dq > 1, "variables", "variable")), ""))
            cat(msg)
        }
    }
    cat("\n")
}

summary.datacggm <- function(object, n, quantile.type = 7L, digits = 3L, quote = FALSE, ...) {
    if (missing(n)) n <- Inf
    Y <- getMatrix(object, "Y")
    X <- getMatrix(object, "X")
    R <- event(object)
    ynm <- colNames(object)$Y
    p <- nresp(object)
    q <- npred(object)
    lo <- lower(object)
    up <- upper(object)
    out <- list(Y = NULL, X.numeric = NULL, X.categorical = NULL)
    lcs <- paste0(format(apply(R == -1, 2L, function(x) mean(x) * 100), digits = digits), "%")
    rcs <- paste0(format(apply(R == +1, 2L, function(x) mean(x) * 100), digits = digits), "%")
    nas <- paste0(format(apply(R == +9, 2L, function(x) mean(x) * 100), digits = digits), "%")
    qq <- matrix(0, nrow = p, ncol = 6L)
    for (m in seq_len(p)) {
        obs <- Y[R[, m] == 0, m]
        qq[m, -4L] <- stats::quantile(obs, names = FALSE, type = quantile.type)
        qq[m, 4L] <- mean(obs)
    }
    tbl <- apply(cbind(lo, qq, up), 2L, format, digits = digits)
    tbl <- cbind(tbl, nas, lcs, rcs)
    colnames(tbl) <- c("Lower", "Min.", "1st Qu.", "Median", "Mean", "3rd Qu.", "Max.", "Upper", "NA%", "LC%", "RC%")
    rownames(tbl) <- ynm
    attr(tbl, "class") <- "table"
    out$Y <- tbl
    if (q > 0L) {
        id.numeric <- sapply(X, is.numeric)
        if (any(id.numeric)) {
            tbl <- apply(X[, id.numeric, drop = FALSE], 2L, function(obs) {
                            qq <- stats::quantile(obs, names = FALSE, type = quantile.type)
                            c(qq[1L:3L], mean(obs), qq[4L:5L])})
            if (q == 1L) tbl <- format(t(tbl), digits = digits)
            else tbl <- apply(t(tbl), 2L, format, digits = digits)
            colnames(tbl) <- c("Min.", "1st Qu.", "Median", "Mean", "3rd Qu.", "Max.")
            attr(tbl, "class") <- "table"
            out$X.numeric <- tbl
        }
        if (any(!id.numeric)) {
            tbl <- lapply(X[, !id.numeric, drop = FALSE], FUN = function(x) {
                            Freq.x <- table(x)
                            Perc.x <- prop.table(Freq.x) * 100
                            out <- cbind(Freq = format(Freq.x, digits = digits), Perc = paste0(format(Perc.x, digits = digits), "%"))
                            attr(out, "class") <- "table"
                            out
                            })
#            tbl <- apply(X[, !id.numeric, drop = FALSE], 2L, FUN = function(x) {
#                            Freq.x <- table(x)
#                            Perc.x <- prop.table(Freq.x) * 100
#                            out <- cbind(Freq = format(Freq.x, digits = digits), Perc = paste0(format(Perc.x, digits = digits), "%"))
#                            attr(out, "class") <- "table"
#                            out
#                            }, simplify = FALSE)
            out$X.categorical <- tbl
        }
    }
    dvar <- dim(out$Y)[1L] - n
    cat("Y:\n")
    print(head(out$Y, n), quote = quote, ...)
    if (dvar > 0)
        cat(paste0("# with ", ifelse(dvar > 0, paste0(dvar, " more ", ifelse(dvar > 1, "variables", "variable")))))
    if (q == 0L) cat("\n\nX-numeric: NULL\n\nX-categorical: NULL\n\n")
    else {
        if (any(id.numeric)) {
            dvar <- dim(out$X.numeric)[1L] - n
            cat("\n\nX-numeric:\n")
            print(head(out$X.numeric, n), quote = quote, ...)
            if (dvar > 0)
                cat(paste0("# with ", ifelse(dvar > 0, paste0(dvar, " more ", ifelse(dvar > 1, "variables", "variable")))))
        }
        if (any(!id.numeric)) {
            n <- ifelse(n >= length(out$X.categorical), length(out$X.categorical), n)
            dvar <- length(out$X.categorical) - n
            #dvar <- ifelse(dvar <= 0, length(out$X.categorical), dvar)
            cat("\nX-categorical:\n")
            print.listof(out$X.categorical[seq_len(n)], quote = quote, ...)
            if (dvar > 0)
                cat(paste0("# with ", ifelse(dvar > 0, paste0(dvar, " more ", ifelse(dvar > 1, "variables", "variable")))))
        }
    }
    cat("\n")
    invisible(out)
}

hist.datacggm <- function(x, breaks = "Sturges", include.lowest = TRUE, right = TRUE, nclass = NULL,
                          which, max.hist = 1L, save.hist = FALSE, grdev = pdf, grdev.arg,
                          polygon.col = adjustcolor("grey", alpha.f = 0.25), polygon.border = NA, segments.lwd = 4L, segments.lty = 2L,
                          segments.col = "gray40", points.pch = c(4L, 1L), points.cex = 1.8, points.col = rep("black", 2L),
                          legend = TRUE, ...) {
    dots <- list(...)
    # testing 'which'
    p <- nresp(x)
    if (missing(which)) which <- seq_len(p)
    else {
        if (!is.vector(which)) stop (sQuote("which"), "is not a vector")
        if (any(abs(as.integer(which) - which) > 0)) stop(sQuote("which"), " is not an object of type ", dQuote("integer"))
        if (min(which) <= 0) stop("some entry in ", sQuote("which"), " is not a positive integer")
        if (any(which > p)) stop("some entry in ", sQuote("which"), " is larger than ", sQuote(p))
    }
    # testing 'max.hist'
    if (!is.vector(max.hist)) stop(sQuote("which"), "is not a vector")
    if (length(max.hist) != 1L) stop(sQuote("max.hist"), "is not a vector of length ", sQuote("1"))
    if (abs(as.integer(max.hist) - max.hist) > 0) stop(sQuote("max.hist"), " is not an object of type ", dQuote("integer"))
    if (max.hist <= 0L) stop(sQuote("max.hist"), " is not a positive integer")
    # testing 'save.hist'
    if (length(save.hist) != 1L) stop(sQuote(save.hist), " is not an object of length ", sQuote("1"))
    if (!inherits(save.hist, c("logical", "character"))) stop(sQuote(save.hist), " is not an object of type ", sQuote("logical"), "or ", , sQuote("character"))
    if(inherits(save.hist, "character")) {
        oldPath <- getwd()
        newPath <- save.hist
        setwd(newPath)
        on.exit(setwd(oldPath), add = TRUE, after = TRUE)
        save.hist <- TRUE
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
    if (missing(grdev.arg)) grdev.arg <- NULL
    else {
        if (!is.list(grdev.arg)) stop(sQuote("grdev.arg"), " is not an object of type ", sQuote("list"))
    }
    R <- event(x)
    Y <- getMatrix(x, "Y")
    lo <- lower(x)
    up <- upper(x)
    ynames <- colNames(x)$Y
    ymu <- ColMeans(x)$Y
    ysd <- sqrt(ColVars(x)$Y)
    nwhich <- length(which)
    ngrp <- ceiling(nwhich / max.hist)
    op <- par(no.readonly = TRUE)
    op$new <- FALSE
    on.exit(par(op), add = TRUE)
    if(save.hist) {
        file.names <- paste0("histogram_group_", formatC(seq_len(ngrp), width = nchar(ngrp),
                             format = "d", flag = "0"), ".", grdev.name)
        cat("Exporting plots\n")
        pb <- txtProgressBar(min = 0L, max = nwhich, style = 3L)
        on.exit(close(pb), add = TRUE, after = TRUE)
        hh <- 0
    }
    if (ngrp > 1L) devAskNewPage(TRUE)
    grp <- rep(seq_len(ngrp), each = max.hist, length = nwhich)
    if(legend) par(xpd = TRUE, mar = par()$mar + c(0, 0, 0, 7.5))
    for(h in seq_len(ngrp)) {
        if(save.hist) {
            if (is.null(grdev.arg)) grdev(file = file.names[h])
            else do.call(function(...) grdev(file = file.names[h], ...), grdev.arg)
            if(legend) par(xpd = TRUE, mar = par()$mar + c(0, 0, 0, 7.5))
        }
        if(ngrp > 1L) dev.hold()
        nhist <- length(which[grp == h])
        if(nwhich > 1) par(mfrow = rev(n2mfrow(nhist)))
        for (i in which[grp == h]) {
            yo <- Y[R[, i] == 0L, i]
            n <- sum(R[, i] != 9L)
            flc <- mean(R[, i] == -1L)
            frc <- mean(R[, i] == +1L)
            plc <- pnorm(q = (lo[i] - ymu[i]) / ysd[i], lower.tail = TRUE)
            prc <- pnorm(q = (up[i] - ymu[i]) / ysd[i], lower.tail = FALSE)
            out.hist <- hist(x = yo, breaks = breaks, include.lowest = include.lowest, right = right, plot = FALSE)
            brks <- out.hist$breaks
            nbreaks <- length(brks)
            xlim <- double(2L)
            if (is.finite(lo[i])) {
                brks[1L] <- max(brks[1L], lo[i])
                xlim[1L] <- min(brks[1L], lo[i])
                out.hist$mids[1L] <- (brks[2L] + brks[1L]) / 2
            } else xlim[1L] <- brks[1L]
            if (is.finite(up[i])) {
                brks[nbreaks] <- min(brks[nbreaks], up[i])
                xlim[2L] <- max(brks[nbreaks], up[i])
                out.hist$mids[nbreaks] <- (brks[nbreaks] + brks[nbreaks - 1L]) / 2
            } else xlim[2L] <- brks[nbreaks]
            dbreaks <- diff(brks)
            out.hist$breaks <- brks
            out.hist$density <- out.hist$counts / (n * dbreaks)
            out.hist$equidist <- diff(range(dbreaks)) < 1e-07 * mean(dbreaks)
            out.hist$xname <- ynames[i]
            seqx <- seq(xlim[1L], xlim[2L], length.out = 100L)
            seqy <- dnorm(seqx, mean = ymu[i], sd = ysd[i])
            seqx <- c(seqx, rev(seqx))
            seqy <- c(seqy, rep(0, 100L))
            ylim <- c(0, max(out.hist$density, seqy, flc, frc, plc, prc))
            do.call(function(...) plot(out.hist, freq = FALSE, xlim = xlim, ylim = ylim, ...), dots)
            polygon(seqx, seqy, col = polygon.col, border = polygon.border)
            if (flc != 0) {
                segments(x0 = lo[i], y0 = 0, x1 = lo[i], y1 = plc,
                         lwd = segments.lwd, lty = segments.lty, col = segments.col)
                points(x = rep(lo[i], 2L), y = c(flc, plc),
                       pch = points.pch, cex = points.cex, col = points.col)
            }
            if (frc != 0) {
                segments(x0 = up[i], y0 = 0, x1 = up[i], y1 = prc,
                         lwd = segments.lwd, lty = segments.lty, col = segments.col)
                points(x = rep(up[i], 2L), y = c(frc, prc),
                       pch = points.pch, cex = points.cex, col = points.col)
                
            }
            if ((flc != 0 | frc != 0) & legend)
                legend(xlim[2L], ylim[2L], pch = c(4L, 1L), bty = "n", pt.cex = 1.8,
                       c("Proportion of\ncensored values\n", "Probability of a\ncensored value"))
            if(save.hist) {
                hh <- hh + 1
                setTxtProgressBar(pb, hh)
            }
        }
        if(ngrp > 1L) dev.flush()
        if(save.hist) dev.off()
    }
}
