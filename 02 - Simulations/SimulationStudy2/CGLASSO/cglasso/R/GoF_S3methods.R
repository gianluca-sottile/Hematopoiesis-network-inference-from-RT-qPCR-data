print.QFun <- function (x, digits = 3L, ...){
    nrho <- x$nrho
    rho <- x$rho
    nlambda <- x$nlambda
    lambda <- x$lambda
    if (nrho == 1L | nlambda == 1L) {
        value <- drop(x$value)
        df <- drop(x$df)
    }
    else {
        value <- as.vector(t(x$value))
        df <- as.vector(t(x$df))
    }
    if (nrho > 1L & nlambda > 1L) {
        lambda <- rep(lambda, each = nrho)
        rho <- rep(rho, times = nlambda)
    }
    dots <- list(...)
    if (is.null(dots$print.gap)) dots$print.gap <- 2L
    if (is.null(dots$quote)) dots$quote <- FALSE
    if (is.null(dots$row.names)) dots$row.names <- FALSE
    tbl <- data.frame(lambda = lambda, rho = rho, df = df, value = value)
    names(tbl)[4L] <- "Q-Values"
    cat("\nQ-values of the fitted", sQuote(x$model[1L]), ifelse(nrho == 1L & nlambda == 1L, "model", "models"))
    cat("\n\nDetails:\n")
    if (nrho == 1L | nlambda == 1L) {
        if (x$q == 0L) tbl <- tbl[, -1L, drop = FALSE]
        do.call(function(...) print.data.frame(tbl, digits = digits, ...), dots)
        cat("\n")
    }
    else {
        if (nlambda <= nrho) f <- rep(seq_len(nlambda), each = nrho)
        else f <- rep(seq_len(nrho), times = nlambda)
        tbl.list <- split(tbl, f = f, drop = FALSE)
        do.call(function(...) print.listof(tbl.list, digits = digits, ...), dots)
    }
    invisible(tbl)
}

print.GoF <- function (x, digits = 3L, ...){
    type <- x$type
    nrho <- x$nrho
    rho <- x$rho
    nlambda <- x$nlambda
    lambda <- x$lambda
    if (nrho == 1L | nlambda == 1L) {
        val <- drop(x$value_gof)
        df <- drop(x$df)
    }
    else {
        val <- as.vector(t(x$value_gof))
        df <- as.vector(t(x$df))
    }
    if (nrho > 1L & nlambda > 1L) {
        lambda <- rep(lambda, each = nrho)
        rho <- rep(rho, times = nlambda)
    }
    tbl <- data.frame(lambda = lambda, rho = rho, df = df, val = val)
    names(tbl)[4L] <- type
    cat("\nSequence of", sQuote(type), "values of the fitted", sQuote(x$model[1L]), ifelse(nrho == 1L & nlambda == 1L, "model", "models"))
    cat("\n\nDetails:\n")
    dots <- list(...)
    if (is.null(dots$print.gap)) dots$print.gap <- 2L
    if (is.null(dots$quote)) dots$quote <- FALSE
    if (is.null(dots$row.names)) dots$row.names <- FALSE
    if (nrho == 1L | nlambda == 1L) {
        if (x$q == 0L) tbl <- tbl[, -1L, drop = FALSE]
        do.call(function(...) print.data.frame(tbl, digits = digits, ...), dots)
        cat("\n")
    }
    else {
        if (nlambda <= nrho) f <- rep(seq_len(nlambda), each = nrho)
        else f <- rep(seq_len(nrho), times = nlambda)
        tbl.list <- split(tbl, f = f, drop = FALSE)
        do.call(function(...) print.listof(tbl.list, digits = digits, ...), dots)
    }
    invisible(tbl)
}

plot.GoF <- function(x, add.line = TRUE, arg.line = list(lty = 2L, lwd = 2L, col = "red"), add.text = FALSE, arg.text = list(side = 3L), arg.points = list(pch = 2L), ...) {
    if (!is.vector(add.line)) stop(sQuote("add.line"), " is not a vector")
    if (length(add.line) != 1L) stop(sQuote("add.line"), " is not a vector of length ", sQuote("1"))
    if (!is.logical(add.line)) stop(sQuote("add.line"), " is not an object of type ", sQuote("logical"))
    if (!is.list(arg.line)) stop(sQuote("arg.line"), " is not a list")
    if (!is.vector(add.text)) stop(sQuote("add.text"), " is not a vector")
    if (length(add.text) != 1L) stop(sQuote("add.text"), " is not a vector of length ", sQuote("1"))
    if (!is.logical(add.text)) stop(sQuote("add.text"), " is not an object of type ", sQuote("logical"))
    if (!is.list(arg.text)) stop(sQuote("arg.text"), " is not a list")
    nrho <- x$nrho
    nlambda <- x$nlambda
    if (nrho == 1L & nlambda == 1L) stop("plot method is not available because ", sQuote("nlambda = 1"), " and ", sQuote("nrho = 1"))
    dots <- list(...)
    if (is.null(arg.line$lty)) arg.line$lty <- 2L
    if (is.null(arg.line$lwd)) arg.line$lwd <- 2L
    if (is.null(arg.line$col)) arg.line$col <- "red"
    if (is.null(arg.text$side)) arg.text$side <- 3L
    if (is.null(arg.points$pch)) arg.points$pch <- 2L
    if (nrho == 1L | nlambda == 1L) {
        val <- drop(x$value_gof)
        pen <- if(nlambda == 1L) x$rho else x$lambda
        df <- drop(x$df)
        minval <- which.min(val)
        if (is.null(dots$main)) dots$main <- "Tuning Parameter Selection"
        if (is.null(dots$sub)) {
            dots$sub <- switch(x$type,
                "AIC" = "Akaike Information Criterion",
                "BIC" = "Bayesian Information Criterion",
                "GoF" = "Measure of Goodness of Fit",
                "eBIC_FD" = "extended Bayesian Information Criterion",
                "eBIC_CC" = "extended Bayesian Information Criterion")
        }
        if (is.null(dots$xlab)) dots$xlab <- ifelse(nlambda == 1L, expression(rho), expression(lambda))
        if (is.null(dots$ylab)) dots$ylab <- "Values"
        if (is.null(dots$type)) dots$type <- "b"
        do.call(function(...) plot(x = pen, y = val, ...), dots)
        if (add.line) {
            do.call(function(...) abline(v = pen[minval], ...), arg.line)
            if (add.text) {
                if (is.null(arg.text$text)) arg.text$text <- paste0("df = ", df[minval])
                do.call(function(...) mtext(at = pen[minval], ...), arg.text)
            }
        }
    } else {
        rho <- rev(x$rho)
        lambda <- rev(x$lambda)
        val <- x$value_gof[nlambda:1, nrho:1]
        minval <- min(val)
        rc.cord <- drop(which(val == minval, arr.ind = TRUE))
        if (is.matrix(rc.cord)) rc.cord <- rc.cord[dim(rc.cord)[1L], ]
        if (is.null(dots$main)) dots$main <- "Tuning Parameters Selection"
        if (is.null(dots$xlab)) dots$xlab <- expression(lambda)
        if (is.null(dots$ylab)) dots$ylab <- expression(rho)
        if (is.null(dots$sub)) {
            dots$sub <- switch(x$type,
                "AIC" = "Akaike Information Criterion",
                "BIC" = "Bayesian Information Criterion",
                "GoF" = "Measure of Goodness of Fit",
                "eBIC_FD" = "extended Bayesian Information Criterion",
                "eBIC_CC" = "extended Bayesian Information Criterion")
        }
        do.call(function(...) contour(x = lambda, y = rho, z = val, ...), dots)
        if (add.line) {
            do.call(function(...) abline(h = rho[rc.cord[2L]], ...), arg.line)
            do.call(function(...) abline(v = lambda[rc.cord[1L]], ...), arg.line)
        }
        do.call(function(...) points(x = lambda[rc.cord[1L]], y = rho[rc.cord[2L]], ...), arg.points)
    }
}
