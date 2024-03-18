print.cggm <- function (x, digits = 3L, ...){
    p <- nresp(x$Z)
    q <- npred(x$Z)
    dots <- list(...)
    if (is.null(dots$print.gap)) dots$print.gap <- 2L
    if (is.null(dots$quote)) dots$quote <- FALSE
    if (is.null(dots$row.names)) dots$row.names <- FALSE
    df.B <- x$dfB[p + 1L, , ] +  p
    df.Tht <- drop(x$dfTht) + p
    df <- df.B + df.Tht
    df.max <- p * (p + 1) / 2 + (q + 1L) * p
    df.per <- formatC(round(df / df.max * 100, digits = digits), format = 'f', digits = digits)
    df.per <- paste("(", df.per, "%)", sep = "")
    ncomp <- drop((x$InfoStructure$ncomp))
    tbl <- data.frame(df, df.per, ncomp)
    names(tbl) <- c("df", "", "N. Comp.")
    cat("\nCall:  ", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
    do.call(function(...) print.data.frame(tbl, digits = digits, ...), dots)
    cat("\n")
}

predict.cggm <- function(object, X.new, ...) {
    q <- npred(object$Z)
    if (q == 0L) stop(sQuote("X.new"), " can not be used because no predictors are available in ", sQuote("object"))
    if (missing(X.new)) stop(sQuote("X.new"), " is missing")
    else {
        if (is.vector(X.new)) X.new <- as.matrix(X.new)
        if (!is.matrix(X.new)) stop(sQuote("X.new"), " is not a matrix")
        if (dim(X.new)[2L] != q) stop("wrong dimension in ", sQuote("X.new"),". The number of columns in not equal to ", sQuote(q))
    }
    B <- coef(object, type = "B", drop = TRUE)
    mu <- cbind(1, X.new) %*% B
    mu
}

plot.cggm <- function(x, type, weighted = FALSE, simplify = TRUE, ...) {
    out <- to_graph(x, weighted = weighted, simplify = simplify)
    if(missing(type)) type <- ifelse(is.null(out$Gxy), "Gyy", "both")
    dots <- list(...)
    do.call(function(...) plot(out, type = type, ...), dots)
    invisible(NULL)
}


