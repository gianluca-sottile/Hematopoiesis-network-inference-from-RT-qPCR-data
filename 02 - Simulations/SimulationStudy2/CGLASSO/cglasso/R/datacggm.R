datacggm <- function(Y, lo = -Inf, up = +Inf, X = NULL, control = list(maxit = 1.0E+4, thr = 1.0E-4)) {
    big <- .Machine$double.xmax
    thr <- big / 2
    # checking 'Y'
    if (missing(Y)) stop(sQuote("Y"), " is missing")
    if (!is.matrix(Y)) stop(sQuote("Y"), " is not a matrix")
    if (any(!is.finite(Y[!is.na(Y)]))) stop("some element in 'Y' is not finite")
    n <- dim(Y)[1L]
    p <- dim(Y)[2L]
    if (is.null(colnames(Y))) {
        vnames <- paste0("Y", 1:p, sep = "")
        colnames(Y) <- vnames
    } else vnames <- colnames(Y)
    if (is.null(rownames(Y))) rownames(Y) <- seq_len(n)
    # checking 'lo' and 'up'
    if (missing(lo)) lo <- rep(-big, p)
    else {
        if (!is.vector(lo)) stop(sQuote("lo"), " is not a vector")
        if (length(lo) == 1) lo <- rep(lo, p)
        id <- -big < lo & lo <= -thr
        if (any(id)) {
            lo[id] <- -big
            message("message: some entries in 'lo' are below the tolerance. These values are treated as -Inf")
        }
        id <- lo == -Inf
        if (any(id)) lo[id] <- -big
    }
    names(lo) <- vnames
    if (missing(up)) up <- rep(big, p)
    else {
        if (!is.vector(up)) stop(sQuote("up"), " is not a vector")
        if (length(up) == 1) up <- rep(up, p)
        id <- thr <= up & up < big
        if (any(id)) {
            up[id] <- big
            message("message: some entries in 'up' are over the tolerance. These values are treated as +Inf")
        }
        id <- up == Inf
        if(any(id)) up[id] <- big
    }
    names(up) <- vnames
    if (!all(lo < up)) stop(sQuote("lo"), " is not less than ", sQuote("up"))
    if (!is.null(X)) {
#        if (!is.numeric(X)) stop(sQuote("X"), " is not numeric")
        if (any(is.na(X))) stop("Missing values in ", sQuote("X"), " are not allowed")
        if (is.vector(X)) X <- matrix(X, ncol = 1)
        if (length(dim(X)) != 2L) stop(sQuote("X"), "is not a matrix like object")
        if (dim(Y)[1L] != dim(X)[1L]) stop(sQuote("X"), " and ", sQuote("Y"), " have a different number of rows")
        if (is.null(colnames(X))) colnames(X) <- paste0("X", seq_len(dim(X)[2L]), sep = "")
        if (is.null(rownames(X))) rownames(X) <- seq_len(n)
        if (!is.data.frame(X)) X <- as.data.frame(X)
        q <- dim(X)[2L]
    } else q <- 0L
    # Testing control
    if (!is.list(control)) stop(sQuote("control"), "is not an object of type ", sQuote("list"))
    if (is.null(names(control))) stop(sQuote("control"), "is not a named list")
    names(control) <- sapply(names(control), match.arg, choice = c("maxit", "thr"))
    maxit <- control$maxit
    thr <- control$thr
    # Testing 'maxit'
    if (!is.vector(maxit)) stop(sQuote("maxit"), " is not a vector")
    if (length(maxit) != 1) stop(sQuote("maxit"), " is not an object of length ", sQuote(1))
    if (abs(as.integer(maxit) - maxit) > 0) stop(sQuote("maxit"), " is not an object of type ", dQuote("integer"))
    if (maxit <= 0) stop(sQuote("maxit"), " is not a positive integer")
    # Testing 'thr'
    if (!is.vector(thr)) stop(sQuote("thr"), " is not a vector")
    if (length(thr) != 1) stop(sQuote("thr"), " is not an object of length ", sQuote(1))
    if (thr <= 0 ) stop(sQuote("thr"), " is not a positive value")
    # testing, for each column, if the entries belong to the interval [lo[m], up[m]]
    for (m in seq_len(p)) {
        lo.id <- which(Y[, m] < lo[m])
        up.id <- which(Y[, m] > up[m])
        if (length(lo.id) != 0L) {
            message("message: some entries in column ", sQuote(vnames[m]), " are below the lower censoring value. These entries are replaced with ", sQuote(lo[m]))
            Y[lo.id, m] <- lo[m]
        }
        if (length(up.id) != 0L) {
            message("message: some entries in column ", sQuote(vnames[m]), " are over the upper censoring value. These entries are replaced with ", sQuote(up[m]))
            Y[up.id, m] <- up[m]
        }
    }
    ######################
    # starting main code #
    ######################
    Yna <- -9999
    Y[is.na(Y)] <- Yna
    # storage.mode
    storage.mode(n) <- "integer"
    storage.mode(p) <- "integer"
    storage.mode(Y) <- "double"
    storage.mode(lo) <- "double"
    storage.mode(up) <- "double"
    storage.mode(Yna) <- "double"
    Info <- .Fortran(C_setup, n = n, p = p, Y = Y, lo = lo, up = up, Yna = Yna, R = matrix(0L, nrow = n + 1L, ncol = p + 1L),
                     startmis = 0L, order = seq_len(n))
    # if (Info$startmis == 0L) message("message: Matrix ", sQuote("Y"), " is fully observed")
    rownames(Info$Y) <- rownames(Info$Y)[Info$order]
    Info$startmis <- NULL
    Info$n <- NULL
    Info$p <- NULL
    Info$Yna <- NULL
    Info$R <- Info$R[-1L, ]
    id.pattern <- Info$R[, 1L]
    storage.mode(id.pattern) <- "logical"
    Info$R <- Info$R[, -1L]
    dimnames(Info$R) <- dimnames(Info$Y)
    if (any(id.pattern)) {
        Pattern <- t(apply(Info$R[id.pattern, , drop = FALSE], 1L, order))
        Pattern <- cbind(which(id.pattern), Pattern, t(apply(Info$R[id.pattern, , drop = FALSE], 1L, function(x) c(sum(x == 0), sum(x == 1), sum(x == 2)))))
    } else Pattern <- matrix(c(n + 1, seq_len(p), 0, 0, 0), nrow = 1L)
    colnames(Pattern) <- c("i", paste0("V", 1:p), "nmar", "nlc", "nrc")
    Info$Pattern <- Pattern
    Info$R[Info$R == 0] <- 9L
    Info$R[Info$R == 1] <- -1L
    Info$R[Info$R == 2] <- 1L
    Info$R[Info$R == 3] <- 0L
    storage.mode(maxit) <- "integer"
    storage.mode(thr) <- "double"
    out <- .Fortran(C_fitmcgm, n = n, p = p, Y = Info$Y , lo = Info$lo, up = Info$up,
                    R = Info$R, nstp = maxit, eps = thr, ym = double(p), yv = double(p),
                    conv = integer(1L))
    if (out$conv != 0) stop("Subroutine ", sQuote("fitmcgm"), " does not converge with code ", out$conv)
    Info$Y[Info$Y == Yna] <- NA
    Info$ym <- out$ym
    Info$yv <- out$yv
    names(Info$ym) <- vnames
    names(Info$yv) <- vnames
    row.order <- order(Info$order)
    if (q != 0) Z <- list(Y = Info$Y[row.order, , drop = FALSE], X = X)
    else Z <- list(Y = Info$Y[row.order, , drop = FALSE], X = NULL)
    Info$Y <- NULL
    Info$R <- Info$R[row.order, , drop = FALSE]
    Info$n <- n
    Info$p <- p
    Info$q <- q
    Z$Info <- Info
    class(Z) <- "datacggm"
    Z
}
