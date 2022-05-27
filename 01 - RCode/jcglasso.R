##### datajcggm #####
datajcggm <- function (Y, lo = -Inf, up = +Inf, X = NULL, 
                       control = list(maxit = 10000, thr = 0.0001)) {
  big <- .Machine$double.xmax
  thr <- big/2
  if (missing(Y)) stop(sQuote("Y"), " is missing")
  if (!is.matrix(Y)) stop(sQuote("Y"), " is not a matrix")
  if (any(!is.finite(Y[!is.na(Y)]))) stop("some element in 'Y' is not finite")
  n <- dim(Y)[1L]
  p <- dim(Y)[2L]
  if (is.null(colnames(Y))) {
    vnames <- paste0("Y", 1:p, sep = "")
    colnames(Y) <- vnames
  }
  else vnames <- colnames(Y)
  if (is.null(rownames(Y))) rownames(Y) <- seq_len(n)
  if (missing(lo)) loy <- rep(-big, p)
  else {
    if (!is.vector(lo)) stop(sQuote("lo"), " is not a vector")
    if (length(lo) == 1) loy <- rep(lo, p)
    else loy <- lo[1:p]
    id <- -big < loy[1:p] & loy[1:p] <= -thr
    if (any(id)) {
      loy[1:p][id] <- -big
      message("message: some entries in 'lo' are below the tolerance. These values are treated as -Inf")
    }
    id <- loy[1:p] == -Inf
    if (any(id)) loy[1:p][id] <- -big
  }
  names(loy)[1:p] <- vnames
  if (missing(up)) upy <- rep(big, p)
  else {
    if (!is.vector(up)) stop(sQuote("up"), " is not a vector")
    if (length(up) == 1) upy <- rep(up, p)
    else upy <- up[1:p]
    id <- thr <= upy[1:p] & upy[1:p] < big
    if (any(id)) {
      upy[1:p][id] <- big
      message("message: some entries in 'up' are over the tolerance. These values are treated as +Inf")
    }
    id <- upy[1:p] == Inf
    if (any(id)) upy[1:p][id] <- big
  }
  names(upy)[1:p] <- vnames
  if (!is.null(X)) {
    # if (any(is.na(X))) stop("Missing values in ", sQuote("X"), " are not allowed")
    if (is.vector(X)) X <- matrix(X, ncol = 1)
    if (length(dim(X)) != 2L) stop(sQuote("X"), "is not a matrix like object")
    if (dim(Y)[1L] != dim(X)[1L]) stop(sQuote("X"), " and ", sQuote("Y"), " have a different number of rows")
    if (is.null(colnames(X))) colnames(X) <- paste0("X", seq_len(dim(X)[2L]), sep = "")
    if (is.null(rownames(X))) rownames(X) <- seq_len(n)
    # if (!is.data.frame(X)) X <- as.data.frame(X)
    if (!is.matrix(X)) stop(sQuote("X"), " is not a matrix")
    q <- dim(X)[2L]
    xnames <- colnames(X)
    if (missing(lo)) lox <- c(loy, rep(-big, q))
    else {
      if (!is.vector(lo)) stop(sQuote("lo"), " is not a vector")
      if (length(lo) == 1) lox <- c(loy, rep(lo, q))
      else lox <- c(loy, lo[(p+1):(p+q)])
      id <- -big < lox[(p+1):(p+q)] & lox[(p+1):(p+q)] <= -thr
      if (any(id)) {
        lox[(p+1):(p+q)][id] <- -big
        message("message: some entries in 'lo' are below the tolerance. These values are treated as -Inf")
      }
      id <- lox[(p+1):(p+q)] == -Inf
      if (any(id)) lox[(p+1):(p+q)][id] <- -big
    }
    names(lox)[(p+1):(p+q)] <- xnames
    if (missing(up)) upx <- c(upy, rep(big, q))
    else {
      if (!is.vector(up)) stop(sQuote("up"), " is not a vector")
      if (length(up) == 1) upx <- c(upy, rep(up, q))
      else upx <- c(upy, up[(p+1):(p+q)])
      id <- thr <= upx[(p+1):(p+q)] & upx[(p+1):(p+q)] < big
      if (any(id)) {
        upx[(p+1):(p+q)][id] <- big
        message("message: some entries in 'up' are over the tolerance. These values are treated as +Inf")
      }
      id <- upx[(p+1):(p+q)] == Inf
      if (any(id)) upx[(p+1):(p+q)][id] <- big
    }
    names(upx)[(p+1):(p+q)] <- xnames
  }
  else q <- 0L
  if (!is.list(control)) stop(sQuote("control"), "is not an object of type ", sQuote("list"))
  if (is.null(names(control))) stop(sQuote("control"), "is not a named list")
  names(control) <- sapply(names(control), match.arg, choice = c("maxit", "thr"))
  maxit <- control$maxit
  thr <- control$thr
  if (!is.vector(maxit)) stop(sQuote("maxit"), " is not a vector")
  if (length(maxit) != 1) stop(sQuote("maxit"), " is not an object of length ", sQuote(1))
  if (abs(as.integer(maxit) - maxit) > 0) stop(sQuote("maxit"), " is not an object of type ", dQuote("integer"))
  if (maxit <= 0) stop(sQuote("maxit"), " is not a positive integer")
  if (!is.vector(thr)) stop(sQuote("thr"), " is not a vector")
  if (length(thr) != 1) stop(sQuote("thr"), " is not an object of length ", sQuote(1))
  if (thr <= 0) stop(sQuote("thr"), " is not a positive value")
  Z <- if(q != 0) cbind(Y, X) else Y
  znames <- if(q != 0) c(vnames, xnames) else vnames
  lo <- if(q != 0) lox else loy
  up <- if(q != 0) upx else upy
  if (!all(lo < up)) stop(sQuote("lo"), " is not less than ", sQuote("up"))
  for (m in seq_len(p + q)) {
    lo.id <- which(Z[, m] < lo[m])
    up.id <- which(Z[, m] > up[m])
    if (length(lo.id) != 0L) {
      message("message: some entries in column ", sQuote(znames[m]), 
              " are below the lower censoring value. These entries are replaced with ", 
              sQuote(lo[m]))
      Z[lo.id, m] <- lo[m]
    }
    if (length(up.id) != 0L) {
      message("message: some entries in column ", sQuote(znames[m]), 
              " are over the upper censoring value. These entries are replaced with ", 
              sQuote(up[m]))
      Z[up.id, m] <- up[m]
    }
  }
  Zna <- -9999
  Z[is.na(Z)] <- Zna
  storage.mode(n) <- "integer"
  storage.mode(p) <- "integer"
  storage.mode(Y) <- "double"
  storage.mode(lo) <- "double"
  storage.mode(up) <- "double"
  storage.mode(Zna) <- "double"
  Info <- .Fortran(cglasso:::C_setup, n = n, p = p + q, Y = Z, lo = lo, up = up, 
                   Yna = Zna, R = matrix(0L, nrow = n + 1L, ncol = p + q + 1L), 
                   startmis = 0L, order = seq_len(n))
  names(Info)[3] <- "Z"
  rownames(Info$Z) <- rownames(Info$Z)[Info$order]
  Info$startmis <- NULL
  Info$n <- NULL
  Info$p <- NULL
  Info$Yna <- NULL
  Info$R <- Info$R[-1L, ]
  id.pattern <- Info$R[, 1L]
  storage.mode(id.pattern) <- "logical"
  Info$R <- Info$R[, -1L]
  dimnames(Info$R) <- dimnames(Info$Z)
  if (any(id.pattern)) {
    Pattern <- t(apply(Info$R[id.pattern, , drop = FALSE], 1L, order))
    Pattern <- cbind(which(id.pattern), Pattern, 
                     t(apply(Info$R[id.pattern, , drop = FALSE], 1L, function(x) c(sum(x == 0), sum(x == 1), sum(x == 2)))))
  }
  else Pattern <- matrix(c(n + 1, seq_len(p + q), 0, 0, 0), nrow = 1L)
  colnames(Pattern) <- c("i", paste0("V", 1:(p + q) ), "nmar", "nlc", "nrc")
  Info$Pattern <- Pattern
  Info$R[Info$R == 0] <- 9L
  Info$R[Info$R == 1] <- -1L
  Info$R[Info$R == 2] <- 1L
  Info$R[Info$R == 3] <- 0L
  storage.mode(maxit) <- "integer"
  storage.mode(thr) <- "double"
  out <- .Fortran(cglasso:::C_fitmcgm, n = n, p = p + q, Y = Info$Z, lo = Info$lo, 
                  up = Info$up, R = Info$R, nstp = maxit, eps = thr, ym = double(p + q), 
                  yv = double(p + q), conv = integer(1L))
  if (out$conv != 0) 
    stop("Subroutine ", sQuote("fitmcgm"), " does not converge with code ", out$conv)
  Info$Z[Info$Z == Zna] <- NA
  Info$zm <- out$ym
  Info$zv <- out$yv
  names(Info$zm) <- znames
  names(Info$zv) <- znames
  row.order <- order(Info$order)
  if (q != 0) 
    Z <- list(Y = Info$Z[row.order, 1:p, drop = FALSE], X = Info$Z[row.order, (p+1):(p+q), drop = FALSE])
  else Z <- list(Y = Info$Z[row.order, , drop = FALSE], X = NULL)
  Info$Z <- NULL
  Info$R <- Info$R[row.order, , drop = FALSE]
  Info$n <- n
  Info$p <- p
  Info$q <- q
  Z$Info <- Info
  class(Z) <- "datajcggm"
  Z
}

is.datajcggm <- function (x) inherits(x, "datajcggm")

print.datajcggm <- function (x, digits = 3L, n = 10L, width = getOption("width"), ...) {
  if (!is.infinite(n)) {
    if (abs(as.integer(n) - n) > 0) 
      stop(sQuote("n"), " is not a positive integer value")
    if (n <= 0) 
      stop(sQuote("n"), " is not a strictly positive integer value")
  }
  if (!is.infinite(width)) {
    if (abs(as.integer(width) - width) > 0) 
      stop(sQuote("width"), " is not a positive integer value")
    if (width <= 0) 
      stop(sQuote("width"), " is not a strictly positive integer value")
  }
  out <- as.character.datajcggm(x, digits = digits, ...)
  out$Y <- cglasso:::fit_matrix(out$Y, n = n, width = width)
  p <- nresp(x); q <- npred(x)
  cat("Printing", sQuote("datajcggm"), "object\n\n")
  cat("Y:", nobs(x), "x", p, "matrix\n\n")
  print(out$Y, quote = FALSE)
  dn <- nobs(x) - n
  dp <- p - dim(out$Y)[2L]
  if (dn > 0 | dp > 0) {
    msg <- paste0("# with ", ifelse(dn > 0, paste0(dn, " more ", ifelse(dn > 1, "rows", "row")), ""), 
                  ifelse(dn > 0 & dp > 0, ", and ", ""), 
                  ifelse(dp > 0, paste0(dp, " more ", ifelse(dp > 1, "variables", "variable")), ""))
    cat(msg)
  }
  if (q == 0L) 
    cat("\n\nX: NULL")
  else {
    # out <- format(getMatrix(x, "X"), digits = digits, ...)
    out$X <- fit_matrix(out$X, n = n, width = width)
    cat("\n\nX:", nobs(x), "x", npred(x), "matrix\n\n")
    print(out$X, quote = FALSE)
    dq <- npred(x) - dim(out$X)[2L]
    if (dn > 0 | dq > 0) {
      msg <- paste0("# with ", ifelse(dn > 0, paste0(dn, " more ", ifelse(dn > 1, "rows", "row")), ""), 
                    ifelse(dn > 0 & dq > 0, ", and ", ""), 
                    ifelse(dq > 0, paste0(dq, " more ", ifelse(dq > 1, "variables", "variable")), ""))
      cat(msg)
    }
  }
  cat("\n")
}

summary.datajcggm <- function (object, n, quantile.type = 7L, digits = 3L, quote = FALSE, ...) {
  if (missing(n)) n <- Inf
  Y <- getMatrix2(object, "Y")
  X <- getMatrix2(object, "X")
  R <- event2(object)
  ynm <- colNames2(object)$Y
  p <- nresp(object)
  q <- npred(object)
  lo <- unlist(lower2(object))
  up <- unlist(upper2(object))
  out <- list(Y = NULL, X = NULL)
  lcs <- paste0(format(apply(R == -1, 2L, function(x) mean(x) * 100), digits = digits), "%")
  rcs <- paste0(format(apply(R == +1, 2L, function(x) mean(x) * 100), digits = digits), "%")
  nas <- paste0(format(apply(R == +9, 2L, function(x) mean(x) * 100), digits = digits), "%")
  qq <- matrix(0, nrow = p, ncol = 6L)
  for (m in seq_len(p)) {
    obs <- Y[R[, m] == 0, m]
    qq[m, -4L] <- stats::quantile(obs, names = FALSE, type = quantile.type)
    qq[m, 4L] <- mean(obs)
  }
  tbl <- apply(cbind(lo[1:p], qq, up[1:p]), 2L, format, digits = digits)
  tbl <- cbind(tbl, nas[1:p], lcs[1:p], rcs[1:p])
  colnames(tbl) <- c("Lower", "Min.", "1st Qu.", "Median", 
                     "Mean", "3rd Qu.", "Max.", "Upper", 
                     "NA%", "LC%", "RC%")
  rownames(tbl) <- ynm
  attr(tbl, "class") <- "table"
  out$Y <- tbl
  if (q > 0L) {
    xnm <- colNames2(object)$X
    qq <- matrix(0, nrow = q, ncol = 6L)
    for (m in seq_len(q)) {
      obs <- X[R[, p + m] == 0, m]
      qq[m, -4L] <- stats::quantile(obs, names = FALSE, type = quantile.type)
      qq[m, 4L] <- mean(obs)
    }
    tbl <- apply(cbind(lo[(p+1):(p+q)], qq, up[(p+1):(p+q)]), 2L, format, digits = digits)
    tbl <- cbind(tbl, nas[(p+1):(p+q)], lcs[(p+1):(p+q)], rcs[(p+1):(p+q)])
    colnames(tbl) <- c("Lower", "Min.", "1st Qu.", "Median", 
                       "Mean", "3rd Qu.", "Max.", "Upper", 
                       "NA%", "LC%", "RC%")
    rownames(tbl) <- xnm
    attr(tbl, "class") <- "table"
    out$X <- tbl
  }
  
  dvar <- dim(out$Y)[1L] - n
  cat("Y:\n")
  print(head(out$Y, n), quote = quote, ...)
  if (dvar > 0) 
    cat(paste0("# with ", ifelse(dvar > 0, paste0(dvar, " more ", ifelse(dvar > 1, "variables", "variable")))))
  if (q == 0L) 
    cat("\n\nX: NULL\n\n")
  else {
    dvar <- dim(out$X)[1L] - n
    cat("\n\nX:\n")
    print(head(out$X, n), quote = quote, ...)
    if (dvar > 0) 
      cat(paste0("# with ", ifelse(dvar > 0, paste0(dvar, " more ", ifelse(dvar > 1, "variables", "variable")))))
  }
  cat("\n")
  invisible(out)
}

dim.datajcggm <- function (x) {
  out <- list(Y = dim(getMatrix2(x, "Y")), X = dim(getMatrix2(x, "X")))
  out
}

nresp.datajcggm <- function (object, ...) object$Info$p

nobs.datajcggm <- function (object, ...) object$Info$n

npred.datajcggm <- function (object, ...) object$Info$q

dimnames.datajcggm <- function (x) {
  Y <- getMatrix2(x, "Y")
  out <- list(Y = dimnames(Y), X = NULL)
  if (npred(x) != 0) {
    X <- getMatrix2(x, "X")
    out$X <- dimnames(X)
  }
  out
}

rowNames2 <- function (x) {
  if (!is.datajcggm(x)) 
    stop(sQuote(x), " is not an object of class ", sQuote("datajcggm"))
  Y <- getMatrix2(x, "Y")
  X <- getMatrix2(x, "X")
  if (is.null(x$X)) 
    out <- list(Y = rownames(Y), X = NULL)
  else out <- list(Y = rownames(Y), X = rownames(X))
  out
}

colNames2 <- function (x) {
  if (!is.datajcggm(x)) 
    stop(sQuote(x), " is not an object of class ", sQuote("datajcggm"))
  Y <- getMatrix2(x, "Y")
  X <- getMatrix2(x, "X")
  if (is.null(x$X)) 
    out <- list(Y = colnames(Y), X = NULL)
  else out <- list(Y = colnames(Y), X = colnames(X))
  out
}

getMatrix2 <- function (x, name = c("Y", "X", "both"), ordered = FALSE) {
  if (!is.datajcggm(x)) 
    stop(sQuote(x), " is not an object of class ", sQuote("datajcggm"))
  if (!is.logical(ordered)) 
    stop(sQuote(ordered), " id not an object of type ", sQuote("logical"))
  name <- match.arg(name)
  if (ordered) {
    row.order <- x$Info$order
    Y <- x$Y[row.order, , drop = FALSE]
    X <- x$X[row.order, , drop = FALSE]
  }
  else {
    Y <- x$Y
    X <- x$X
  }
  out <- switch(name, Y = Y, X = X, both = list(Y = Y, X = X))
  out
}

event2 <- function (x, ordered = FALSE) {
  if (!is.datajcggm(x)) 
    stop("'x' is not an object of class 'datajcggm'")
  if (!is.logical(ordered)) 
    stop(sQuote(ordered), " id not an object of type ", sQuote("logical"))
  R <- x$Info$R
  if (ordered) 
    R <- R[x$Info$order, , drop = FALSE]
  R
}

lower2 <- function (x) {
  if (!is.datajcggm(x)) 
    stop(sQuote(x), " is not an object of class ", sQuote("datajcggm"))
  small <- -.Machine$double.xmax
  p <- nresp(x); q <- npred(x)
  out <- x$Info$lo
  out[out == small] <- -Inf
  if (q != 0L) {
    out <- list(Y = out[1:p], X = out[(p+1):(p+q)])
  } else out <- list(Y = out[1:p], X = NULL)
  out
}

upper2 <- function (x) {
  if (!is.datajcggm(x)) 
    stop(sQuote(x), " is not an object of class ", sQuote("datajcggm"))
  big <- .Machine$double.xmax
  p <- nresp(x); q <- npred(x)
  out <- x$Info$up
  out[out == big] <- +Inf
  if (q != 0L) {
    out <- list(Y = out[1:p], X = out[(p+1):(p+q)])
  } else out <- list(Y = out[1:p], X = NULL)
  out
}

ColMeans2 <- function (x) {
  if (!is.datajcggm(x)) 
    stop(sQuote("x"), " is not a ", sQuote("datajcggm"), " object")
  p <- nresp(x); q <- npred(x)
  if (q != 0L) {
    out <- list(Y = x$Info$zm[1:p], X = x$Info$zm[(p+1):(p+q)])
  } else out <- list(Y = x$Info$zm, X = NULL)
  out
}

ColVars2 <- function (x) {
  if (!is.datajcggm(x)) 
    stop(sQuote("x"), " is not a ", sQuote("datajcggm"), " object")
  p <- nresp(x); q <- npred(x)
  if (q != 0L) {
    out <- list(Y = x$Info$zv[1:p], X = x$Info$zv[(p+1):(p+q)])
  } else out <- list(Y = x$Info$zv, X = NULL)
  out
}

as.character.datajcggm <- function (x, digits, ...) {
  p <- nresp(x); q <- npred(x)
  if(q != 0) {
    Z <- getMatrix2(x, "both")
    Z <- format(cbind(Z$Y, Z$X), digits = digits, ...)
  } else Z <- format(getMatrix2(x, "Y"), digits = digits, ...)
  marker <- event2(x)
  Z[marker == 1] <- paste0(Z[marker == 1], "+")
  Z[marker == -1] <- paste0(Z[marker == -1], "-")
  if(q != 0) 
    list(Y = Z[, 1:p, drop = FALSE], X = Z[, (p+1):(p+q), drop = FALSE])
  else list(Y = Z[, 1:p, drop = FALSE], X = NULL)
}

##### jcglasso functions ######
jcglasso <- function (formula, data, subset, contrasts = NULL, diagonal = FALSE, weights.B = NULL, 
                      weights.Tht = NULL, penalty = c("group", "fused"), 
                      nrho, rho.min.ratio, rho, nlambda, lambda.min.ratio, lambda, nu = NULL, alpha = 0.5,
                      maxit.em = 1E+5, thr.em = 1E-3, maxit.bcd = 1E+5, thr.bcd = 1E-04, 
                      trace = 0L, offset = NULL, covar2corr = FALSE, truncate = 1E-6) {
  this.call <- match.call()
  penalty <- match.arg(penalty)
  zero <- 1e-06
  
  # if(missing(rho1) | missing(rho2)) stop("Please, insert non-negative values for rho1 and rho2!!")
  if(is.datajcggm(data)) data <- list(data)
  K <- length(data)
  
  X <- Y <- Z <- vector(mode = "list", length = K)
  for(k in seq_len(K)) {
    # testing 'formula'
    if (missing(formula)) formula <- . ~ . # stop(sQuote("formula"), " is missing")
    # testing 'data'
    if (!is.datajcggm(data[[k]])) stop(sQuote("data"), " is not an object of class ", sQuote("datajcggm"))
    
    ##### formula2datacggm #####
    # testing LHS 'formula'
    fmlTerms <- function(fml, pos) paste(deparse(fml[[pos]], width.cutoff = 500L), collapse = " ")
    if (fmlTerms(formula, 2L) == ".")
      formula <- formula(paste0(paste0("cbind(", paste(colNames2(data[[k]])$Y, collapse = ", "), ")"), " ~ ", fmlTerms(formula, 3L)))
    if (as.character(formula[[2L]])[1L] != "cbind")
      stop("Please use ", sQuote("cbind"), " to specify LHS in ", sQuote("formula"), " object")
    Y.LHS <- unlist(lapply(formula[[2L]][-1L], deparse))
    if (any(table(Y.LHS) > 1L))
      stop("repeated response variables are not permitted in ", sQuote("formula"), " object")
    noVars <- !is.element(Y.LHS, colNames2(data[[k]])$Y)
    if (any(noVars)) stop("Following variables are not stored as rensponse variables: ", Y.LHS[noVars])
    if (is.null(getMatrix2(data[[k]], name = "X"))) {
      if (fmlTerms(formula, 3L) == ".") formula <- update(formula, . ~ 1)
      mt <- terms(formula)
      if (length(attr(mt, which = "term.labels")) != 0L)
        stop("Predictors are not stored in ", sQuote("data"))
      if (!as.logical(attr(mt, "intercept"))) {
        warning("Current version does not fit models without intercept term, thus it is added to the current formula.")
        formula <- update(formula, . ~ 1)
        mt <- terms(formula)
      }
      data.df <- data.frame(getMatrix2(data[[k]], name = "Y")[, Y.LHS, drop = FALSE])
      nResp <- length(Y.LHS)
      nPred <- 0L
    } 
    else {
      FML.RHS <- attr(terms(formula, data = data.frame(getMatrix2(data[[k]], name = "X"))), "term.labels")
      if (length(FML.RHS) > 0L) {
        FML.RHS <- paste(FML.RHS, collapse = " + ")
        formula <- formula(paste0(fmlTerms(formula, 2L), " ~ ", FML.RHS))
        X.RHS <- all.vars(formula[[3L]])
        noVars <- !is.element(X.RHS, colNames2(data[[k]])$X)
        if (any(noVars)) stop("Following variables are not stored as predictors: ", X.RHS[noVars])
        data.df <- data.frame(getMatrix2(data[[k]], name = "Y")[, Y.LHS, drop = FALSE], data.frame(getMatrix2(data[[k]], name = "X")[, X.RHS, drop = FALSE]))
        nPred <- length(X.RHS)
      } else {
        data.df <- data.frame(getMatrix2(data[[k]], name = "Y")[, Y.LHS, drop = FALSE])
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
    
    # creating model.frame #
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
    Y[[k]] <- model.response(mf, type = "numeric")
    lo <- lower2(data[[k]])$Y[colnames(Y[[k]])]
    up <- upper2(data[[k]])$Y[colnames(Y[[k]])]
    if (length(attr(mt, which = "term.labels")) == 0L) {
      Z[[k]] <- datajcggm(Y = Y[[k]], lo = lo, up = up)
    } 
    else {
      X[[k]] <- model.matrix(mt, mf, contrasts)[, -1L, drop = FALSE]
      lo <- c(lo, lower2(data[[k]])$X[colnames(X[[k]])])
      up <- c(up, upper2(data[[k]])$X[colnames(X[[k]])])
      Z[[k]] <- datajcggm(Y = Y[[k]], X = X[[k]], lo = lo, up = up)
    }
  }
  
  X.null <- all(sapply(X, is.null))
  n <- sapply(Z, nobs)
  p <- sapply(Z, nresp)
  if(all(p == p[1])) { p <- p[1] } else { stop("variabili diverse!??!") }
  q <- sapply(Z, npred)
  if(all(q == q[1])) { q <- q[1] } else { stop("variabili diverse!??!") }
  if (p == 1L) stop("number of response variables is equal to ", sQuote(1))
  
  xnames <- colNames2(Z[[1]])$X
  ynames <- colNames2(Z[[1]])$Y
  
  if (!is.vector(diagonal)) stop(sQuote("diagonal"), " is not a vector")
  if (length(diagonal) > 1L) stop(sQuote("diagonal"), " is not an object of length ", sQuote(1))
  if (!is.logical(diagonal)) stop(sQuote("diagonal"), " is not a logical object")
  if (!is.null(weights.B)) {
    for(k in seq_len(K)) {
      if (X.null) stop("Argument ", sQuote("weights.B"), " can not be used because ", sQuote("X"), " is missing")
      if (!is.matrix(weights.B[[k]])) stop(sQuote("weights.B"), " is not a matrix")
      weights.B.dim <- c(q, p)
      if (any(dim(weights.B[[k]]) != weights.B.dim)) 
        stop(sQuote("weights.B"), " is not a matrix with 'dim' attribute equal to ", 
             weights.B.dim[1L], " ", weights.B.dim[2L])
      if (any(weights.B[[k]] < 0)) 
        stop("negative weights in ", sQuote("weights.B "), " are not allowed")
      if (all(weights.B[[k]] == 0)) 
        stop("all entries in ", sQuote("weights.B "), " are equal to zero")
      weights.B[[k]][weights.B[[k]] == +Inf] <- .Machine$double.xmax
      rownames(weights.B[[k]]) <- xnames
      colnames(weights.B[[k]]) <- ynames
    }
  }
  else {
    if (!X.null) {
      weights.B <- vector(mode = "list", length = K)
      for(k in seq_len(K)) 
        weights.B[[k]] <- matrix(1, nrow = q, ncol = p, dimnames = list(xnames, ynames))
    }
  }
  if (!is.null(weights.Tht)) {
    for(k in seq_len(K)) {
      if (!is.matrix(weights.Tht[[k]])) stop(sQuote("weights.Tht"), " is not a matrix")
      weights.Tht.dim <- c(p, p)
      if (any(dim(weights.Tht[[k]]) != weights.Tht.dim))
        stop(sQuote("weights.Tht"), " is not a matrix with dim attribute equal to ",
             weights.Tht.dim[1L], " ", weights.Tht.dim[2L])
      if (any(weights.Tht[[k]] < 0))
        stop("negative weights in ", sQuote("weights.Tht "), " are not allowed")
      if (all(weights.Tht[[k]] == 0))
        stop("all entries in ", sQuote("weights.Tht "), " are equal to zero")
      weights.Tht[[k]][weights.Tht[[k]] == +Inf] <- .Machine$double.xmax
      if (!all(weights.Tht[[k]] == t(weights.Tht[[k]]))) stop(sQuote("weights.Tht"), " is not a symmetric matrix")
      if (!diagonal & any(diag(weights.Tht[[k]]) != 0))
        stop(sQuote("diagonal = FALSE"), " but some diagonal entry in ",
             sQuote("weights.Tht"), " is not zero")
      if (diagonal & all(diag(weights.Tht[[k]]) == 0))
        stop(sQuote("diagonal = TRUE"), " but all diagonal entries in ",
             sQuote("weights.Tht"), " are zero")
      rownames(weights.Tht[[k]]) <- ynames
      colnames(weights.Tht[[k]]) <- ynames
    }
  }
  else {
    weights.Tht <- vector(mode = "list", length = K)
    for(k in seq_len(K)) {
      weights.Tht[[k]] <- matrix(1, nrow = p, ncol = p)
      diag(weights.Tht[[k]]) <- ifelse(diagonal, 1, 0)
      rownames(weights.Tht[[k]]) <- ynames
      colnames(weights.Tht[[k]]) <- ynames
    }
  }
  
  if (!missing(nlambda) & !missing(lambda)) 
    warning("Argumnet ", sQuote("nlambda"), " is overwritten using length(lambda)")
  if (!missing(lambda.min.ratio) & !missing(lambda)) 
    warning("Argumnet ", sQuote("lambda.min.ratio"), " is overwritten using lambda")
  if (!missing(nlambda)) {
    if (X.null)
      stop("Argument ", sQuote("nlambda"), " can not be used because ", sQuote("X"), " is missing")
    if (!is.vector(nlambda)) 
      stop(sQuote("nlambda"), " is not a vector")
    if (length(nlambda) != 1) 
      stop(sQuote("nlambda"), " is not an object of length ", sQuote(1))
    if (abs(as.integer(nlambda) - nlambda) > 0) 
      stop(sQuote("nlambda"), " is not an integer value")
    if (nlambda <= 0) 
      stop(sQuote("nrnlambdaho1"), " is not a strictly positive integer value")
  }
  else nlambda <- ifelse(X.null, 1L, 10L)
  if (!missing(lambda.min.ratio)) {
    if (X.null)
      stop("Argument ", sQuote("lambda.min.ratio"), " can not be used because ", sQuote("X"), " is missing")
    if (!is.vector(lambda.min.ratio)) 
      stop(sQuote("lambda.min.ratio"), " is not a vector")
    if (length(lambda.min.ratio) != 1) 
      stop(sQuote("lambda.min.ratio"), " is not an object of length ", sQuote(1))
    if (lambda.min.ratio < 0 | lambda.min.ratio > 1) 
      stop(sQuote("lambda.min.ratio"), " does not belong to the closed interval [0, 1]")
  }
  else {
    if(X.null) lambda.min.ratio <- 0
    else lambda.min.ratio <- ifelse(all(p < n), zero, 0.01)
  }
  if (!missing(lambda)) {
    if (X.null)
      stop("Argument ", sQuote("lambda"), " can not be used because ", sQuote("X"), " is missing")
    if (!is.vector(lambda)) 
      stop(sQuote("lambda"), " is not a vector")
    if (any(lambda < 0)) 
      stop("some entry in ", sQuote("lambda"), " is negative")
    lambda <- sort(lambda, decreasing = TRUE)
    id <- lambda <= zero
    if (any(id)) 
      lambda <- c(lambda[!id], zero)
    nlambda <- length(lambda)
    lambda.min.ratio <- lambda[nlambda]/lambda[1L]
  }
  else {
    if(X.null) lambda <- 0
    else lambda <- double(nlambda)
  }
  
  if (!missing(nrho) & !missing(rho)) 
    warning("Argumnet ", sQuote("nrho"), " is overwritten using length(rho)")
  if (!missing(rho.min.ratio) & !missing(rho)) 
    warning("Argumnet ", sQuote("rho.min.ratio"), " is overwritten using rho")
  if (!missing(nrho)) {
    if (!is.vector(nrho)) 
      stop(sQuote("nrho"), " is not a vector")
    if (length(nrho) != 1) 
      stop(sQuote("nrho"), " is not an object of length ", sQuote(1))
    if (abs(as.integer(nrho) - nrho) > 0) 
      stop(sQuote("nrho"), " is not an integer value")
    if (nrho <= 0) 
      stop(sQuote("nrho"), " is not a strictly positive integer value")
  }
  else nrho <- 10L
  if (!missing(rho.min.ratio)) {
    if (!is.vector(rho.min.ratio)) 
      stop(sQuote("rho.min.ratio"), " is not a vector")
    if (length(rho.min.ratio) != 1) 
      stop(sQuote("rho.min.ratio"), " is not an object of length ", sQuote(1))
    if (rho.min.ratio < 0 | rho.min.ratio > 1) 
      stop(sQuote("rho.min.ratio"), " does not belong to the closed interval [0, 1]")
  }
  else rho.min.ratio <- ifelse(all(p < n), zero, 0.01)
  if (!missing(rho)) {
    if (!is.vector(rho)) 
      stop(sQuote("rho"), " is not a vector")
    if (any(rho < 0)) 
      stop("some entry in ", sQuote("rho"), " is negative")
    rho <- sort(rho, decreasing = TRUE)
    id <- rho <= zero
    if (any(id)) 
      rho <- c(rho[!id], zero)
    nrho <- length(rho)
    rho.min.ratio <- rho[nrho]/rho[1L]
  }
  else rho <- double(nrho)
  
  if (!is.vector(maxit.em)) stop(sQuote("maxit.em"), " is not a vector")
  if (length(maxit.em) != 1) 
    stop(sQuote("maxit.em"), " is not an object of length ", sQuote(1))
  if (abs(as.integer(maxit.em) - maxit.em) > 0) 
    stop(sQuote("maxit.em"), " is not an object of type ", dQuote("integer"))
  if (maxit.em <= 0) stop(sQuote("maxit.em"), " is not a positive integer")
  if (!is.vector(thr.em)) stop(sQuote("thr.em"), " is not a vector")
  if (length(thr.em) != 1) 
    stop(sQuote("thr.em"), " is not an object of length ", sQuote(1))
  if (thr.em <= 0) stop(sQuote("thr.em"), " is not a positive value")
  if (!is.vector(maxit.bcd)) stop(sQuote("maxit.bcd"), " is not a vector")
  if (length(maxit.bcd) != 1) 
    stop(sQuote("maxit.bcd"), " is not an object of length ", sQuote(1))
  if (abs(as.integer(maxit.bcd) - maxit.bcd) > 0) 
    stop(sQuote("maxit.bcd"), " is not an object of type ", dQuote("integer"))
  if (maxit.bcd <= 0) stop(sQuote("maxit.bcd"), " is not a positive integer")
  if (!is.vector(thr.bcd)) stop(sQuote("thr.bcd"), " is not a vector")
  if (length(thr.bcd) != 1) 
    stop(sQuote("thr.bcd"), " is not an object of length ", sQuote(1))
  if (thr.bcd <= 0) stop(sQuote("thr.bcd"), " is not a positive value")
  if (!is.vector(trace)) stop(sQuote("trace"), " is not a vector")
  if (length(trace) != 1) 
    stop(sQuote("trace"), " is not an object of length ", sQuote(1))
  if (is.logical(trace)) 
    stop(sQuote("trace"), " is not an object of type ", dQuote("integer"))
  if (abs(as.integer(trace) - trace) > 0) 
    stop(sQuote("trace"), " is not an object of type ", dQuote("integer"))
  if (!is.element(trace, c(0L, 1L, 2L))) 
    stop("not allowed value in ", sQuote("trace"), ". Please, choice ", 
         sQuote("0"), ", ", sQuote("1"), " or ", sQuote("2"))
  out <- jcglasso.fit(Z = Z, diagonal = diagonal, weights.B = weights.B, weights.Tht = weights.Tht, 
                      penalty = penalty, nrho = nrho, rho.min.ratio = rho.min.ratio, rho = rho, 
                      nlambda = nlambda, lambda.min.ratio = lambda.min.ratio, lambda = lambda, 
                      nu = nu, alpha = alpha, maxit.em = maxit.em, thr.em = thr.em, maxit.bcd = maxit.bcd, 
                      thr.bcd = thr.bcd, trace = trace, offset = offset, covar2corr = covar2corr, truncate = truncate)
  if (out$conv != "Ok") {
    msg <- paste("jcglasso does not converge. Subroutine",
                 sQuote(out$subrout), "returns message:\n", sQuote(out$conv))
    warning(msg)
  }
  for(k in seq_len(K)) {
    if (!X.null) weights.B[[k]][weights.B[[k]] == .Machine$double.xmax] <- +Inf
    out$wTht[[k]][out$wTht[[k]] == .Machine$double.xmax] <- +Inf
  }
  InfoStructure <- list(Adj = out$Adj, 
                        Id = out$Id, np = out$nP, InfoP = out$InfoP,
                        ncomp = out$ncomp, Ck = out$Ck, pk = out$pk,
                        lo = out$lo, up = out$up, zm = out$zm, zv = out$zv,
                        id_X = out$id_X, id_Y = out$id_Y)
  if (all(sapply(seq_len(K), function(k) Z[[k]]$Info$Pattern[1L, "i"] > n[k])))
    model <- "joint glasso"
  else {
    id.mar <- any(sapply(seq_len(K), function(k) Z[[k]]$Info$Pattern[1L, "nmar"] > 0))
    id.cens <- any(sapply(seq_len(K), function(k) Z[[k]]$Info$Pattern[1L, c("nlc", "nrc")] > 0))
    if (id.mar & !id.cens)
      model <- "joint missglasso"
    if (!id.mar & id.cens)
      model <- "joint censored glasso"
    if (id.mar & id.cens)
      model <- "joint hybrid glasso"
  }
  if (q > 0L) model <- paste("conditional", model)
  
  out.jcglasso <- list(call = this.call, Zipt = out$Zipt, B = out$B,
                       mu = out$mu, R = out$R, S = out$S, Sgm = out$Sgm, Tht = out$Tht,
                       Omega = out$Omega, dfB = out$dfB, dfTht = out$dfTht, InfoStructure = InfoStructure,
                       nit = out$nit, Z = Z, diagonal = diagonal, 
                       weights.B = weights.B, weights.Tht = out$wTht,
                       nlambda = nlambda, lambda.min.ratio = lambda.min.ratio, lambda = out$lambda, 
                       nrho = nrho, rho.min.ratio = rho.min.ratio, rho = out$rho, 
                       nu = nu, alpha = alpha, lambda2 = out$lambda2, rho2 = out$rho2, 
                       connected = out$connected, penalty = penalty,
                       model = model, maxit.em = maxit.em, thr.em = thr.em,
                       maxit.bcd = maxit.bcd, thr.bcd = thr.bcd, conv = out$conv,
                       offset = out$offset, covar2corr = covar2corr, truncate = truncate,
                       subrout = out$subrout, trace = trace, nobs = n, nresp = p, npred = q)
  
  # out.cjglasso <- list(call = this.call, Yipt = out$Yipt, Xipt = out$Xipt, B = out$B,
  #                      mu = out$mu, R = out$R, S = out$S, Sgm = out$Sgm, Tht = out$Tht,
  #                      Sxx = out$Sxx, Sxy = out$Sxy, Sgmxx = out$Sgmxx, Sgmxy = out$Sgmxy, 
  #                      Thtxx = out$Thtxx, Thtxy = out$Thtxy,
  #                      dfB = out$dfB, dfTht = out$dfTht, InfoStructure = InfoStructure,
  #                      nit = out$nit, Z = Z, diagonal = diagonal, weights.B = weights.B,
  #                      weights.Tht = weights.Tht,
  #                      nlambda = nlambda, lambda.min.ratio = lambda.min.ratio, lambda = out$lambda, 
  #                      nrho = nrho, rho.min.ratio = rho.min.ratio, rho = out$rho, rho.x = rho.x, alpha = alpha,
  #                      lambda2 = out$lambda2, rho2 = out$rho2, 
  #                      connected = out$connected, penalty = penalty,
  #                      model = model, maxit.em = maxit.em, thr.em = thr.em,
  #                      maxit.bcd = maxit.bcd, thr.bcd = thr.bcd, conv = out$conv,
  #                      offset = out$offset, covar2corr = covar2corr, truncate = truncate,
  #                      subrout = out$subrout, trace = trace, nobs = n, nresp = p, npred = q)
  class(out.jcglasso) <- "jcglasso"
  out.jcglasso
}

jcglasso.fit <- function (Z, diagonal, weights.B, weights.Tht, penalty, 
                          nrho, rho.min.ratio, rho, nlambda, lambda.min.ratio, lambda, nu, alpha,
                          maxit.em, thr.em, maxit.bcd, thr.bcd, trace, 
                          offset = NULL, covar2corr = FALSE, truncate = 1e-6) {
  
  storage.mode(diagonal) <- "integer"
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
  storage.mode(trace) <- "integer"
  
  if(is.datajcggm(Z)) Z <- list(Z)
  K <- length(Z)
  storage.mode(K) <- "integer"

  n <- sapply(Z, nobs)
  storage.mode(n) <- "integer"
  weights <- n / sum(n)
  storage.mode(weights) <- "double"
  
  p <- sapply(Z, nresp)
  if(all(p == p[1])) { p <- p[1] } else { stop("variabili diverse!??!") }
  storage.mode(p) <- "integer"
  
  q <- sapply(Z, npred)
  if(all(q == q[1])) { q <- q[1] } else { stop("variabili diverse!??!") }
  storage.mode(q) <- "integer"
  
  k_lab <- paste0("class", seq_len(K))
  rho_lab <- paste0("rho_", seq_len(nrho))
  lambda_lab <- paste0("lambda_", seq_len(nlambda))
  
  ynames <- dimnames(Z[[1]])$Y[[2]]
  # yrnames <- dimnames(Z[[1]])$Y[[1]]
  xnames <- dimnames(Z[[1]])$X[[2]]
  znames <- c(ynames, xnames)
  
  X.null <- q == 0L
  
  seq_k <- seq_len(K)
  id_Y <- seq_len(p)
  id_X <- if(!X.null) seq_len(q) + p else NULL
  dim_Z <- p + q
  id_Z <- seq_len(dim_Z)
  storage.mode(dim_Z) <- "integer"
  
  ##### computing output and working objects #####
  Zmat <- vector(mode = "list", length = K)
  zm <- zv <- vector(mode = "list", length = K)
  loz <- upz <- vector(mode = "list", length = K)
  Id <- InfoP <- nP <- vector(mode = "list", length = K)
  Zipt <- vector(mode = "list", length = K)
  B <- mu <- R <- S <- Sgm <- Tht <- vector(mode = "list", length = K)
  B_n <- mu_n <- R_n <- S_n <- Sgm_n <- Tht_n <- vector(mode = "list", length = K)
  Adj <- dfB <- dfTht <- ncomp <- Ck <- pk <- vector(mode = "list", length = K)
  Zipt_lo <- Zipt_up <- Zipt_n <- vector(mode = "list", length = K)
  T1o <- T2o <- vector(mode = "list", length = K)
  T2 <- vector(mode = "list", length = K)
  Thtxx <- vector(mode = "list", length = K)
  for(k in seq_k){
    Zmat[[k]] <- cbind(getMatrix2(Z[[k]], name = "Y", ordered = TRUE), 
                       getMatrix2(Z[[k]], name = "X", ordered = TRUE))
    Zmat[[k]][is.na(Zmat[[k]])] <- 0
    storage.mode(Zmat[[k]]) <- "double"
    
    zm[[k]] <- unlist(ColMeans2(Z[[k]]), use.names = FALSE)
    storage.mode(zm[[k]]) <- "double"
    
    zv[[k]] <- unlist(ColVars2(Z[[k]]), use.names = FALSE)
    storage.mode(zv[[k]]) <- "double"
    
    loz[[k]] <- Z[[k]]$Info$lo
    storage.mode(loz[[k]]) <- "double"
    
    upz[[k]] <- Z[[k]]$Info$up
    storage.mode(upz[[k]]) <- "double"
    
    Id[[k]] <- event2(Z[[k]], ordered = TRUE)
    storage.mode(Id[[k]]) <- "integer"
    
    InfoP[[k]] <- Z[[k]]$Info$Pattern
    storage.mode(InfoP[[k]]) <- "integer"
    
    nP[[k]] <- dim(InfoP[[k]])[1L]
    storage.mode(nP[[k]]) <- "integer"
    
    storage.mode(weights.Tht[[k]]) <- "double"
    
    Zipt[[k]] <- array(0, dim = c(n[[k]], dim_Z, nlambda, nrho), 
                       dimnames = list(NULL, response = znames, lambda = lambda_lab, rho = rho_lab))
    storage.mode(Zipt[[k]]) <- "double"
    
    B[[k]] <- array(0, dim = c(q + 1L, p, nlambda, nrho), 
                    dimnames = list(coef = c("Int.", xnames), response = ynames, lambda = lambda_lab, rho = rho_lab))
    storage.mode(B[[k]]) <- "double"
    B_n[[k]] <- array(0, dim = c(q + 1L, p), 
                      dimnames = list(coef = c("Int.", xnames), response = ynames))
    storage.mode(B_n[[k]]) <- "double"
    
    mu[[k]] <- array(0, dim = c(n[[k]], dim_Z, nlambda, nrho), 
                     dimnames = list(NULL, response = znames, lambda = lambda_lab, rho = rho_lab))
    storage.mode(mu[[k]]) <- "double"
    mu_n[[k]] <- mu[[k]][, , 1L, 1L]
    storage.mode(mu_n[[k]]) <- "double"
    
    R[[k]] <- array(0, dim = c(n[[k]], p, nlambda, nrho), 
                    dimnames = list(NULL, response = ynames, lambda = lambda_lab, rho = rho_lab))
    storage.mode(R[[k]]) <- "double"
    R_n[[k]] <- R[[k]][, , 1L, 1L]
    storage.mode(R_n[[k]]) <- "double"
    
    S[[k]] <- array(0, dim = c(dim_Z, dim_Z, nlambda, nrho), 
                    dimnames = list(response = znames, response = znames, lambda = lambda_lab, rho = rho_lab))
    storage.mode(S[[k]]) <- "double"
    S_n[[k]] <- S[[k]][, , 1L, 1L]
    storage.mode(S_n[[k]]) <- "double"
    
    Sgm[[k]] <- S[[k]]
    storage.mode(Sgm[[k]]) <- "double"
    Sgm_n[[k]] <- Sgm[[k]][, , 1L, 1L]
    storage.mode(Sgm_n[[k]]) <- "double"
    
    Tht[[k]] <- S[[k]]
    storage.mode(Tht[[k]]) <- "double"
    Tht_n[[k]] <- Tht[[k]][, , 1L, 1L]
    storage.mode(Tht_n[[k]]) <- "double"
    
    Adj[[k]] <- S[[k]]
    storage.mode(Adj[[k]]) <- "integer"
    
    dfB[[k]] <- array(0L, dim = c(p + 1L, nlambda, nrho), 
                      dimnames = list(df = c(ynames, "Tot."), lambda = lambda_lab, rho = rho_lab))
    storage.mode(dfB[[k]]) <- "integer"
    
    dfTht[[k]] <- matrix(0L, nrow = nlambda, ncol = nrho, 
                         dimnames = list(lambda = lambda_lab, rho = rho_lab))
    storage.mode(dfTht[[k]]) <- "integer"
    
    ncomp[[k]] <- dfTht[[k]]
    storage.mode(ncomp[[k]]) <- "integer"
    
    Ck[[k]] <- array(0L, dim = c(dim_Z, nlambda, nrho), 
                     dimnames = list(NULL, lambda = lambda_lab, rho = rho_lab))
    storage.mode(Ck[[k]]) <- "integer"
    
    pk[[k]] <- Ck[[k]]
    storage.mode(pk[[k]]) <- "integer"
    
    Zipt_lo[[k]] <- double(dim_Z)
    storage.mode(Zipt_lo[[k]]) <- "double"
    
    Zipt_up[[k]] <- Zipt_lo[[k]]
    storage.mode(Zipt_up[[k]]) <- "double"
    
    Zipt_n[[k]] <- Zipt[[k]][, , 1L, 1L]
    storage.mode(Zipt_n[[k]]) <- "double"
    
    T1o[[k]] <- Zipt_lo[[k]]
    storage.mode(T1o[[k]]) <- "double"
    
    T2o[[k]] <- S_n[[k]]
    storage.mode(T2o[[k]]) <- "double"
    
    T2[[k]] <- S_n[[k]]
    storage.mode(T2[[k]]) <- "double"
    
    if(!X.null) 
      Thtxx[[k]] <- array(0, dim = c(q, q, nlambda, nrho), 
                          dimnames = list(response = xnames, response = xnames, lambda = lambda_lab, rho = rho_lab))
  }

  if(!X.null) {
    wThtxx <- attr(weights.Tht, "Thtxx")
    if(is.null(wThtxx)) {
      wThtxx <- vector(mode = "list", length = K)
      for(k in seq_k) {
        wThtxx[[k]] <- matrix(1.0, nrow = q, ncol = q)
        diag(wThtxx[[k]]) <- 0
        dimnames(wThtxx[[k]]) <- list(xnames, xnames)
        storage.mode(wThtxx[[k]]) <- "double"
      }
    }
  }
  
  nit.tot <- array(0L, dim = c(2L, nlambda, nrho), 
                   dimnames = list(steps = c("EM", "nit"), lambda = lambda_lab, rho = rho_lab))
  storage.mode(nit.tot) <- "integer"
  
  cnnctd <- array(0L, dim = c(dim_Z, nlambda, nrho), 
                  dimnames = list(response = znames, lambda = lambda_lab, rho = rho_lab))
  conv <- integer(1)
  subrout <- integer(1)
  
  is.null.offset <- FALSE
  m_offset <- double(K)
  if(is.null(offset)) {
    is.null.offset <- TRUE
    offset <- sapply(n, function(.n) {
      temp <- rep(0, .n)
      storage.mode(temp) <- "double"
      temp
    }, simplify = FALSE)
  }
  
  for(k in seq_k) {
    ##### computing starting values #####
    T1o[[k]] <- sapply(id_Z, function(j) { id <- Id[[k]][, j] == 0; sum(Zmat[[k]][id, j]) })
    for(i in id_Z) {
      for(j in i:dim_Z) {
        if(any(id <- Id[[k]][, i] == 0 & Id[[k]][, j] == 0)) {
          T2o[[k]][i, j] <- sum(Zmat[[k]][id, i] * Zmat[[k]][id, j])
          T2o[[k]][j, i] <- T2o[[k]][i, j]
        }
      }
    }
    
    Zipt_lo[[k]] <- zm[[k]] - 3 * sqrt(zv[[k]])
    Zipt_up[[k]] <- zm[[k]] + 3 * sqrt(zv[[k]])
    Zipt_n[[k]] <- Zmat[[k]]
    
    ##### computing Yipt_n, Yipt_lo and Yipt_up #####
    for(j in id_Z) {
      if(any(id <- Id[[k]][, j] == 9)) Zipt_n[[k]][id, j] <- zm[[k]][j]
      if(any(id <- Id[[k]][, j] == -1)) {
        z <- (loz[[k]][j] - zm[[k]][j]) / sqrt(zv[[k]][j])
        tmean <- zm[[k]][j] - sqrt(zv[[k]][j]) * dnorm(z) / pnorm(z)
        tmean <- ifelse(tmean <= Zipt_lo[[k]][j], Zipt_lo[[k]][j], tmean)
        Zipt_n[[k]][id, j] <- tmean
      }
      if(any(id <- Id[[k]][, j] == 1)) {
        z <- (upz[[k]][j] - zm[[k]][j]) / sqrt(zv[[k]][j])
        tmean <- zm[[k]][j] + sqrt(zv[[k]][j]) * dnorm(z) / pnorm(z)
        tmean <- ifelse(tmean >= Zipt_up[[k]][j], Zipt_up[[k]][j], tmean)
        Zipt_n[[k]][id, j] <- tmean
      }
    }
    
    zm[[k]] <- colMeans(Zipt_n[[k]])
    m_offset[k] <- mean(offset[[k]])
    B_n[[k]][1L, ] <- zm[[k]][id_Y] - m_offset[k]
    
    ##### computing R_n and xtr_n #####
    mu_n[[k]][, id_Y] <- outer(offset[[k]], B_n[[k]][1L, ], FUN = "+")
    R_n[[k]] <- Zipt_n[[k]][, id_Y] - mu_n[[k]][, id_Y]
    
    ##### identify connected components and computing S #####
    S_n[[k]][id_Y, id_Y] <- cov(R_n[[k]]) * (n[k] - 1) / n[k]
    diag(S_n[[k]]) <- pmax(1E-6, diag(S_n[[k]]))
    if(!X.null) {
      S_n[[k]][id_X, id_X] <- crossprod(Zipt_n[[k]][, id_X]) / n[k]
      S_n[[k]][id_X, id_Y] <- crossprod(Zipt_n[[k]][, id_X], R_n[[k]]) / n[k]
      S_n[[k]][id_Y, id_X] <- t(S_n[[k]][id_X, id_Y]) #xtr_n[[k]]
      mu_n[[k]][, id_X] <- rep(zm[[k]][id_X], each = n[k])
    }
    if(covar2corr) S_n[[k]] <- cov2cor(S_n[[k]]) 
    Sgm_n[[k]] <- diag(diag(S_n[[k]]))
    Tht_n[[k]] <- diag(1 / diag(S_n[[k]]))
  }
  
  U <- outer(1:p, 1:p, "<")
  # U1 <- outer(1:p, 1:p, "<=")
  if(!X.null) U2 <- outer(1:q, 1:q, "<")
  
  if(all(rho == 0)) {
    if(penalty == "group"){
      rho_max <- max(sqrt(rowSums(sapply(seq_k, function(k) pmax(abs(weights[k]*S_n[[k]][id_Y, id_Y][U]) - min(rho), 0))^2)))
    }else{
      rho_max <- max(sapply(seq_k, function(k) max(abs(weights[k]*S_n[[k]][id_Y, id_Y][U]) - min(rho))))
    }
    rho_min <- rho_max * rho.min.ratio
    rho <- exp(seq(from = log(rho_max), to = log(rho_min), length = nrho))
  }
  rho2 <- (1 - alpha) * rho
  rho <- alpha * rho
  if(all(lambda == 0) & !X.null) {
    lambda_max <- max(sqrt(rowSums(sapply(seq_k, function(k) pmax(weights[k]*abs(S_n[[k]][id_X, id_Y]) - min(lambda), 0))^2)))
    lambda_min <- lambda_max * lambda.min.ratio
    lambda <- exp(seq(from = log(lambda_max), to = log(lambda_min), length = nlambda))
  }
  lambda2 <- (1 - alpha) * lambda
  lambda <- alpha * lambda
  if(!X.null) {
    rho.xx <- if(is.null(nu)) rho / alpha else rep(nu, nrho)
    rho2.xx <- (1 - alpha) * rho.xx
    rho.xx <- alpha * rho.xx
  }
  
  for(h in seq_len(nlambda)) {
    for(o in seq_len(nrho)) {
      nit <- double(2L)
      
      rho1mat <- vector(mode = "list", length = K)
      for(k in seq_k){
        rho1mat[[k]] <- rho[o]*weights.Tht[[k]]
        if(!X.null) rho1mat[[k]] <- as.matrix(Matrix:::bdiag(rho1mat[[k]], rho.xx[o]*wThtxx[[k]]))
      }
      rho2mat <- rho2[o] * matrix(1, nrow = p, ncol = p)
      if(!X.null) rho2mat <- as.matrix(Matrix:::bdiag(rho2mat, rho2.xx[o] * matrix(1, nrow = q, ncol = q)))
      if(penalty != "fused") diag(rho2mat) <- 0
      
      T1 <- double(dim_Z)
      
      ##### starting EM algorithm #####
      if(trace == 2) {
        if(!X.null){
          cat("\n*************************************************************************\n",
              "Fitting jcglasso model number ", 
              formatC(o + nrho * (h - 1), digits = 0, width = 5, format = "d"),
              "\n\t\t       nu = ", 
              formatC(rho.xx[o]+rho2.xx[o], digits = 6, width = 9, format = "f"),
              "\n\t\t      rho = ",
              formatC(rho[o]+rho2[o], digits = 6, width = 9, format = "f"),
              "\n\t\t   lambda = ",
              formatC(lambda[h]+lambda2[h], digits = 6, width = 9, format = "f"),
              "\n\t\t    alpha = ",
              formatC(alpha, digits = 3, width = 8, format = "f"),
              "\n")
        } 
        else {
          cat("\n*************************************************************************\n",
              "Fitting jcglasso model number ", 
              formatC(o + nrho * (h - 1), digits = 0, width = 5, format = "d"),
              "\n\t\t      rho = ",
              formatC(rho[o]+rho2[o], digits = 6, width = 9, format = "f"),
              "\n\t\t    alpha = ",
              formatC(alpha, digits = 3, width = 6, format = "f"),
              "\n")
        }
      }
      for(ii in 1:maxit.em) {
        B_o <- B_n
        Tht_o <- Tht_n
        
        ##### computing E step and Multilasso step #####
        for(k in seq_k) {
          temp <- .Fortran(cglasso:::C_e_step_v2, n = n[k], p = dim_Z, Y = Zmat[[k]], 
                           lo = loz[[k]], up = upz[[k]], nP = nP[[k]], 
                           InfoP = InfoP[[k]], T1o = T1o[[k]], T2o = T2o[[k]], 
                           mu_n = mu_n[[k]], Sgm = Sgm_n[[k]], Tht = Tht_n[[k]], 
                           Yipt_lo = Zipt_lo[[k]], Yipt_up = Zipt_up[[k]], 
                           Yipt_n = Zipt_n[[k]], T1 = T1, T2 = T2[[k]], conv = conv)
          if(temp$conv != 0) { subrout <- 1; stop("error in E step!") }
          if(trace == 2) cat("\n\tE-step for group ",
                             formatC(k, digits = 0, width = 3, format = "d"), 
                             ": completed!")
          T2[[k]] <- temp$T2
          zm[[k]] <- temp$T1 / n[k]
          Zipt_n[[k]] <- temp$Yipt_n
          B_n[[k]][1L, ] <- zm[[k]][id_Y] - m_offset[[k]]
          mu_n[[k]][, id_Y] <- outer(offset[[k]], B_n[[k]][1L, ], FUN = "+")
          R_n[[k]] <- Zipt_n[[k]][, id_Y] - mu_n[[k]][, id_Y]
          if(!X.null) {
            S_n[[k]][id_X, id_X] <- T2[[k]][id_X, id_X] / n[k]
            S_n[[k]][id_X, id_Y] <- crossprod(Zipt_n[[k]][, id_X], R_n[[k]]) / n[k]
            S_n[[k]][id_Y, id_X] <- t(S_n[[k]][id_X, id_Y]) #xtr_n[[k]]
            mu_n[[k]][, id_X] <- rep(zm[[k]][id_X], each = n[k])
            
            # Sgm_n[[k]][id_X, id_X] <- T2[[k]][id_X, id_X] / n[k]
            # Thtxx_n[[k]] <- solve(Sgmxx_n[[k]]) #diag(1 / diag(Sxx_n[[k]]))
            # xtr_n[[k]] <- crossprod(Xipt_n[[k]], R_n[[k]]) / n[k]
            # mu_n[[k]] <- matrix(rep(xm[[k]], each = n[k]), nrow = n[k], ncol = q)
          }
        }
        
        ##### fitting multilasso model #####
        if(!X.null) {
          if(trace == 2) cat("\n\tM-step:\n\t       fitting multivariate sparse group lasso model with ",
                             "\n\t       lambda = ", formatC(lambda[h]+lambda2[h], digits = 6, width = 9, format = "f"), 
                             "and alpha = ", formatC(alpha, digits = 3, width = 4, format = "f"), 
                             "\n\n")
          # tempmul <- admm.iters.3(p, q, ym, xm, xtx_n, xtr_n, B_n, Tht_n, weights.B, lambda[h], lambda2[h], maxit.bcd, 
          #                         thr.bcd, rho = 2, rho.increment = 1, 0L)
          
          tempmul <- admm.iters.3(p, q, 
                                  ym = lapply(zm, function(x) x[id_Y]), 
                                  xm = lapply(zm, function(x) x[id_X]), 
                                  xtx = lapply(S_n, function(x) x[id_X, id_X]), 
                                  xtr = lapply(S_n, function(x) x[id_X, id_Y]), B = B_n, 
                                  Tht = lapply(Tht_n, function(x) x[id_Y, id_Y]), wB = weights.B, 
                                  lambda[h], lambda2[h], maxit.bcd, thr.bcd, weights, 
                                  rho = 2, rho.increment = 1.0, trace = 0L)
          
          if(tempmul$conv != 0) {
            conv <- tempmul$conv
            subrout <- 2
            stop("error in multilasso step!")
          }
          if(trace == 2) {
            cat("\t\tADMM step\t||B_new(k) - B_old(k)||_1/||B_old(k)||_1\n")
            cat("\t\t   ", formatC(tempmul$nit[1], digits = 0, width = 5, format = "d"),
                "\t\t\t", formatC(tempmul$diff, digits = 8, width = 10, format = "f"),"\n")
          }
          for(k in seq_k){
            B_n[[k]] <- tempmul$B[[k]]
            dfB[[k]][, h, o] <- tempmul$df[[k]]
            nit[2] <- nit[2] + tempmul$nit[1]
            
            mu_n[[k]][, id_Y] <- outer(offset[[k]], B_n[[k]][1L, ], FUN = "+")
            mu_n[[k]][, id_Y] <- mu_n[[k]][, id_Y] + Zipt_n[[k]][, id_X] %*% B_n[[k]][-1L, , drop = FALSE]
            R_n[[k]] <- Zipt_n[[k]][, id_Y] - mu_n[[k]][, id_Y]
            YM <- crossprod(Zipt_n[[k]][, id_Y], mu_n[[k]][, id_Y])
            # temp$T2 == crossprod(Yipt_n[[k]])
            S_n[[k]][id_Y, id_Y] <- (T2[[k]][id_Y, id_Y] + crossprod(mu_n[[k]][, id_Y]) - YM - t(YM)) / n[k]
            if(covar2corr) S_n[[k]][id_Y, id_Y] <- cov2cor(S_n[[k]][id_Y, id_Y]) #S_n[[k]] / outer(sqrt(diag(S_n[[k]])), sqrt(diag(S_n[[k]])))
            # Sgm_n[[k]] <- S_n[[k]] + rho1mat[[k]] * sign(Tht_n[[k]])
          }
        }
        else{
          for(k in seq_k){
            YM <- crossprod(Zipt_n[[k]][, id_Y], mu_n[[k]][, id_Y])
            # temp$T2 == crossprod(Yipt_n[[k]])
            # T2 <- crossprod(Yipt_n[[k]])
            S_n[[k]][id_Y, id_Y] <- (T2[[k]][id_Y, id_Y] + crossprod(mu_n[[k]][, id_Y]) - YM - t(YM)) / n[k]
            if(covar2corr) S_n[[k]][id_Y, id_Y] <- cov2cor(S_n[[k]][id_Y, id_Y])
          }
        }
        
        ##### checking connected components #####
        
        if(penalty == "fused") {
          crit1 <- vector(mode = "list", length = K)
          if(K == 2) {
            S.sum <- matrix(0, p, p)
            for(k in seq_k) { 
              temp1 <- weights[k] * S_n[[k]][id_Y, id_Y]
              S.sum <- S.sum + temp1
              crit1[[k]] <- abs(temp1) > rho[o] + rho2[o]
            }  
            crit2 <- (abs(S.sum) > 2 * rho[o])
          }
          else {
            for(k in seq_k) { crit1[[k]] <- abs(weights[k] * S_n[[k]][id_Y, id_Y]) > rho[o] }  
            crit2 <- 0
          }
          
          critboth <- crit2
          for(k in seq_k) { critboth <- critboth + crit1[[k]] }
          critboth <- (critboth != 0)				
          diag(critboth) <- 1
        } 
        else {
          tempsum <- 0
          for(k in seq_k) { tempsum <- tempsum + (pmax(abs(weights[k] * S_n[[k]][id_Y, id_Y]) - rho[o], 0))^2 }    
          critboth <- (tempsum > rho2[o]^2)
          diag(critboth) <- 1
          if(!X.null){
            tempsum2 <- 0
            for(k in seq_k) { tempsum2 <- tempsum2 + (pmax(abs(weights[k] * S_n[[k]][id_X, id_X]) - rho.xx[o], 0))^2 }    
            critboth2 <- (tempsum2 > rho2.xx[o]^2)
            diag(critboth2) <- 1
            critboth <- as.matrix(Matrix:::bdiag(critboth, critboth2))
          }
        }
        
        ##### identify block structure using igraph #####
        g1 <- graph.adjacency(critboth)	
        cout <- clusters(g1)
        
        ##### identify unconnected elements, and get blocks #####
        if(min(cout$membership) == 0){ cout$membership <- cout$membership + 1 }
        unconnected <- cout$membership %in% setdiff(seq_len(cout$no), which(cout$csize > 1))
        connected <- (!unconnected)
        blocklist <- sapply(which(cout$csize > 1), function(.bl) which(cout$membership == .bl), simplify = FALSE)
        
        whole.S <- lapply(S, function(x) x[, , h, o])
        whole.theta <- lapply(Tht, function(x) x[, , h, o])
        
        ##### computing theta and S values for unconnected components #####
        if(sum(unconnected) != 0) {
          if(covar2corr) {
            for(k in seq_k) {
              diag(whole.theta[[k]])[unconnected] <- rep(1, sum(unconnected))
              diag(whole.S[[k]])[unconnected] <- rep(1, sum(unconnected))
            }
          } 
          else {
            S.uncon <- lapply(S_n, function(.S_n) diag(.S_n)[unconnected])
            for(k in seq_k) {
              diag(whole.theta[[k]])[unconnected] <- 1 / S.uncon[[k]]
              diag(whole.S[[k]])[unconnected] <- S.uncon[[k]]
            }
            rm(list = c("S.uncon"))
          }
        }
        
        ##### computing theta and S values for connected components #####
        if(trace == 2) {
          if(X.null) {
            cat("\n")
            cat("\tM-step:\n\t       fitting joint glasso model with",
                " rho = ", formatC(rho[o]+rho2[o], digits = 6, width = 9, format = "f"),
                "and alpha = ", formatC(alpha, digits = 3, width = 4, format = "f"),
                "\n\n")
          } 
          else {
            cat("\tM-step:\n\t       fitting joint glasso model with",
                "\n\t       nu = ", formatC(rho.xx[o]+rho2.xx[o], digits = 6, width = 9, format = "f"), 
                ", rho = ", formatC(rho[o]+rho2[o], digits = 6, width = 9, format = "f"),
                "and alpha = ", formatC(alpha, digits = 3, width = 4, format = "f"),
                "\n\n")
          }
        }
        if(length(blocklist) > 0){
          for(i in seq_len(length(blocklist))) {
            
            # the variables in the block
            bl <- blocklist[[i]] 
            
            # run the admm algorithm on the block:
            Thetabl <- admm.iters.2(S = lapply(S_n, function(.S_n) .S_n[bl, bl]), 
                                    theta = lapply(Tht_n, function(.Tht_n) .Tht_n[bl, bl]), 
                                    lam1 = lapply(rho1mat, function(x) x[bl, bl]), 
                                    lam2 = rho2mat[bl, bl], 
                                    penalty = penalty, weights = weights, rho = 1, rho.increment = 1.0,
                                    penalize.diagonal = diagonal, maxiter = maxit.bcd, 
                                    tol = thr.bcd, covar2corr = covar2corr, trace = 0L)
            if(trace == 2) {
              cat("\tconnected components number ", 
                  formatC(i, digits = 0, width = 4, format = "d"), "\n")
              cat("\t\tADMM step\t||Tht_new(k) - Tht_old(k)||_1/||Tht_old(k)||_1\n")
              cat("\t\t   ", formatC(Thetabl$iter, digits = 0, width = 5, format = "d"),
                  "\t\t\t", formatC(Thetabl$diff, digits = 8, width = 10, format = "f"),"\n")
            }
            if(Thetabl$conv) {
              conv <- 1
              subrout <- 2
              stop("error in M step (connected components)!")
            }
            nit[2] <- nit[2] + Thetabl$iters
            
            # update results:
            for(k in seq_k) { 
              Thetabl$Z[[k]][abs(Thetabl$Z[[k]]) <= truncate] <- 0
              whole.theta[[k]][bl, bl] <- Thetabl$Z[[k]]
              whole.S[[k]][bl, bl] <- Thetabl$S[[k]]
            }   
          }
          rm(list = c("Thetabl"))
        }
        Sgm_n <- whole.S
        Tht_n <- whole.theta
        nit[1] <- ii
        
        # dB <- sum(sapply(seq_len(K), function(k) sqrt(sum((B_n[[k]] - B_o[[k]])^2)) / (p + dfB[[k]][p + 1, h, o])))
        # dSgm <- sum(sapply(seq_len(K), function(k) sqrt(sum((Tht_n[[k]] - Tht_o[[k]])[id_Y, id_Y][U1]^2)) / sum(Tht_n[[k]][id_Y, id_Y][U1] != 0)))
        
        # dB <- sum(sapply(seq_len(K), function(k) sum(abs(B_n[[k]] - B_o[[k]])) / sum(abs(B_o[[k]]))))
        # dTht <- sum(sapply(seq_len(K), function(k) sum(abs((Tht_n[[k]][id_Y, id_Y] - Tht_o[[k]][id_Y, id_Y]))) / sum(abs(Tht_o[[k]][id_Y, id_Y]))))

        dB <- sum(sapply(seq_len(K), function(k) norm(B_n[[k]] - B_o[[k]], type = "F") / (p + dfB[[k]][p + 1, h, o])))
        dTht <- sum(sapply(seq_len(K), function(k) norm((Tht_n[[k]][id_Y, id_Y] - Tht_o[[k]][id_Y, id_Y]), type = "F") / (sum(Tht_n[[k]][id_Y, id_Y][U] != 0) + p)))

        if(trace == 2) cat("\tM-step: completed!\n")
        
        if(!X.null){
          for(k in seq_k) {
            Adj[[k]][, , h, o] <- 1*(Tht_n[[k]] != 0)
            Thtxx[[k]][, , h, o] <- Tht_n[[k]][id_X, id_X]
            Tht_n[[k]][id_X, id_X] <- Tht_n[[k]][id_X, id_X] + B_n[[k]][-1L, ] %*% Tht_n[[k]][id_Y, id_Y] %*% t(B_n[[k]][-1L, ])
            Tht_n[[k]][id_X, id_Y] <- -(B_n[[k]][-1L, ] %*% Tht_n[[k]][id_Y, id_Y])
            Tht_n[[k]][id_Y, id_X] <- t(Tht_n[[k]][id_X, id_Y])
          }
        } 
        else {
          for(k in seq_k) Adj[[k]][, , h, o] <- 1*(Tht_n[[k]] != 0)
        }
        
        if(trace == 2) {
          cat("\n\tChecking convergence criterion (threshold = ", 
              formatC(thr.em, digits = 6, width = 8, format = "f"),
              ")\n\t||B_old - B_new|_F / dfB  = ",
              formatC(dB, digits = 7, width = 9, format = "f"),
              "\n\t||Tht_old - Tht_new||_F / dfTht = ",
              formatC(dTht, digits = 7, width = 9, format = "f"), "\n")
          if(max(dB, dTht) <= thr.em) cat("\tConvergence criterion is met!\n\n")
        } 
        if(max(dB, dTht) <= thr.em) break
      }
      if(ii >= maxit.em) {
        conv <- 1
        subrout <- 3
        stop("error in M step (maximum number of iteration attained)!")
      }
      if(trace == 1) {
        if(!X.null){
          cat("\njcglasso model number ", 
              formatC(o + nrho * (h - 1), digits = 0, width = 5, format = "d"),
              "\n\t         nu = ",
              formatC(rho.xx[o]+rho2.xx[o], digits = 6, width = 9, format = "f"),
              "\n\t        rho = ",
              formatC(rho[o]+rho2[o], digits = 6, width = 9, format = "f"),
              "\n\t     lambda = ",
              formatC(lambda[h]+lambda2[h], digits = 6, width = 9, format = "f"),
              "\n\t      alpha = ",
              formatC(alpha, digits = 3, width = 6, format = "f"),
              "\n\t   EM-steps = ",
              formatC(nit[1], digits = 0, width = 5, format = "d"),
              "\n\t      steps = ",
              formatC(nit[2], digits = 0, width = 5, format = "d"),
              "\n"
          )
        } 
        else {
          cat("\njcglasso model number ", 
              formatC(o + nrho * (h - 1), digits = 0, width = 5, format = "d"),
              "\n\t        rho = ",
              formatC(rho[o]+rho2[o], digits = 6, width = 9, format = "f"),
              "\n\t      alpha = ",
              formatC(alpha, digits = 3, width = 6, format = "f"),
              "\n\t   EM-steps = ",
              formatC(nit[1], digits = 0, width = 5, format = "d"),
              "\n\t      steps = ",
              formatC(nit[2], digits = 0, width = 5, format = "d"),
              "\n"
          )
        }
      }
      
      ##### buiding output components #####
      for(k in seq_k) {
        Zipt[[k]][, , h, o] <- Zipt_n[[k]]
        B[[k]][, , h, o] <- B_n[[k]]
        mu[[k]][, , h, o] <- mu_n[[k]]
        R[[k]][, , h, o] <- R_n[[k]]
        S[[k]][, , h, o] <- S_n[[k]]
        Sgm[[k]][, , h, o] <- Sgm_n[[k]]
        Tht[[k]][, , h, o] <- Tht_n[[k]]
        ncomp[[k]][h, o] <- cout$no
        Ck[[k]][, h, o] <- cout$membership
        pk[[k]][1:length(cout$csize), h, o] <- cout$csize
        # Adj[[k]][, , h, o] <- 1*(Tht[[k]][, , h, o] != 0)
        diag(Adj[[k]][, , h, o]) <- 0
        # dimnames(Adj[[k]])[1:2] <- list(ynames, ynames)
        dfTht[[k]][h, o] <- sum(Adj[[k]][id_Y, id_Y, h, o][U])
        if(!X.null) {
          Adj[[k]][id_X, id_Y, h, o] <- 1*(B[[k]][-1, , h, o] != 0)
          Adj[[k]][id_Y, id_X, h, o] <- t(Adj[[k]][id_X, id_Y, h, o])
          dfTht[[k]][h, o] <- dfTht[[k]][h, o] + sum(Adj[[k]][id_X, id_X, h, o][U2])
        }
        cnnctd[, h, o] <- connected
        
        row.order <- order(Z[[k]]$Info$order)
        Zipt[[k]][, , h, o] <- Zipt[[k]][, , h, o][row.order, , drop = FALSE]
        mu[[k]][, , h, o] <- mu[[k]][, , h, o][row.order, , drop = FALSE]
        R[[k]][, , h, o] <- R[[k]][, , h, o][row.order, , drop = FALSE]
      }
      nit.tot[, h, o] <- nit
    }
  }
  
  names(Zipt) <- names(B) <- names(mu) <- names(R) <- k_lab
  names(S) <- names(Sgm) <- names(Tht) <- k_lab
  names(ncomp) <- names(Ck) <- names(pk) <- k_lab
  names(Adj) <- names(dfTht) <- names(dfB) <- k_lab
  names(Z) <- names(Thtxx) <- k_lab
  
  if(!X.null){
    for(k in seq_k) weights.Tht[[k]] <- as.matrix(Matrix:::bdiag(weights.Tht[[k]], wThtxx[[k]]))
  }
  
  out <- list(n = n, p = p, q = q, id_X = id_X, id_Y = id_Y, Z = Zmat, 
              Id = Id, nP = nP, InfoP = InfoP, ncomp = ncomp, Ck = Ck, pk = pk, 
              lo = loz, up = upz, zm = zm, zv = zv, wTht = weights.Tht, 
              nrho = nrho, rhoratio = rho.min.ratio, rho = rho,
              nlambda = nlambda, lambdaratio = lambda.min.ratio, lambda = lambda, 
              alpha = alpha, rho2 = rho2, lambda2 = lambda2, 
              pendiag = diagonal, connected = cnnctd, 
              maxit_em = maxit.em, thr_em = thr.em, maxit_bcd = maxit.bcd, thr_bcd = thr.bcd, 
              Zipt = Zipt, B = B, mu = mu, R = R, S = S, Sgm = Sgm, Tht = Tht, 
              Adj = Adj, Omega = Thtxx, dfB = dfB, dfTht = dfTht,
              nit = nit.tot, conv = conv, offset = offset,
              subrout = subrout, trace = trace)
        
  # out <- list(n = n, p = p, Y = lapply(Zmat, function(x) x[, id_Y]), 
  #             Id = Id, nP = nP, InfoP = InfoP, 
  #             lo = lapply(loz, function(x) x[id_Y]), 
  #             up = lapply(upz, function(x) x[id_Y]), 
  #             pendiag = diagonal, wTht = weights.Tht, 
  #             ym = lapply(zm, function(x) x[id_Y]), 
  #             yv = lapply(zv, function(x) x[id_Y]), 
  #             nrho = nrho, rhoratio = rho.min.ratio, rho = rho,
  #             nlambda = nlambda, lambdaratio = lambda.min.ratio, lambda = lambda, 
  #             alpha = alpha,
  #             rho2 = rho2, lambda2 = lambda2,
  #             maxit_em = maxit.em, thr_em = thr.em, maxit_bcd = maxit.bcd, thr_bcd = thr.bcd, 
  #             connected = cnnctd, 
  #             Xipt = if(!X.null) lapply(Zipt, function(x) x[, id_X, , , drop = FALSE]) else NULL, 
  #             Yipt = lapply(Zipt, function(x) x[, id_Y, , , drop = FALSE]), 
  #             B = B, mu = lapply(mu, function(x) x[, id_Y, , , drop = FALSE]), 
  #             R = R, S = lapply(S, function(x) x[id_Y, id_Y, , , drop = FALSE]), 
  #             Sgm = lapply(Sgm, function(x) x[id_Y, id_Y, , , drop = FALSE]), 
  #             Tht = lapply(Tht, function(x) x[id_Y, id_Y, , , drop = FALSE]), 
  #             Sxx = if(!X.null) lapply(S, function(x) x[id_X, id_X, , , drop = FALSE]) else NULL, 
  #             Sxy = if(!X.null) lapply(S, function(x) x[id_X, id_Y, , , drop = FALSE]) else NULL, 
  #             Sgmxx = if(!X.null) lapply(Sgm, function(x) x[id_X, id_X, , , drop = FALSE]) else NULL, 
  #             Sgmxy = if(!X.null) lapply(Sgm, function(x) x[id_X, id_Y, , , drop = FALSE]) else NULL, 
  #             Thtxx = if(!X.null) lapply(Tht, function(x) x[id_X, id_X, , , drop = FALSE]) else NULL, 
  #             Thtxy = if(!X.null) lapply(Tht, function(x) x[id_X, id_Y, , , drop = FALSE]) else NULL,
  #             Adj_yy = lapply(Adj, function(x) x[id_Y, id_Y, , , drop = FALSE]), 
  #             Adj_xy = if(!X.null) lapply(B, function(x) 1*(x[-1L, , , , drop = FALSE] != 0)) else NULL, 
  #             Adj_xx = if(!X.null) lapply(Adj, function(x) 1*(x[id_X, id_X, , , drop = FALSE] != 0)) else NULL, 
  #             dfB = dfB, dfTht = dfTht, ncomp = ncomp, Ck = Ck, pk = pk, 
  #             nit = nit.tot, conv = conv, offset = offset,
  #             subrout = subrout, trace = trace)
  out$conv <- switch(as.character(out$conv),
                     "-1" = "memory allocation error",
                     "0" = "Ok",
                     "1" = "maximum number of iterations has been exceeded",
                     "2" = "error in E-step",
                     "3" = "matrix inversion failed")
  
  return(out)
}

##### jcglasso tools #####
admm.iters.2 <- function (S, theta, lam1, lam2, penalty = "fused", rho = 1, rho.increment = 1.0, 
                          weights, penalize.diagonal, maxiter = 1000, tol = 1e-04, 
                          trace = FALSE, covar2corr = TRUE) {
  K <- length(S)
  p <- dim(S[[1]])[2]
  n <- weights
  
  A <- W <- Z <- vector(mode = "list", length = K)
  for(k in seq_len(K)) A[[k]] <- W[[k]] <- Z[[k]] <- matrix(0, p, p)
  iter <- 1
  diff_value <- 10
  # eps1 <- eps2 <- 10
  # criterion <- TRUE
  while ((iter <= maxiter) && diff_value > tol) {
    if(trace) cat("\n", paste0("iter = ", formatC(iter, digits = 0, width = 5, format = "d")), 
                  paste0(" rho = ", formatC(rho, digits = 2, width = 8, format = "f")), 
                  paste0(" diff = ", formatC(diff_value, digits = 5, width = 8, format = "f")))
    # paste0(" diff = ", formatC(eps1, digits = 5, width = 8, format = "f"), 
    #        " - diff2 = ", formatC(eps2, digits = 5, width = 8, format = "f")))
    theta.prev <- Z
    for (k in seq_len(K)) {
      edecomp <- eigen(S[[k]] - rho * Z[[k]] / n[k] + rho * W[[k]] / n[k])
      D <- edecomp$values
      V <- edecomp$vectors
      D2 <- n[k] / (2 * rho) * (-D + sqrt(D^2 + 4 * rho / n[k]))
      theta[[k]] <- V %*% diag(D2) %*% t(V)
      A[[k]] <- theta[[k]] + W[[k]]
    }
    
    if (penalty == "lasso") {
      Z <- lasso(A, rho, lam1, penalize.diagonal = penalize.diagonal)
    } 
    else {
      if (penalty == "fused") {
        if (K == 2) {
          Z <- flsa2(A, rho, lam1, lam2, penalize.diagonal = penalize.diagonal)
        }
        else {
          Z <- flsa2.general(A, rho, lam1, lam2, penalize.diagonal = penalize.diagonal)
        }
      }
      else {
        Z <- dsgl2(A, rho, lam1, lam2, penalize.diagonal = penalize.diagonal)
      }
    }
    
    for (k in seq_len(K)) W[[k]] <- W[[k]] + (theta[[k]] - Z[[k]])
    iter <- iter + 1
    
    # s <- sum(sapply(seq_len(K), function(k) rho*norm(Z[[k]] - theta.prev[[k]], "F")))
    # r <- sum(sapply(seq_len(K), function(k) rho*norm(theta[[k]] - Z[[k]], "F")))
    # if (r > 10*s){
    #   rho <- rho * 1.2
    # }
    # if (s > 10*r){
    #   rho <- rho / 1.2
    # }
    # eps1 <- p * tol + tol * sum(sapply(seq_len(K), function(k) max(norm(theta[[k]], "F"), norm(Z[[k]], "F"))))
    # eps2 <- p * tol + tol * sum(sapply(seq_len(K), function(k) norm(W[[k]], "F")))
    # criterion <- (r >= eps1 | s >= eps2)
    
    diff_value <- 0
    for (k in seq_len(K)) {
      diff_value <- diff_value + sum(abs(Z[[k]] - theta.prev[[k]])) / sum(abs(theta.prev[[k]]))
    }
    if(iter %% 10 == 0) rho <- rho * rho.increment
  }
  # diff <- 0
  for (k in 1:K) {
    # diff <- diff + sum(abs(theta[[k]] - Z[[k]]))
    S[[k]] <- solve(Z[[k]])
    # id <- (Z[[k]] == 0)
    if(covar2corr) {
      S[[k]] <- cov2cor(S[[k]]) #S[[k]] / outer(sqrt(diag(S[[k]])), sqrt(diag(S[[k]])))
      Z[[k]] <- solve(S[[k]])
    }
    # Z[[k]][id] <- 0
  }
  out <- list(theta = theta, Z = Z, S = S, diff = diff_value, conv = 1*(iter == maxiter), iters = iter)
  return(out)
}

flsa2 <- function (A, L, lam1, lam2, penalize.diagonal) {
  S1 <- abs(A[[1]] - A[[2]]) <= 2 * lam2/L
  X1 <- (A[[1]] + A[[2]])/2
  Y1 <- X1
  S2 <- (A[[1]] > A[[2]] + 2 * lam2/L)
  X2 <- A[[1]] - lam2/L
  Y2 <- A[[2]] + lam2/L
  S3 <- (A[[2]] > A[[1]] + 2 * lam2/L)
  X3 <- A[[1]] + lam2/L
  Y3 <- A[[2]] - lam2/L
  X <- soft(a = S1 * X1 + S2 * X2 + S3 * X3, lam = lam1[[1]]/L, penalize.diagonal = penalize.diagonal)
  Y <- soft(a = S1 * Y1 + S2 * Y2 + S3 * Y3, lam = lam1[[2]]/L, penalize.diagonal = penalize.diagonal)
  return(list(X, Y))
}
flsa2.general <- function (A, L, lam1, lam2, penalize.diagonal) {
  trueA <- A
  K <- length(A)
  if (is.matrix(A[[1]])) { 
    p <- dim(A[[1]])[1] 
    fusions <- array(FALSE, dim = c(K, K, p, p))
  }
  if (is.vector(A[[1]])) { 
    p <- length(A[[1]]) 
    fusions = array(FALSE, dim = c(K, K, p, 1))
  }
  newc <- list()
  for (k in seq_len(K)) {
    newc[[k]] = A[[k]] * 0
    others <- setdiff(1:K, k)
    others.smaller.k <- 1:(k - 1)
    for (o in others) {
      newc[[k]] <- newc[[k]] + (A[[o]] - A[[k]] < -1e-04) - (A[[o]] - A[[k]] > 1e-04)
    }
  }
  for (iter in seq_len(K - 1)) {
    ordermats <- list()
    for (k in seq_len(K)) {
      ordermats[[k]] = A[[k]] * 0
      others <- setdiff(1:K, k)
      others.smaller.k <- 1:(k - 1)
      for (o in others) {
        ordermats[[k]] <- ordermats[[k]] + (A[[k]] - A[[o]] > 1e-04)
      }
      if (k > 1) {
        for (o in others.smaller.k) {
          ordermats[[k]] <- ordermats[[k]] + (abs(A[[o]] - A[[k]]) < 1e-04)
        }
      }
      ordermats[[k]] <- ordermats[[k]] + 1
    }
    betas.g <- A
    for (k in seq_len(K)) {
      betas.g[[k]] = A[[k]] - lam2 / L * newc[[k]]
    }
    new.ordermats <- list()
    for (k in seq_len(K)) {
      new.ordermats[[k]] = A[[k]] * 0
      others <- setdiff(1:K, k)
      others.smaller.k <- 1:(k - 1)
      for (o in others) {
        new.ordermats[[k]] <- new.ordermats[[k]] + (betas.g[[k]] - betas.g[[o]] > 1e-04)
      }
      if (k > 1) {
        for (o in others.smaller.k) {
          new.ordermats[[k]] <- new.ordermats[[k]] + (abs(betas.g[[o]] - betas.g[[k]]) < 1e-04)
        }
      }
      new.ordermats[[k]] <- new.ordermats[[k]] + 1
    }
    for (k in seq_len(K)) {
      for (kp in seq_len(K)) {
        fusions[k, kp, , ] <- fusions[k, kp, , ] + 
          ((ordermats[[k]] - 1 == ordermats[[kp]]) & (new.ordermats[[k]] < new.ordermats[[kp]])) + 
          ((ordermats[[k]] + 1 == ordermats[[kp]]) & (new.ordermats[[k]] > new.ordermats[[kp]])) + 
          (abs(A[[k]] - A[[kp]]) < 1e-04)
        fusions <- (fusions > 0) * 1
      }
    }
    for (k in seq_len(K)) {
      for (kp in seq_len(K)) {
        others <- setdiff(1:K, c(k, kp))
        for (o in others) {
          bothfused <- fusions[k, o, , ] & fusions[kp, o, , ]
          fusions[k, kp, , ] <- fusions[k, kp, , ] | bothfused
        }
      }
    }
    for (k in seq_len(K)) {
      others <- setdiff(1:K, k)
      fusemean <- trueA[[k]]
      denom <- A[[k]] * 0 + 1
      for (o in others) {
        fusemean <- fusemean + fusions[k, o, , ] * trueA[[o]]
        denom <- denom + fusions[k, o, , ]
      }
      A[[k]] <- fusemean / denom
    }
    newc <- list()
    for (k in seq_len(K)) {
      newc[[k]] <- A[[k]] * 0
      others <- setdiff(1:K, k)
      others.smaller.k <- 1:(k - 1)
      for (o in others) {
        newc[[k]] <- newc[[k]] + (A[[o]] - A[[k]] < -1e-04) - (A[[o]] - A[[k]] > 1e-04)
      }
    }
  }
  for (k in seq_len(K)) {
    betas.g[[k]] <- A[[k]] - lam2 / L * newc[[k]]
  }
  X <- list()
  for (k in seq_len(K)) {
    X[[k]] <- soft(betas.g[[k]], lam = lam1[[k]] / L, penalize.diagonal = penalize.diagonal)
  }
  return(X)
}
dsgl2 <- function (A, L, lam1, lam2, penalize.diagonal) {
  
  lam2 <- lam2 * 1/L
  if (is.matrix(A[[1]])) { p <- dim(A[[1]])[1] }
  if (is.vector(A[[1]])) { p <- length(A[[1]]) }
  K <- length(A)
  softA <- A
  normsoftA <- A[[1]] * 0
  for (k in 1:K) {
    softA[[k]] <- soft(A[[k]], lam1[[k]] / L, penalize.diagonal = penalize.diagonal)
    normsoftA <- normsoftA + (softA[[k]])^2
  }
  normsoftA <- sqrt(normsoftA)
  notshrunk <- (normsoftA > lam2) * 1
  normsoftA <- normsoftA + (1 - notshrunk)
  out <- A
  for (k in 1:K) {
    out[[k]] <- softA[[k]] * (1 - lam2 / normsoftA)
    # if(any(id <- is.na(out[[k]]))) out[[k]][id] <- 0
    out[[k]] <- out[[k]] * notshrunk
  }
  return(out)
}
lasso <- function (A, L, lam1, penalize.diagonal) {
  if (is.matrix(A[[1]])) { p <- dim(A[[1]])[1] }
  if (is.vector(A[[1]])) { p <- length(A[[1]]) }
  K <- length(A)
  softA <- A
  normsoftA <- A[[1]] * 0
  for (k in 1:K) {
    softA[[k]] <- soft(A[[k]], lam1[[k]] * 1/L, penalize.diagonal = penalize.diagonal)
  }
  out <- softA
  return(out)
}
soft <- function (a, lam, penalize.diagonal) {
  out <- sign(a) * pmax(0, abs(a) - lam)
  if (!penalize.diagonal) diag(out) <- diag(a)
  return(out)
}

admm.iters.3 <- function (p, q, ym, xm, xtx, xtr, B, Tht, wB, lmb1, lmb2, maxit, thr, weights, rho = 1, rho.increment = 1, trace) {
  
  softl1 <- function(A, l1) pmax(0, abs(A) - l1) * sign(A)
  softl2 <- function(A, normA, l2) {
    A <- pmax(0, 1 - l2 / normA) * A
    if(any(id <- is.na(A))) A[id] <- 0
    A
  }
  
  K <- length(xtx)
  pq <- p * q
  
  u <- x <- vector(mode = "list", length = K)
  XtRTht <- THTXtX <- W <- invTHTXtX <- vector(mode = "list", length = K)
  update <- 1*(rho.increment != 1)
  Iqp <- diag(1, nrow = pq, ncol = pq)
  for(k in seq_len(K)) {
    u[[k]] <- matrix(0, nrow = q, ncol = p)
    x[[k]] <- B[[k]][-1, , drop = FALSE]
    # xtr[[k]] <- xtr[[k]] * n[k]^2 / sum(n)
    # xtx[[k]] <- xtx[[k]] * n[k]^2 / sum(n)
    XtRTht[[k]] <- xtr[[k]] %*% Tht[[k]] # q x p #2 * weights[k] * 
    THTXtX[[k]] <- Tht[[k]] %x% xtx[[k]] # qp x qp #2 * weights[k] * 
    W[[k]] <- wB[[k]] * rep(diag(Tht[[k]]), each = q)
    if(!update) invTHTXtX[[k]] <- chol2inv(chol(THTXtX[[k]] + rho * Iqp))
    # if(!update) invTHTXtX[[k]] <- chol(THTXtX[[k]] + rho * Iqp)
  }
  
  step <- 0
  while ((step <- step + 1) <= maxit) {
    Bold <- B
    eps <- 0
    # browser()
    
    normsoftB <- B[[1]][-1, , drop = FALSE] * 0
    for(k in seq_len(K)){
      if(update) invTHTXtX[[k]] <- chol2inv(chol(THTXtX[[k]] + rho * Iqp))
      # if(update) invTHTXtX[[k]] <- chol(THTXtX[[k]] + rho * Iqp)
      tempY <- c(XtRTht[[k]] + rho * (B[[k]][-1, , drop = FALSE] - u[[k]]))
      
      x[[k]] <- matrix(invTHTXtX[[k]] %*% tempY, nrow = q, ncol = p)
      # x[[k]] <- matrix(forwardsolve(invTHTXtX[[k]],
      #                               backsolve(invTHTXtX[[k]], tempY, transpose = TRUE),
      #                               upper.tri = TRUE),
      #                  nrow = q, ncol = p)
      B[[k]][-1, ] <- softl1(u[[k]] + x[[k]], lmb1 * W[[k]] / rho)
      normsoftB <- normsoftB + (B[[k]][-1, , drop = FALSE])^2
      B[[k]][1, ] <- ym[[k]] - crossprod(B[[k]][-1, , drop = FALSE], xm[[k]])
      u[[k]] <- u[[k]] + (x[[k]] - B[[k]][-1, , drop = FALSE])
    }
    normsoftB <- sqrt(normsoftB)
    notshrunk <- (normsoftB > (lmb2 / rho)) * 1
    normsoftB <- normsoftB + (1 - notshrunk)
    normsoftB <- (1 - (lmb2 / rho) / normsoftB) * notshrunk
    for (k in 1:K) {
      B[[k]][-1, ] <- B[[k]][-1, , drop = FALSE] * normsoftB
      # B[[k]][-1, ] <- softl2(B[[k]][-1, , drop = FALSE], normsoftB, lmb2 / rho)
      # u[[k]] <- u[[k]] + x[[k]] - B[[k]][-1, , drop = FALSE]
      eps <- eps + sum(abs(B[[k]] - Bold[[k]])) / sum(abs(Bold[[k]]))
    }
    if(trace) print(data.frame(step, rho, eps))
    if(eps <= thr) break
    if(update) {if(step %% 10 == 0) rho <- rho * rho.increment}
  }
  df <- vector(mode = "list", length = K)
  for(k in seq_len(K)) {
    df[[k]] <- apply(abs(B[[k]][-1, , drop = FALSE]) != 0, 2, sum)
    df[[k]] <- c(df[[k]], sum(df[[k]]))
  }
  
  return(list(iter = step, B = B, df = df, nit = step, conv = 1*(step == maxit), diff = eps))
}

list2array <- function(x) {
  if(!is.list(x)) stop("x is not a list!")
  K <- length(x)
  if(is.matrix(x[[1]])) {
    n <- sapply(x, nrow)
    p <- ncol(x[[1]])
    y <- array(0, dim = c(max(n), p, K), dimnames = list(NULL, colnames(x[[1]]), names(x)))
    for(k in seq_len(K)){ y[seq_len(n[k]), , k] <- x[[k]] }
  } else {
    n <- length(x[[1]])
    y <- matrix(unlist(x), nrow = n, ncol = K, dimnames = list(names(x[[1]]), names(x)))
  }
  y
}
array2list <- function(x) {
  if(!is.array(x)) stop("x is not a list!")
  K <- tail(dim(x), 1)
  y <- vector(mode = "list", length = K)
  for(k in seq_len(K)) y[[k]] <- if(length(dim(x)) == 2) x[, k] else na.omit(x[,,k])
  names(y) <- tail(dimnames(x), 1)[[1]]
  y
}

##### jcglasso S3 methods #####
print.jcglasso <- function (x, digits = 3L, ...){
  K <- length(x$Z)
  p <- x$nresp
  q <- x$npred
  nrho <- x$nrho
  rho1 <- x$rho
  rho2 <- x$rho2
  nlambda <- x$nlambda
  lambda1 <- x$lambda
  lambda2 <- x$lambda2
  alpha <- x$alpha
  rho.x <- x$nu
  
  dots <- list(...)
  if (is.null(dots$print.gap)) dots$print.gap <- 2L
  if (is.null(dots$quote)) dots$quote <- FALSE
  if (is.null(dots$row.names)) dots$row.names <- FALSE
  
  if (nrho == 1L | nlambda == 1L) {
    df.B <- lapply(x$dfB, function(x) x[p + 1L, , ] +  p)
    df.Tht <- lapply(x$dfTht, function(x) drop(x) + p + q)
    df <- sapply(seq_len(K), function(i) df.B[[i]] + df.Tht[[i]])
    if(is.vector(df)) df <- t(df)
    df.max <- p * (p + 1) / 2 + (q + 1L) * p + q * (q + 1) / 2
    df.per <- formatC(round(df / df.max * 100, digits = digits), format = 'f', digits = digits)
    df.per <- apply(df.per, 2, function(.df.per) paste("(", .df.per, "%)", sep = ""))
    if(is.vector(df.per)) df.per <- t(df.per)
    # ncomp <- sapply(x$InfoStructure$ncomp, function(x) as.vector(t(x)))
    # colnames(ncomp) <- paste0("N. Comp. Group", seq_len(K))
    colnames(df.per) <- paste0("(df%) Group", seq_len(K))
    colnames(df) <- paste0("df Group", seq_len(K))
    tbl <- data.frame("nu" = if(is.null(rho.x)) rho1+rho2 else rho.x, 
                      "rho" = rho1+rho2, "lambda" = lambda1+lambda2, "alpha" = alpha,
                      df, df.per, fix.empty.names = FALSE, check.names = FALSE)
    # names(tbl) <- c("lambda", "rho1", "rho2", "df", "(df%)", "N. Comp.")
    if (q == 0L) tbl <- tbl[, -c(1L,3L), drop = FALSE]
    cat("\nCall:  ", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
    do.call(function(...) print.data.frame(tbl, digits = digits, ...), dots)
    cat("\n---\n")
  } else {
    lambda1 <- rep(lambda1, each = nrho)
    lambda2 <- rep(lambda2, each = nrho)
    rho1 <- rep(rho1, times = nlambda)
    rho2 <- rep(rho2, times = nlambda)
    df.B <- lapply(x$dfB, function(x) as.vector(t(x[p + 1L, , ])) +  p)
    df.Tht <- lapply(x$dfTht, function(x) as.vector(t(x)) + p + q)
    df <- sapply(seq_len(K), function(i) df.B[[i]] + df.Tht[[i]])
    df.max <- p * (p + 1) / 2 + (q + 1L) * p + q * (q + 1) / 2
    df.per <- formatC(round(df / df.max * 100, digits = digits), format = 'f', digits = digits)
    df.per <- apply(df.per, 2, function(x) paste("(", x, "%)", sep = ""))
    # ncomp <- sapply(x$InfoStructure$ncomp, function(x) as.vector(t(x)))
    # colnames(ncomp) <- paste0("N. Comp. Group", seq_len(K))
    colnames(df.per) <- paste0("(df%) Group", seq_len(K))
    colnames(df) <- paste0("df Group", seq_len(K))
    tbl <- data.frame("nu" = if(is.null(rho.x)) rho1+rho2 else rho.x, 
                      "rho" = rho1+rho2, "lambda" = lambda1+lambda2, "alpha" = alpha,
                      df, df.per, fix.empty.names = FALSE, check.names = FALSE)
    # names(tbl) <- c("rho1", "rho2", "df", "", "N. Comp.")
    if (q == 0L) tbl <- tbl[, -c(1L,3L), drop = FALSE]
    cat("\nCall:  ", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
    #tbl.list <- split(tbl, f = rep(seq_len(nlambda), each = nrho), drop = FALSE)
    if (nlambda <= nrho) f <- rep(seq_len(nlambda), each = nrho)
    else f <- rep(seq_len(npen1), times = npen2)
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

coef.jcglasso <- function(object, type = c("all", "B", "Theta", "Omega"), class.id, rho.id, lambda.id, drop = TRUE, ...) {
  type <- match.arg(type)
  K <- length(object$nobs)
  nrho <- object$nrho
  nlambda <- object$nlambda
  id_X <- object$InfoStructure$id_X
  id_Y <- object$InfoStructure$id_Y
  if (!is.logical(drop)) stop(sQuote("drop"), " is not an object of type ", sQuote("integer"))
  if (missing(class.id)) class.id <- seq_len(K)
  else {
    if(!is.vector(class.id)) stop(sQuote("class.id"), " is not a vector")
    if(any(abs(as.integer(class.id) - class.id) > 0)) stop(sQuote("class.id"), " is not an object of type ", dQuote("integer"))
    if(min(class.id) <= 0) stop("some entry in ", sQuote("class.id"), " is not a positive integer")
    if(max(class.id) > K) stop("some entry in ", sQuote("class.id"), " is larger than ", sQuote(class.id))
  }
  
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

  if (type == "all") {
    out.coef <- list(B = object$B[class.id],
                     # Sigma = object$Sgm[class.id],
                     Omega = object$Omega[class.id],
                     Theta = lapply(object$Tht, function(x) x[id_Y, id_Y, , , drop = FALSE])[class.id])
    out.coef$B <- lapply(out.coef$B, function(x) x[, , lambda.id, rho.id])
    # out.coef$Sigma <- lapply(out.coef$Sigma, function(x) x[, , lambda.id, rho.id])
    out.coef$Omega <- lapply(out.coef$Omega, function(x) x[, , lambda.id, rho.id])
    out.coef$Theta <- lapply(out.coef$Theta, function(x) x[, , lambda.id, rho.id])
    if(drop & length(class.id) == 1L) out.coef <- lapply(out.coef, function(x) x[[1]])
  } else {
    out.coef <- switch(type,
                       B = object$B[class.id],
                       # Sigma = object$Sgm[class.id],
                       Omega = object$Omega[class.id],
                       Theta = lapply(object$Tht, function(x) x[id_Y, id_Y, , , drop = FALSE])[class.id])
    out.coef <- lapply(out.coef, function(x) x[, , lambda.id, rho.id])
    if(drop & length(class.id) == 1L) out.coef <- out.coef[[1]]
  }
  out.coef
}

fitted.jcglasso <- function(object, class.id, rho.id, lambda.id, drop = TRUE, ...) {
  K <- length(object$nobs)
  nrho <- object$nrho
  nlambda <- object$nlambda
  id_X <- object$InfoStructure$id_X
  id_Y <- object$InfoStructure$id_Y
  if (!is.logical(drop)) stop(sQuote("drop"), " is not an object of type ", sQuote("integer"))
  if (missing(class.id)) class.id <- seq_len(K)
  else {
    if(!is.vector(class.id)) stop(sQuote("class.id"), " is not a vector")
    if(any(abs(as.integer(class.id) - class.id) > 0)) stop(sQuote("class.id"), " is not an object of type ", dQuote("integer"))
    if(min(class.id) <= 0) stop("some entry in ", sQuote("class.id"), " is not a positive integer")
    if(max(class.id) > K) stop("some entry in ", sQuote("class.id"), " is larger than ", sQuote(class.id))
  }
  
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
  out.fitted <- object$mu[class.id]
  out.fitted <- lapply(out.fitted, function(x) x[, id_Y, lambda.id, rho.id])
  if(drop & length(class.id) == 1L) out.fitted <- out.fitted[[1]]
  out.fitted
}

residuals.jcglasso <- function(object, type = c("observed", "working"), class.id, rho.id, lambda.id, drop = TRUE, ...) {
  type <- match.arg(type)
  K <- length(object$nobs)
  nrho <- object$nrho
  nlambda <- object$nlambda
  
  if (!is.logical(drop)) stop(sQuote("drop"), " is not an object of type ", sQuote("integer"))
  if (missing(class.id)) class.id <- seq_len(K)
  else {
    if(!is.vector(class.id)) stop(sQuote("class.id"), " is not a vector")
    if(any(abs(as.integer(class.id) - class.id) > 0)) stop(sQuote("class.id"), " is not an object of type ", dQuote("integer"))
    if(min(class.id) <= 0) stop("some entry in ", sQuote("class.id"), " is not a positive integer")
    if(max(class.id) > K) stop("some entry in ", sQuote("class.id"), " is larger than ", sQuote(class.id))
  }
  
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
  
  out.residuals <- object$R
  if (type == "observed") {
    for(k in class.id){
      R <- event2(object$Z[[k]])
      residuals.dim <- dim(out.residuals[[k]])
      for (i in 1:residuals.dim[3L]) {
        for (j in 1:residuals.dim[4L]) out.residuals[[k]][, , i, j][R != 0] <- NA
      }
    }
  }
  out.residuals <- out.residuals[class.id]
  out.residuals <- lapply(out.residuals, function(x) x[, , lambda.id, rho.id])
  if(drop & length(class.id) == 1L) out.residuals <- out.residuals[[1]]
  out.residuals
}

summary.jcglasso <- function(object, GoF = AIC, print.all = TRUE, digits = 3L,  ...){
  if (!is.element(class(GoF), c("function", "GoF2")))
    stop (sQuote("GoF"), " is not a goodness-of-fit function (AIC or BIC) neither an object of class ", sQuote("GoF"))
  dots <- list(...)
  n <- object$nobs
  p <- object$nresp
  q <- object$npred
  K <- length(n)
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
  rho1 <- object$rho
  rho2 <- object$rho2
  nrho <- object$nrho
  lambda1 <- object$lambda
  lambda2 <- object$lambda2
  nlambda <- object$nlambda
  rho.x <- object$nu
  alpha <- object$alpha
  if (is.null(dots$print.gap)) dots$print.gap <- 2L
  if (is.null(dots$quote)) dots$quote <- FALSE
  if (is.null(dots$row.names)) dots$row.names <- FALSE
  
  if (nrho == 1L | nlambda == 1L) {
    df.B <- sapply(object$dfB, function(x) as.vector(t(x[p + 1L, , ])) + p)
    if(is.vector(df.B)) df.B <- t(df.B)
    df.B <- rowSums(df.B)
    df.Tht <- sapply(object$dfTht, function(x) as.vector(t(x)) + p + q)
    if(is.vector(df.Tht)) df.Tht <- t(df.Tht)
    df.Tht <- rowSums(df.Tht)
    df <- df.B + df.Tht
    df.max <- p * (p + 1) / 2 + (q + 1L) * p + q * (q + 1) / 2
    df.max <- df.max * K
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
    tbl <- data.frame(nu = if(is.null(rho.x)) rho1+rho2 else rho.x, rho = rho1+rho2, 
                      lambda = lambda1+lambda2, alpha = alpha, df.B, df.Tht, df, df.per, val, Rank = rnk)
    names(tbl)[8L] <- "(df%)"
    names(tbl)[9L] <- GoF$type
    if (q == 0L) tbl <- tbl[, -c(1L,3L), drop = FALSE]
    cat("\nCall:  ", paste(deparse(object$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
    if (print.all) do.call(function(...) print.data.frame(tbl, digits = digits, ...), dots)
  } 
  else {
    lambda1 <- rep(lambda1, each = nrho)
    lambda2 <- rep(lambda2, each = nrho)
    lambda.id <- rep(seq_len(nlambda), each = nrho)
    rho1 <- rep(rho1, times = nlambda)
    rho2 <- rep(rho2, times = nlambda)
    rho.id <- rep(seq_len(nrho), times = nlambda)
    df.B <- sapply(object$dfB, function(x) as.vector(t(x[p + 1L, , ])) + p)
    df.B <- rowSums(df.B)
    df.Tht <- sapply(object$dfTht, function(x) as.vector(t(x)) + p+ q)
    df.Tht <- rowSums(df.Tht)
    df <- df.B + df.Tht
    df.max <- p * (p + 1) / 2 + (q + 1L) * p + q * (q + 1) / 2
    df.max <- df.max * K
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
    tbl <- data.frame(nu = if(is.null(rho.x)) rho1+rho2 else rho.x, rho = rho1+rho2, 
                      lambda = lambda1+lambda2, alpha = alpha, df.B, df.Tht, df, df.per, val, Rank = rnk)
    names(tbl)[8L] <- "(df%)"
    names(tbl)[9L] <- GoF$type
    if (q == 0L) tbl <- tbl[, -c(1L,3L), drop = FALSE]
    cat("\nCall:  ", paste(deparse(object$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
    tbl.list <- split(tbl, f = rep(seq_len(nlambda), each = nrho))
    if (print.all) do.call(function(...) print.listof(tbl.list, digits = digits, ...), dots)
  }
  
  lbl <- c("model", "nClasses", "nObs", "nResp", "nPred", "penalty", "lambda", "lambda.id",
           "rho", "rho.id", "nu", "alpha", GoF$type, "df.B", "df.Tht", "df")  
  
  lbl <- paste("\n", format(lbl, justify = "right"), ":", sep = "")
  cat("\n===============================================================")
  cat("\n\nSummary of the Selected Model\n")
  if (!is.null(object$model)) cat(lbl[1L], sQuote(object$model))
  cat(lbl[2L], K)
  cat(lbl[3L], n)
  cat(lbl[4L], p)
  cat(lbl[5L], q)
  cat(lbl[6L], object$penalty)
  if (q > 0L) {
    cat(lbl[7L], tbl[rnk_min, "lambda"]) 
    cat(lbl[8L], lambda.id)
  }
  cat(lbl[9], tbl[rnk_min, "rho"])
  cat(lbl[10L], rho.id)
  if (q > 0L) cat(lbl[11L], tbl[rnk_min, "nu"])
  cat(lbl[12L], tbl[rnk_min, "alpha"])
  cat(lbl[13L], tbl[rnk_min, GoF$type])
  cat(lbl[14L], tbl[rnk_min, "df.B"])
  cat(lbl[15L], tbl[rnk_min, "df.Tht"])
  cat(lbl[16L], tbl[rnk_min, "df"])
  cat("\n\n===============================================================\n\n")
  invisible(list(table = tbl, rho.id = rho.id, lambda.id = lambda.id))
}

plot.jcglasso <- function(x, what = c("Theta", "diag(Theta)", "b0", "B"), penalty = ifelse(nrho >= nlambda, "rho", "lambda"),
                          given = NULL, GoF = AIC, add.labels, matplot.arg1, matplot.arg2, labels.arg, abline.arg, mtext.arg, save.plot,
                          grdev = pdf, grdev.arg, digits = 4L, ...) {
  p <- x$nresp
  q <- x$npred
  K <- length(x$nobs)
  nrho <- x$nrho
  nlambda <- x$nlambda
  # if (nrho1 == 1L & nrho2 == 1L) stop("plot method is not available because ", sQuote("nrho2 = 1"), " and ", sQuote("nrho1 = 1"))
  if (nrho == 1L & nlambda == 1L) { plot.jcglasso2(x, ...); return(invisible(NULL)) }
  if (inherits(what, "formula")) {
    xy <- function2xy.v2(what)
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
  if (is.null(names(matplot.arg1))) stop(sQuote("matplot.arg1"), " is not a named list")
  # testing 'matplot.arg2'
  if (missing(matplot.arg2)) matplot.arg2 <- vector(mode = "list")
  if (is.null(matplot.arg2$type)) matplot.arg2$type <- "l"
  if (is.null(matplot.arg2$col)) matplot.arg2$col <- "gray70"
  if (is.null(matplot.arg2$lty)) matplot.arg2$lty <- 2L
  if (!is.list(matplot.arg2)) stop(sQuote("matplot.arg2"), " is not an object of type ", dQuote("list"))
  if (is.null(names(matplot.arg2))) stop(sQuote("matplot.arg2"), " is not a named list")
  # testing 'labels.arg'
  if (missing(labels.arg)) labels.arg <- vector(mode = "list")
  if (is.null(labels.arg$pos)) labels.arg$pos <- 4L
  if (!is.list(labels.arg)) stop(sQuote("labels.arg"), " is not an object of type ", dQuote("list"))
  if (is.null(names(labels.arg))) stop(sQuote("labels.arg"), " is not a named list")
  # testing 'abline.arg'
  if (missing(abline.arg)) abline.arg <- vector(mode = "list")
  if (is.null(abline.arg$lwd)) abline.arg$lwd <- 2L
  if (is.null(abline.arg$lty)) abline.arg$lty <- 2L
  if (is.null(abline.arg$col)) abline.arg$col <- 2L
  if (!is.list(abline.arg)) stop(sQuote("abline.arg"), " is not an object of type ", dQuote("list"))
  if (is.null(names(abline.arg))) stop(sQuote("abline.arg"), " is not a named list")
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
  } 
  else {
    if (ngiven > 1L | what == "B") {
      op <- par(no.readonly = TRUE)
      on.exit(par(op), add = TRUE)
      devAskNewPage(TRUE)
    }
  }
  if (what == "Theta") U <- .row(c(p, p)) < .col(c(p, p))
  if (add.labels) {
    if (what == "b0") lbls <- colNames(x$Z[[1]])$Y
    if (what == "B") lbls <- colNames(x$Z[[1]])$X
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
      Y <- coef.jcglasso(x, type = coef.name, rho.id = given[m], drop = TRUE)
      if (!is.null(GoF))  minGoF.id <- which.min(GoF$value_gof[, given[m]])
    }
    else {
      if (q != 0) 
        tp.lab <- bquote(rho ~ " | {" ~ lambda %~~% .(tp.given[m]) ~ "}")
      else tp.lab <- bquote(rho)
      Y <- coef.jcglasso(x, type = coef.name, lambda.id = given[m], drop = TRUE)
      if (!is.null(GoF)) minGoF.id <- which.min(GoF$value_gof[given[m], ])
    }
    ##################
    # plotting 'b0'
    if (what == "b0") {
      hh <- hh + 1L
      if (save.plot) {
        if (is.null(grdev.arg)) grdev(file = file.names[hh])
        else do.call(function(...) grdev(file = file.names[hh], ...), grdev.arg)
      } else dev.hold()
      
      xlab <- if (is.null(matplot.arg1$xlab)) tp.lab else matplot.arg1$xlab
      ylab <- if (is.null(matplot.arg1$ylab))  "" else matplot.arg1$ylab
      main <- if (is.null(matplot.arg1$main)) ifelse(q == 0L, "Expected Values", "Intercepts") else matplot.arg1$main
      arg.name <- setdiff(names(matplot.arg1), c("x", "y", "xlab", "ylab", "main"))
      for(k in seq_len(K)){
        if (q == 0L) yy <- t(Y[[k]])
        else yy <- t(Y[[k]][1L, , ])
        do.call(what = function(...) matplot(x = tp, y = yy, xlab = xlab, ylab = ylab, main = main, ...), args = matplot.arg1[arg.name])
        if (!is.null(GoF)) {
          if(add.labels) {
            labels <- if(is.null(labels.arg$labels)) lbls else labels.arg$labels
            arg.name2 <- setdiff(names(labels.arg), "labels")
            do.call(what = function(...) text(x = tp[minGoF.id], y = yy[minGoF.id, ], labels = labels, ...), args = labels.arg[arg.name2])
          }
          abline.arg$v <- tp[minGoF.id]
          do.call(what = abline, args = abline.arg)
          mtext.arg$text <- GoF$type
          mtext.arg$at <- tp[minGoF.id]
          do.call(what = mtext, args = mtext.arg)
        }
      }
      
      if (save.plot) setTxtProgressBar(pb, hh) else dev.flush()
      if(save.plot) dev.off()
    }
    ##################
    # plotting 'B'
    if (what == "B") {
      varname <- colNames(x$Z[[1L]])$Y
      for (k in seq_len(p)) {
        hh <- hh + 1L
        if (save.plot) {
          if (is.null(grdev.arg)) grdev(file = file.names[hh])
          else do.call(function(...) grdev(file = file.names[hh], ...), grdev.arg)
        } else dev.hold()
        
        xlab <- if (is.null(matplot.arg1$xlab)) tp.lab else matplot.arg1$xlab
        ylab <- if (is.null(matplot.arg1$ylab)) "Regression Coefficients" else matplot.arg1$ylab
        main <- if (is.null(matplot.arg1$main)) paste0("Response Variable ", varname[k]) else matplot.arg1$main
        arg.name <- setdiff(names(matplot.arg1), c("x", "y", "xlab", "ylab", "main"))
        for(o in seq_len(K)) {
          yy <- Y[[o]][-1L, k, ]
          if (!is.vector(yy)) yy <- t(yy)
          if (is.null(GoF)) do.call(what = function(...) matplot(x = tp, y = yy, xlab = xlab, ylab = ylab, main = main, ...), args = matplot.arg1[arg.name])
          else {
            matplot(x = tp, y = yy, xlab = xlab, ylab = ylab, main = main, type = "n")
            A <- abs(yy[minGoF.id, ]) > 0
            if(any(!A)) do.call(what = function(...) matpoints(x = tp, y = yy[, !A], ...), args = matplot.arg2[arg.name])
            if(any(A)) {
              do.call(what = function(...) matpoints(x = tp, y = yy[, A], ...), args = matplot.arg1[arg.name])
              if(add.labels) {
                labels <- if (is.null(labels.arg$labels)) lbls[A] else labels.arg$labels
                arg.name2 <- setdiff(names(labels.arg), "labels")
                do.call(what = function(...) text(x = tp[minGoF.id], y = yy[minGoF.id, A], labels = labels, ...), args = labels.arg[arg.name2])
              }
            }
            abline.arg$v <- tp[minGoF.id]
            do.call(what = abline, args = abline.arg)
            mtext.arg$text <- GoF$type
            mtext.arg$at <- tp[minGoF.id]
            do.call(what = mtext, args = mtext.arg)
          }
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
      
      xlab <- if (is.null(matplot.arg1$xlab)) tp.lab else matplot.arg1$xlab
      ylab <- if (is.null(matplot.arg1$ylab)) "Diagonal Elements" else matplot.arg1$ylab
      main <- if (is.null(matplot.arg1$main)) "Precision Matrix" else matplot.arg1$main
      arg.name <- setdiff(names(matplot.arg1), c("x", "y", "xlab", "ylab", "main"))
      for(k in seq_len(K)){
        yy <- t(apply(Y[[k]], 3L, diag))
        do.call(what = function(...) matplot(x = tp, y = yy, xlab = xlab, ylab = ylab, main = main, ...), args = matplot.arg1[arg.name])
        if (!is.null(GoF)) {
          if(add.labels) {
            labels <- if(is.null(labels.arg$labels)) lbls else labels.arg$labels
            arg.name2 <- setdiff(names(labels.arg), "labels")
            do.call(what = function(...) text(x = tp[minGoF.id], y = yy[minGoF.id, ], labels = labels, ...), args = labels.arg[arg.name2])
          }
          abline.arg$v <- tp[minGoF.id]
          do.call(what = abline, args = abline.arg)
          mtext.arg$text <- GoF$type
          mtext.arg$at <- tp[minGoF.id]
          do.call(what = mtext, args = mtext.arg)
        }
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
      
      xlab <- if (is.null(matplot.arg1$xlab)) tp.lab else matplot.arg1$xlab
      ylab <- if (is.null(matplot.arg1$ylab)) "Off-Diagonal Elements" else matplot.arg1$ylab
      main <- if (is.null(matplot.arg1$main)) "Precision Matrix" else matplot.arg1$main
      arg.name <- setdiff(names(matplot.arg1), c("x", "y", "xlab", "ylab", "main"))
      for(k in seq_len(K)) {
        yy <- t(apply(Y[[k]], 3L, function(M) M[U]))
        if (is.null(GoF)) do.call(what = function(...) matplot(x = tp, y = yy, xlab = xlab, ylab = ylab, main = main, ...), args = matplot.arg1[arg.name])
        else {
          matplot(x = tp, y = yy, xlab = xlab, ylab = ylab, main = main, type = "n")
          A <- abs(yy[minGoF.id, ]) > 0
          if(any(!A)) do.call(what = function(...) matpoints(x = tp, y = yy[, !A], ...), args = matplot.arg2[arg.name])
          if(any(A)) {
            do.call(what = function(...) matpoints(x = tp, y = yy[, A], ...), args = matplot.arg1[arg.name])
            if(add.labels) {
              labels <- if(is.null(labels.arg$labels)) lbls[A] else labels.arg$labels
              arg.name2 <- setdiff(names(labels.arg), "labels")
              do.call(what = function(...) text(x = tp[minGoF.id], y = yy[minGoF.id, A], labels = labels, ...), args = labels.arg[arg.name2])
            }
          }
          abline.arg$v <- tp[minGoF.id]
          do.call(what = abline, args = abline.arg)
          mtext.arg$text <- GoF$type
          mtext.arg$at <- tp[minGoF.id]
          do.call(what = mtext, args = mtext.arg)
        }
      }
      
      if (save.plot) setTxtProgressBar(pb, hh) else dev.flush()
      if(save.plot) dev.off()
    }
  }
  invisible(NULL)
}

function2xy.v2 <- function (frml) {
  x <- match.arg(all.vars(frml)[2L], c("rho", "lambda"))
  y <- as.character(terms(frml)[[2L]])
  if (length(y) == 2L) {
    if (y[1L] != "diag" | y[2L] != "Theta") 
      stop(sQuote(paste0(y[1L], "(", y[2L], ")")), " is not available as response. Please, use ", 
           sQuote("diag(Theta)"))
    y <- "diag(Theta)"
  }
  given <- as.character(terms(frml)[[3L]])
  if (length(given) == 1L) 
    given <- NULL
  else {
    if (given[3L] == "rho") {
      if (x == "rho") 
        stop("You can not condition on ", sQuote(given[3L]), 
             ". Please, use ", sQuote("lambda"))
      return(list(x = x, y = y, given = NULL))
    }
    if (given[3L] == "lambda") {
      if (x == "lambda") 
        stop("You can not condition on ", sQuote(given[3L]), 
             ". Please, use ", sQuote("rho"))
      return(list(x = x, y = y, given = NULL))
    }
    given.terms <- gsub("[[:blank:]]", "", given[3L])
    given.terms <- unlist(strsplit(given.terms, split = "\n", 
                                   fixed = TRUE))
    if (length(given.terms) == 1L) 
      given <- eval(parse(text = given.terms))
    else {
      given.terms <- given.terms[2L]
      given.name <- unlist(strsplit(given.terms, "=", fixed = TRUE))[1L]
      if (!is.element(given.name, c("lambda.id", "rho.id"))) 
        stop("Please, use ", sQuote(ifelse(x == "rho", "lambda.id", "rho.id")), " instead of ", sQuote(given.name))
      if ((x == "rho" & given.name == "rho.id") | (x == "lambda" & given.name == "lambda.id")) 
        stop("You can not condition on ", sQuote(given.name), 
             ". Please, use ", sQuote(ifelse(given.name == "rho.id", "lambda.id", "rho.id")))
      given <- eval(parse(text = given.terms))
    }
  }
  out <- list(penalty = x, what = y, given = given)
  out
}

plot.jcglasso2 <- function(x, type = c("Gyy", "Gxy", "Gxx"), ...) {
  K <- length(x$Z)
  type <- match.arg(type)
  id_X <- x$InfoStructure$id_X
  id_Y <- x$InfoStructure$id_Y
  if (is.null(x$npred > 0L) & is.element(type, c("Gxy", "Gxx")))
    stop(sQuote(type), " is not available. Please, use type = ", dQuote("Gyy"))
  mfrow <- if(K == 1) c(1, 2) else c(2, K)
  par(mfrow = mfrow, pty = "s", omi = c(0.3, 0.3, 0.3, 0.3), mai = c(0.3, 0.3, 0.3, 0.3))
  if(type == "Gyy") {
    lapply(x$InfoStructure$Adj, function(x) {diag(x[id_Y, id_Y, 1L, 1L]) <- 1; image(x[id_Y, id_Y, 1L, 1L], col = gray.colors(256), main = "Adjacency Matrix")})
    g <- lapply(x$InfoStructure$Adj, function(x) graph.adjacency(x[id_Y, id_Y, 1L, 1L], mode = "undirected", diag = FALSE))
    layout.grid = lapply(g, layout.fruchterman.reingold)
    sapply(seq_len(length(g)), function(i) plot(g[[i]], layout = layout.grid[[i]], edge.color = "gray50", 
                                                vertex.color = "red", vertex.size = 3, vertex.label = NA, 
                                                main = "Graph Pattern"))
  } 
  if(type == "Gxx") {
    lapply(x$InfoStructure$Adj, function(x) {diag(x[id_X, id_X, 1L, 1L]) <- 1; image(x[id_X, id_X, 1L, 1L], col = gray.colors(256), main = "Adjacency Matrix")})
    g <- lapply(x$InfoStructure$Adj, function(x) graph.adjacency(x[id_X, id_X, 1L, 1L], mode = "undirected", diag = FALSE))
    layout.grid = lapply(g, layout.fruchterman.reingold)
    sapply(seq_len(length(g)), function(i) plot(g[[i]], layout = layout.grid[[i]], edge.color = "gray50", 
                                                vertex.color = "red", vertex.size = 3, vertex.label = NA, 
                                                main = "Graph Pattern"))
  }
  if(type == "Gxy") {
    lapply(x$InfoStructure$Adj, function(x) image(x[id_X, id_Y, 1L, 1L], col = gray.colors(256), main = "Adjacency Matrix"))
    g <- lapply(x$InfoStructure$Adj, function(x) {
      out <- x[id_X, id_Y, 1L, 1L]
      nmsX <- rownames(out)
      nmsY <- colnames(out)
      eB <- data.frame(from = rep(nmsX, ncol(out)), to = rep(nmsY, each = nrow(out)))[c(out) != 0, ]
      vB <- data.frame(id = c(nmsY, nmsX))
      graph_from_data_frame(d = eB, directed = TRUE, vertices = vB)
    })
    layout.grid = lapply(g, layout.fruchterman.reingold)
    sapply(seq_len(length(g)), function(i) plot(g[[i]], layout = layout.grid[[i]], edge.color = "gray50", 
                                                vertex.color = "red", vertex.size = 3, vertex.label = NA, 
                                                main = "Graph Pattern"))
  }
  
  invisible(NULL)
}

##### jcggm function #####
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
    wTht[[i]][id_Y, id_Y][object$InfoStructure$Adj[[i]][id_Y, id_Y, 1L, 1L] == 0] <- .Machine$double.xmax
    if(q > 0L) {
      wTht[[i]][id_X, id_X][object$InfoStructure$Adj[[i]][id_X, id_X, 1L, 1L] == 0] <- .Machine$double.xmax
      wB[[i]][object$InfoStructure$Adj[[i]][id_X, id_Y, 1L, 1L] == 0] <- .Machine$double.xmax
    }
  }
  
  rho1.true <- object$rho
  rho2.true <- object$rho2
  lambda1.true <- object$lambda
  lambda2.true <- object$lambda2
  alpha <- object$alpha
  model <- paste0("MLE of a ", object$model)
  
  object$weights.B <- wB
  object$weights.Tht <- wTht
  
  if(!refit) {
    out <- jcggm.fit(object, tp.min = tp.min, trace = trace, ...)
  }
  else {
    out <- jcglasso.fit(Z = object$Z, diagonal = object$diagonal, weights.B = wB, weights.Tht = wTht, 
                        penalty = object$penalty, nrho = nrho, rho.min.ratio = object$rho.min.ratio, 
                        rho = tp.min, nlambda = nlambda, lambda.min.ratio = object$lambda.min.ratio, 
                        lambda = tp.min, nu = tp.min, alpha = object$alpha, maxit.em = object$maxit.em, 
                        thr.em = object$thr.em, maxit.bcd = object$maxit.bcd, thr.bcd = object$thr.bcd, 
                        trace = trace, offset = object$offset, covar2corr = object$covar2corr, 
                        truncate = object$truncate)
  }
 
  InfoStructure <- list(Adj = out$Adj, ncomp = out$ncomp, Ck = out$Ck, pk = out$pk,
                        id_X = out$id_X, id_Y = out$id_Y)
  for(i in seq_len(num_class)) {
    wTht[[i]][wTht[[i]] == .Machine$double.xmax] <- +Inf
    if(q > 0L) wB[[i]][wB[[i]] == .Machine$double.xmax] <- +Inf
  }
  out.jcggm <- list(call = this.call, Zipt = out$Zipt, B = out$B,
                    mu = out$mu, R = out$R, S = out$S, Sgm = out$Sgm, Tht = out$Tht,
                    Omega = out$Omega, dfB = out$dfB, dfTht = out$dfTht, InfoStructure = InfoStructure,
                    nit = out$nit, Z = object$Z, diagonal = object$diagonal, 
                    weights.B = wB, weights.Tht = wTht,
                    nlambda = nlambda, lambda.min.ratio = object$lambda.min.ratio, lambda = lambda1.true, 
                    nrho = nrho, rho.min.ratio = object$rho.min.ratio, rho = rho1.true, 
                    nu = object$nu, alpha = object$alpha, lambda2 = lambda2.true, rho2 = rho2.true, 
                    connected = out$connected, penalty = object$penalty,
                    model = model, maxit.em = object$maxit.em, thr.em = object$thr.em,
                    maxit.bcd = object$maxit.bcd, thr.bcd = object$thr.bcd, conv = out$conv,
                    offset = object$offset, covar2corr = object$covar2corr, truncate = object$truncate,
                    subrout = out$subrout, trace = trace, nobs = object$nobs, nresp = object$nresp, npred = object$npred)
  
  # out.cjggm <- list(call = this.call, Yipt = out$Yipt, Xipt = out$Xipt, B = out$B,
  #                   mu = out$mu, R = out$R, S = out$S, Sgm = out$Sgm, Tht = out$Tht,
  #                   Sxx = out$Sxx, Sxy = out$Sxy, Sgmxx = out$Sgmxx, Sgmxy = out$Sgmxy, 
  #                   Thtxx = out$Thtxx, Thtxy = out$Thtxy,
  #                   dfB = out$dfB, dfTht = out$dfTht, InfoStructure = InfoStructure,
  #                   nit = out$nit, Z = object$Z, diagonal = out$diagonal, weights.B = wB,
  #                   weights.Tht = wTht, nlambda = nlambda, lambda.min.ratio = object$lambda.min.ratio, 
  #                   lambda = lambda1.true, nrho = nrho, rho.min.ratio = object$rho.min.ratio, 
  #                   rho = rho1.true, rho.x = object$rho.x, alpha = object$alpha,
  #                   lambda2 = lambda2.true, rho2 = rho2.true, connected = out$connected, 
  #                   penalty = object$penalty, model = model, maxit.em = object$maxit.em, thr.em = object$thr.em,
  #                   maxit.bcd = object$maxit.bcd, thr.bcd = object$thr.bcd, conv = out$conv,
  #                   offset = object$offset, covar2corr = object$covar2corr, truncate = object$truncate,
  #                   subrout = out$subrout, trace = trace, 
  #                   nobs = object$nobs, nresp = object$nresp, npred = object$npred)
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
  
  weights.Tht <- object$weights.Tht
  weights.B <- object$weights.B

  alpha <- object$alpha
  storage.mode(tp.min) <- "double"
  storage.mode(trace) <- "integer"
  
  nrho <- 1L
  nlambda <- 1L
  
  K <- length(Z)
  storage.mode(K) <- "integer"
  
  n <- sapply(Z, nobs)
  storage.mode(n) <- "integer"
  weights <- n / sum(n)
  storage.mode(weights) <- "double"
  
  p <- object$nresp
  storage.mode(p) <- "integer"
  
  q <- object$npred
  storage.mode(q) <- "integer"
  
  k_lab <- paste0("class", seq_len(K))
  rho_lab <- paste0("rho_", seq_len(nrho))
  lambda_lab <- paste0("lambda_", seq_len(nlambda))
  
  ynames <- dimnames(Z[[1]])$Y[[2]]
  # yrnames <- dimnames(Z[[1]])$Y[[1]]
  xnames <- dimnames(Z[[1]])$X[[2]]
  znames <- c(ynames, xnames)
  
  X.null <- q == 0L
  
  seq_k <- seq_len(K)
  id_Y <- object$InfoStructure$id_Y
  id_X <- object$InfoStructure$id_X
  dim_Z <- p + q
  id_Z <- seq_len(dim_Z)
  storage.mode(dim_Z) <- "integer"
  
  ##### computing output and working objects #####
  Zmat <- vector(mode = "list", length = K)
  zm <- vector(mode = "list", length = K)
  zv <- vector(mode = "list", length = K)
  loz <- vector(mode = "list", length = K)
  upz <- vector(mode = "list", length = K)
  Id <- vector(mode = "list", length = K)
  InfoP <- vector(mode = "list", length = K)
  nP <- vector(mode = "list", length = K)
  Zipt <- object$Zipt
  B <- object$B
  mu <- object$mu
  R <- object$R
  S <- object$S
  Sgm <- object$Sgm
  Tht <- object$Tht
  Thtxx <- object$Omega
  B_n <- mu_n <- R_n <- S_n <- Sgm_n <- Tht_n <- vector(mode = "list", length = K)
  Adj <- object$InfoStructure$Adj
  dfB <- object$dfB
  dfTht <- object$dfTht
  ncomp <- object$InfoStructure$ncomp
  Ck <- object$InfoStructure$Ck
  pk <- object$InfoStructure$Ck
  Zipt_lo <- Zipt_up <- Zipt_n <- vector(mode = "list", length = K)
  T1o <- T2o <- vector(mode = "list", length = K)
  T1 <- double(dim_Z)
  T2 <- vector(mode = "list", length = K)
  rho1mat <- vector(mode = "list", length = K)
  for(k in seq_k){
    Zmat[[k]] <- cbind(getMatrix2(Z[[k]], name = "Y", ordered = TRUE), 
                       getMatrix2(Z[[k]], name = "X", ordered = TRUE))
    Zmat[[k]][is.na(Zmat[[k]])] <- 0
    storage.mode(Zmat[[k]]) <- "double"
    
    zm[[k]] <- unlist(ColMeans2(Z[[k]]), use.names = FALSE)
    storage.mode(zm[[k]]) <- "double"

    zv[[k]] <- unlist(ColVars2(Z[[k]]), use.names = FALSE)
    storage.mode(zv[[k]]) <- "double"

    loz[[k]] <- Z[[k]]$Info$lo
    storage.mode(loz[[k]]) <- "double"

    upz[[k]] <- Z[[k]]$Info$up
    storage.mode(upz[[k]]) <- "double"

    Id[[k]] <- event2(Z[[k]], ordered = TRUE)
    storage.mode(Id[[k]]) <- "integer"

    InfoP[[k]] <- Z[[k]]$Info$Pattern
    storage.mode(InfoP[[k]]) <- "integer"

    nP[[k]] <- dim(InfoP[[k]])[1L]
    storage.mode(nP[[k]]) <- "integer"

    row.order <- object$Z[[k]]$Info$order
    
    Zipt[[k]] <- Zipt[[k]][row.order, , , , drop = FALSE] #cbind(object$Yipt[[k]][row.order, , 1L, 1L], object$Xipt[[k]][row.order, , 1L, 1L])
    # dim(Zipt[[k]]) <- c(dim(Zipt[[k]]), 1L, 1L)
    # dimnames(Zipt[[k]]) <- list(NULL, response = znames, 
    #                             lambda = dimnames(object$Yipt[[k]])[[3]], 
    #                             rho = dimnames(object$Yipt[[k]])[[4]])
    storage.mode(Zipt[[k]]) <- "double"
    
    # B[[k]] <- object$B[[k]]
    # storage.mode(B[[k]]) <- "double"
    B_n[[k]] <- B[[k]][, , 1L, 1L]
    if(is.vector(B_n[[k]])) B_n[[k]] <- t(B_n[[k]])
    storage.mode(B_n[[k]]) <- "double"
    
    mu[[k]] <- mu[[k]][row.order, , , , drop = FALSE] #cbind(object$mu[[k]][row.order, , 1L,1L], matrix(rep(zm[[k]][id_X], each = n[k]), nrow = n[k]))
    # dim(mu[[k]]) <- c(dim(mu[[k]]), 1L, 1L)
    # dimnames(mu[[k]]) <- list(NULL, response = znames, 
    #                             lambda = dimnames(object$mu[[k]])[[3]], 
    #                             rho = dimnames(object$mu[[k]])[[4]])
    storage.mode(mu[[k]]) <- "double"
    mu_n[[k]] <- mu[[k]][, , 1L, 1L]
    storage.mode(mu_n[[k]]) <- "double"
    
    R[[k]] <- R[[k]][row.order, , , , drop = FALSE]
    storage.mode(R[[k]]) <- "double"
    R_n[[k]] <- R[[k]][, , 1L, 1L]
    storage.mode(R_n[[k]]) <- "double"
  
    # S[[k]] <- if(!X.null) cbind(rbind(object$S[[k]][, , 1L, 1L], object$Sxy[[k]][, , 1L, 1L]), 
    #                             rbind(t(object$Sxy[[k]][, , 1L, 1L]), object$Sxx[[k]][, , 1L, 1L])) else object$S[[k]][, , 1L, 1L]
    # dim(S[[k]]) <- c(dim(S[[k]]), 1L, 1L)
    # dimnames(S[[k]]) <- list(response = znames, response = znames, 
    #                             lambda = dimnames(object$S[[k]])[[3]], 
    #                             rho = dimnames(object$S[[k]])[[4]])
    # storage.mode(S[[k]]) <- "double"
    S_n[[k]] <- S[[k]][, , 1L, 1L]
    storage.mode(S_n[[k]]) <- "double"
    
    # Sgm[[k]] <- if(!X.null) cbind(rbind(object$Sgm[[k]][, , 1L, 1L], object$Sgmxy[[k]][, , 1L, 1L]), 
    #                               rbind(t(object$Sgmxy[[k]][, , 1L, 1L]), object$Sgmxx[[k]][, , 1L, 1L])) else object$Sgm[[k]][, , 1L, 1L]
    # dim(Sgm[[k]]) <- c(dim(Sgm[[k]]), 1L, 1L)
    # dimnames(Sgm[[k]]) <- list(response = znames, response = znames, 
    #                             lambda = dimnames(object$Sgm[[k]])[[3]], 
    #                             rho = dimnames(object$Sgm[[k]])[[4]])
    # storage.mode(Sgm[[k]]) <- "double"
    Sgm_n[[k]] <- Sgm[[k]][, , 1L, 1L]
    storage.mode(Sgm_n[[k]]) <- "double"
    
    # Tht[[k]] <- if(!X.null) cbind(rbind(object$Tht[[k]][, , 1L, 1L], object$Thtxy[[k]][, , 1L, 1L]), 
    #                               rbind(t(object$Thtxy[[k]][, , 1L, 1L]), object$Thtxx[[k]][, , 1L, 1L])) else object$Tht[[k]][, , 1L, 1L]
    # dim(Tht[[k]]) <- c(dim(Tht[[k]]), 1L, 1L)
    # dimnames(Tht[[k]]) <- list(response = znames, response = znames, 
    #                             lambda = dimnames(object$Tht[[k]])[[3]], 
    #                             rho = dimnames(object$Tht[[k]])[[4]])
    # storage.mode(Tht[[k]]) <- "double"
    Tht_n[[k]] <- Tht[[k]][, , 1L, 1L]
    storage.mode(Tht_n[[k]]) <- "double"
    
    # Adj[[k]] <- 1 * (Tht_n[[k]] != 0)
    # dim(Adj[[k]]) <- c(dim(Adj[[k]]), 1L, 1L)
    # storage.mode(Adj[[k]]) <- "integer"
    
    # dfB[[k]] <- object$dfB[[k]]
    # storage.mode(dfB[[k]]) <- "integer"
    
    # dfTht[[k]] <- object$dfTht[[k]]
    # storage.mode(dfTht[[k]]) <- "integer"
    
    # ncomp[[k]] <- object$InfoStructure$ncomp
    # storage.mode(ncomp[[k]]) <- "integer"
    
    Ck[[k]] <- object$InfoStructure$Ck[[k]][, 1L, 1L]
    storage.mode(Ck[[k]]) <- "integer"
    
    pk[[k]] <- object$InfoStructure$pk[[k]][, 1L, 1L]
    storage.mode(pk[[k]]) <- "integer"
    
    Zipt_lo[[k]] <- zm[[k]] - 3 * sqrt(zv[[k]])
    storage.mode(Zipt_lo[[k]]) <- "double"
    
    Zipt_up[[k]] <- zm[[k]] + 3 * sqrt(zv[[k]])
    storage.mode(Zipt_up[[k]]) <- "double"
    
    Zipt_n[[k]] <- Zipt[[k]][, , 1L, 1L]
    storage.mode(Zipt_n[[k]]) <- "double"
    
    T1o[[k]] <- sapply(id_Z, function(j) { id <- Id[[k]][, j] == 0; sum(Zmat[[k]][id, j]) })
    storage.mode(T1o[[k]]) <- "double"
    
    T2o[[k]] <- 0 * S_n[[k]]
    for(i in id_Z) {
      for(j in i:dim_Z) {
        if(any(id <- Id[[k]][, i] == 0 & Id[[k]][, j] == 0)) {
          T2o[[k]][i, j] <- sum(Zmat[[k]][id, i] * Zmat[[k]][id, j])
          T2o[[k]][j, i] <- T2o[[k]][i, j]
        }
      }
    }
    storage.mode(T2o[[k]]) <- "double"
    
    T2[[k]] <- 0 * S_n[[k]]
    storage.mode(T2[[k]]) <- "double"
    
    rho1mat[[k]] <- tp.min*weights.Tht[[k]]
  }
  
  nit.tot <- object$nit
  storage.mode(nit.tot) <- "integer"
  nit <- double(2L)
  cnnctd <- object$connected[, 1L, 1L]
  conv <- integer(1)
  subrout <- integer(1)
  offset <- object$offset
  m_offset <- sapply(offset, mean)
  
  rho2mat <- matrix(tp.min, nrow = p+q, ncol = p+q)
  # if(!X.null) rho2mat <- as.matrix(Matrix:::bdiag(rho2mat, rho2.x * matrix(1, nrow = q, ncol = q)))
  if(penalty != "fused") diag(rho2mat) <- 0
  
  ##### identify unconnected elements, and get blocks #####
  unconnected <- !cnnctd; connected <- cnnctd
  blocklist <- sapply(which(pk[[1L]][pk[[1L]] != 0] > 1), function(.bl) which(Ck[[1L]] == .bl), simplify = FALSE)
  # blocklist <- if(!X.null) list(id_Y, id_X) else list(id_Y)
  
  U <- outer(1:p, 1:p, "<")
  if(!X.null) U2 <- outer(1:q, 1:q, "<")
  
  ##### starting EM algorithm #####
  if(trace == 2) {
    if(!X.null){
      cat("\n*************************************************************************\n",
          "Fitting jcglasso model number ", 
          formatC(nrho, digits = 0, width = 5, format = "d"),
          "\n\t\t       nu = ", 
          formatC(tp.min, digits = 6, width = 9, format = "f"),
          "\n\t\t      rho = ",
          formatC(tp.min, digits = 6, width = 9, format = "f"),
          "\n\t\t   lambda = ",
          formatC(tp.min, digits = 6, width = 9, format = "f"),
          "\n\t\t    alpha = ",
          formatC(alpha, digits = 3, width = 6, format = "f"),
          "\n")
    } 
    else {
      cat("\n*************************************************************************\n",
          "Fitting jcglasso model number ", 
          formatC(nrho, digits = 0, width = 5, format = "d"),
          "\n\t\t      rho = ",
          formatC(tp.min, digits = 6, width = 9, format = "f"),
          "\n\t\t    alpha = ",
          formatC(alpha, digits = 3, width = 6, format = "f"),
          "\n")
    }
  }
  
  for(ii in 1:maxit.em) {
    B_o <- B_n
    Tht_o <- Tht_n
    
    ##### computing E step #####
    for(k in seq_len(K)) {
      temp <- .Fortran(cglasso:::C_e_step_v2, n = n[k], p = dim_Z, Y = Zmat[[k]], 
                       lo = loz[[k]], up = upz[[k]], nP = nP[[k]], 
                       InfoP = InfoP[[k]], T1o = T1o[[k]], T2o = T2o[[k]], 
                       mu_n = mu_n[[k]], Sgm = Sgm_n[[k]], Tht = Tht_n[[k]], 
                       Yipt_lo = Zipt_lo[[k]], Yipt_up = Zipt_up[[k]], 
                       Yipt_n = Zipt_n[[k]], T1 = T1, T2 = T2[[k]], conv = conv)
      if(temp$conv != 0) { subrout <- 1; stop("error in E step!") }
      if(trace == 2) cat("\n\tE-step for group ",
                         formatC(k, digits = 0, width = 3, format = "d"), 
                         ": completed!")
      T2[[k]] <- temp$T2
      zm[[k]] <- temp$T1 / n[k]
      Zipt_n[[k]] <- temp$Yipt_n
      B_n[[k]][1L, ] <- zm[[k]][id_Y] - m_offset[[k]]
      mu_n[[k]][, id_Y] <- outer(offset[[k]], B_n[[k]][1L, ], FUN = "+")
      R_n[[k]] <- Zipt_n[[k]][, id_Y] - mu_n[[k]][, id_Y]
      if(!X.null) {
        S_n[[k]][id_X, id_X] <- T2[[k]][id_X, id_X] / n[k]
        S_n[[k]][id_X, id_Y] <- crossprod(Zipt_n[[k]][, id_X], R_n[[k]]) / n[k]
        S_n[[k]][id_Y, id_X] <- t(S_n[[k]][id_X, id_Y]) #xtr_n[[k]]
        mu_n[[k]][, id_X] <- rep(zm[[k]][id_X], each = n[k])
      }
    }
    
    ##### fitting multilasso model #####
    if(!X.null) {
      if(trace == 2) cat("\n\tM-step:\n\t       fitting multivariate sparse group lasso model with ",
                         "\n\t       lambda = ", formatC(tp.min, digits = 6, width = 9, format = "f"), 
                         "and alpha = ", formatC(alpha, digits = 3, width = 4, format = "f"), 
                         "\n\n")
      # tempmul <- admm.iters.3(p, q, ym, xm, xtx_n, xtr_n, B_n, Tht_n, weights.B, lambda[h], lambda2[h], maxit.bcd, 
      #                         thr.bcd, rho = 2, rho.increment = 1, 0L)

      tempmul <- admm.iters.3(p, q, 
                              ym = lapply(zm, function(x) x[id_Y]), 
                              xm = lapply(zm, function(x) x[id_X]), 
                              xtx = lapply(S_n, function(x) x[id_X, id_X]), 
                              xtr = lapply(S_n, function(x) x[id_X, id_Y]), B = B_n, 
                              Tht = lapply(Tht_n, function(x) x[id_Y, id_Y]), wB = weights.B, 
                              tp.min, tp.min, maxit.bcd, thr.bcd, weights,
                              rho = 2, rho.increment = 1, trace = 0L)
      
      if(tempmul$conv != 0) {
        conv <- tempmul$conv
        subrout <- 2
        stop("error in multilasso step!")
      }
      if(trace == 2) {
        cat("\t\tADMM step\t||B_new(k) - B_old(k)||_1/||B_old(k)||_1\n")
        cat("\t\t   ", formatC(tempmul$nit[1], digits = 0, width = 5, format = "d"),
            "\t\t\t", formatC(tempmul$diff, digits = 8, width = 10, format = "f"),"\n")
      }
      for(k in seq_k){
        B_n[[k]] <- tempmul$B[[k]]
        dfB[[k]][, 1L, 1L] <- tempmul$df[[k]]
        nit[2] <- nit[2] + tempmul$nit[1]
        
        mu_n[[k]][, id_Y] <- outer(offset[[k]], B_n[[k]][1L, ], FUN = "+")
        mu_n[[k]][, id_Y] <- mu_n[[k]][, id_Y] + Zipt_n[[k]][, id_X] %*% B_n[[k]][-1L, , drop = FALSE]
        R_n[[k]] <- Zipt_n[[k]][, id_Y] - mu_n[[k]][, id_Y]
        YM <- crossprod(Zipt_n[[k]][, id_Y], mu_n[[k]][, id_Y])
        # temp$T2 == crossprod(Yipt_n[[k]])
        S_n[[k]][id_Y, id_Y] <- (T2[[k]][id_Y, id_Y] + crossprod(mu_n[[k]][, id_Y]) - YM - t(YM)) / n[k]
        if(covar2corr) S_n[[k]][id_Y, id_Y] <- cov2cor(S_n[[k]][id_Y, id_Y]) #S_n[[k]] / outer(sqrt(diag(S_n[[k]])), sqrt(diag(S_n[[k]])))
        # Sgm_n[[k]] <- S_n[[k]] + rho1mat[[k]] * sign(Tht_n[[k]])
      }
    }
    else{
      for(k in seq_k){
        YM <- crossprod(Zipt_n[[k]][, id_Y], mu_n[[k]][, id_Y])
        # temp$T2 == crossprod(Yipt_n[[k]])
        # T2 <- crossprod(Yipt_n[[k]])
        S_n[[k]][id_Y, id_Y] <- (T2[[k]][id_Y, id_Y] + crossprod(mu_n[[k]][, id_Y]) - YM - t(YM)) / n[k]
        if(covar2corr) S_n[[k]][id_Y, id_Y] <- cov2cor(S_n[[k]][id_Y, id_Y])
      }
    }
    
    whole.S <- lapply(Sgm_n, function(x) 0 * x)
    whole.theta <- lapply(Tht_n, function(x) 0 * x)
    
    ##### computing theta and S values for unconnected components #####
    if(sum(unconnected) != 0) {
      if(covar2corr) {
        for(k in seq_k) {
          diag(whole.theta[[k]])[unconnected] <- rep(1, sum(unconnected))
          diag(whole.S[[k]])[unconnected] <- rep(1, sum(unconnected))
        }
      }
      else {
        S.uncon <- lapply(S_n, function(.S_n) diag(.S_n)[unconnected])
        for(k in seq_k) {
          diag(whole.theta[[k]])[unconnected] <- 1 / S.uncon[[k]]
          diag(whole.S[[k]])[unconnected] <- S.uncon[[k]]
        }
        rm(list = c("S.uncon"))
      }
    }
    
    ##### computing theta and S values for connected components #####
    if(trace == 2) {
      if(X.null) {
        cat("\n")
        cat("\tM-step:\n\t       fitting joint glasso model with",
            "\n\t       rho = ", formatC(tp.min, digits = 6, width = 9, format = "f"), 
            "and alpha = ", formatC(alpha, digits = 3, width = 4, format = "f"), 
            "\n\n")
      } 
      else {
        cat("\tM-step:\n\t       fitting joint glasso model with",
            "\n\t       nu = ", formatC(tp.min, digits = 6, width = 9, format = "f"), 
            ", rho = ", formatC(tp.min, digits = 6, width = 9, format = "f"), 
            "and alpha = ", formatC(alpha, digits = 3, width = 4, format = "f"), 
            "\n\n")
      }
    }
    if(length(blocklist) > 0){
      for(i in seq_len(length(blocklist))) {
        
        # the variables in the block
        bl <- blocklist[[i]] 
        
        # run the admm algorithm on the block:
        Thetabl <- admm.iters.2(S = lapply(S_n, function(.S_n) .S_n[bl, bl]), 
                                theta = lapply(Tht_n, function(.Tht_n) .Tht_n[bl, bl]), 
                                lam1 = lapply(rho1mat, function(x) x[bl, bl]), 
                                lam2 = rho2mat[bl, bl], 
                                penalty = penalty, weights = weights, rho = 2, rho.increment = 1.0,
                                penalize.diagonal = diagonal, maxiter = maxit.bcd, 
                                tol = thr.bcd, covar2corr = covar2corr, trace = 0L)
        if(trace == 2) {
          cat("\tconnected components number ", 
              formatC(i, digits = 0, width = 4, format = "d"), "\n")
          cat("\t\tADMM step\t||Tht_new(k) - Tht_old(k)||_1/||Tht_old(k)||_1\n")
          cat("\t\t   ", formatC(Thetabl$iter, digits = 0, width = 5, format = "d"),
              "\t\t\t", formatC(Thetabl$diff, digits = 8, width = 10, format = "f"),"\n")
        }
        if(Thetabl$conv) {
          conv <- 1
          subrout <- 2
          stop("error in M step (connected components)!")
        }
        nit[2] <- nit[2] + Thetabl$iters
        
        # update results:
        for(k in seq_k) { 
          Thetabl$Z[[k]][abs(Thetabl$Z[[k]]) <= truncate] <- 0
          whole.theta[[k]][bl, bl] <- Thetabl$Z[[k]]
          whole.S[[k]][bl, bl] <- Thetabl$S[[k]]
        }   
      }
      rm(list = c("Thetabl"))
    }
    Sgm_n <- whole.S
    Tht_n <- whole.theta
    nit[1] <- ii
    
    # dB <- sum(sapply(seq_len(K), function(k) sqrt(sum((B_n[[k]] - B_o[[k]])^2)) / (p + dfB[[k]][p + 1, 1L, 1L])))
    # dSgm <- sum(sapply(seq_len(K), function(k) sqrt(sum((Tht_n[[k]] - Tht_o[[k]])[id_Y, id_Y][U1]^2)) / sum(Tht_n[[k]][id_Y, id_Y][U1] != 0)))
    
    # dB <- sum(sapply(seq_len(K), function(k) sum(abs(B_n[[k]] - B_o[[k]])) / sum(abs(B_o[[k]]))))
    # dSgm <- sum(sapply(seq_len(K), function(k) sum(abs((Tht_n[[k]][id_Y, id_Y] - Tht_o[[k]][id_Y, id_Y]))) / sum(abs(Tht_o[[k]][id_Y, id_Y]))))
    
    dB <- sum(sapply(seq_len(K), function(k) norm(B_n[[k]] - B_o[[k]], type = "F") / (p + dfB[[k]][p + 1, 1L, 1L])))
    dTht <- sum(sapply(seq_len(K), function(k) norm((Tht_n[[k]][id_Y, id_Y] - Tht_o[[k]][id_Y, id_Y]), type = "F") / (sum(Tht_n[[k]][id_Y, id_Y][U] != 0) + p)))
    
    if(!X.null){
      for(k in seq_k) {
        Adj[[k]][, , 1L, 1L] <- 1*(Tht_n[[k]] != 0)
        Thtxx[[k]][, , 1L, 1L] <- Tht_n[[k]][id_X, id_X]
        Tht_n[[k]][id_X, id_X] <- Tht_n[[k]][id_X, id_X] + B_n[[k]][-1L, ] %*% Tht_n[[k]][id_Y, id_Y] %*% t(B_n[[k]][-1L, ])
        Tht_n[[k]][id_X, id_Y] <- -(B_n[[k]][-1L, ] %*% Tht_n[[k]][id_Y, id_Y])
        Tht_n[[k]][id_Y, id_X] <- t(Tht_n[[k]][id_X, id_Y])
      }
    } 
    else {
      Adj[[k]][, , 1L, 1L] <- 1*(Tht_n[[k]] != 0)
    }
    
    if(trace == 2) {
      cat("\n\tChecking convergence criterion (threshold = ", 
          formatC(thr.em, digits = 6, width = 8, format = "f"),
          ")\n\t||B_old - B_new||_F / dfB = ",
          formatC(dB, digits = 7, width = 9, format = "f"),
          "\n\t||Tht_old - Tht_new||_F / dfTht = ",
          formatC(dTht, digits = 7, width = 9, format = "f"), "\n")
      if(max(dB, dTht) <= thr.em) cat("\tConvergence criterion is met!\n\n")
    }
    if(max(dB, dTht) <= thr.em) break
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
          "\n\t         nu = ",
          formatC(tp.min, digits = 6, width = 9, format = "f"),
          "\n\t        rho = ",
          formatC(tp.min, digits = 6, width = 9, format = "f"),
          "\n\t     lambda = ",
          formatC(tp.min, digits = 6, width = 9, format = "f"),
          "\n\t      alpha = ",
          formatC(alpha, digits = 3, width = 6, format = "f"),
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
          formatC(tp.min, digits = 6, width = 9, format = "f"),
          "\n\t      alpha = ",
          formatC(alpha, digits = 3, width = 6, format = "f"),
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
    Zipt[[k]][, , 1L, 1L] <- Zipt_n[[k]]
    B[[k]][, , 1L, 1L] <- B_n[[k]]
    mu[[k]][, , 1L, 1L] <- mu_n[[k]]
    R[[k]][, , 1L, 1L] <- R_n[[k]]
    S[[k]][, , 1L, 1L] <- S_n[[k]]
    Sgm[[k]][, , 1L, 1L] <- Sgm_n[[k]]
    Tht[[k]][, , 1L, 1L] <- Tht_n[[k]]
    ncomp[[k]] <- object$InfoStructure$ncomp[[k]]
    Ck[[k]] <- object$InfoStructure$Ck[[k]]
    pk[[k]] <- object$InfoStructure$pk[[k]]
    # Adj[[k]][, , 1L, 1L] <- 1*(Tht[[k]][, , 1L, 1L] != 0)
    diag(Adj[[k]][, , 1L, 1L]) <- 0
    # dimnames(Adj[[k]])[1:2] <- list(ynames, ynames)
    dfTht[[k]][1L, 1L] <- sum(Adj[[k]][id_Y, id_Y, 1L, 1L][U])
    if(!X.null) {
      Adj[[k]][id_X, id_Y, 1L, 1L] <- 1*(B[[k]][-1, , 1L, 1L] != 0)
      Adj[[k]][id_Y, id_X, 1L, 1L] <- t(Adj[[k]][id_X, id_Y, 1L, 1L])
      dfTht[[k]][1L, 1L] <- dfTht[[k]][1L, 1L] + sum(Adj[[k]][id_X, id_X, 1L, 1L][U2])
    }
    
    row.order <- order(Z[[k]]$Info$order)
    Zipt[[k]][, , 1L, 1L] <- Zipt[[k]][, , 1L, 1L][row.order, , drop = FALSE]
    mu[[k]][, , 1L, 1L] <- mu[[k]][, , 1L, 1L][row.order, , drop = FALSE]
    R[[k]][, , 1L, 1L] <- R[[k]][, , 1L, 1L][row.order, , drop = FALSE]
  }
  nit.tot[, 1L, 1L] <- nit
  cnnctd <- object$connected
  
  names(Zipt) <- names(B) <- names(mu) <- names(R) <- k_lab
  names(S) <- names(Sgm) <- names(Tht) <- k_lab
  names(ncomp) <- names(Ck) <- names(pk) <- k_lab
  names(Adj) <- names(dfTht) <- names(dfB) <- k_lab
  names(Z) <- names(Thtxx) <- k_lab
  
  out <- list(n = n, p = p, q = q, id_X = id_X, id_Y = id_Y, Z = Zmat, 
              Id = Id, nP = nP, InfoP = InfoP, ncomp = ncomp, Ck = Ck, pk = pk, 
              lo = loz, up = upz, zm = zm, zv = zv, wTht = weights.Tht, 
              nrho = nrho, rhoratio = object$rho.min.ratio, rho = tp.min,
              nlambda = nlambda, lambdaratio = object$lambda.min.ratio, lambda = tp.min, 
              nu = tp.min, alpha = alpha, rho2 = tp.min, lambda2 = tp.min, 
              pendiag = diagonal, connected = cnnctd, 
              maxit_em = maxit.em, thr_em = thr.em, maxit_bcd = maxit.bcd, thr_bcd = thr.bcd, 
              Zipt = Zipt, B = B, mu = mu, R = R, S = S, Sgm = Sgm, Tht = Tht, 
              Adj = Adj, Omega = Thtxx, dfB = dfB, dfTht = dfTht,
              nit = nit.tot, conv = conv, offset = offset,
              subrout = subrout, trace = trace)
  
  # out <- list(n = n, p = p, Y = lapply(Zmat, function(x) x[, id_Y]), 
  #             Id = Id, nP = nP, InfoP = InfoP, 
  #             lo = lapply(loz, function(x) x[id_Y]), 
  #             up = lapply(upz, function(x) x[id_Y]), 
  #             pendiag = diagonal, wTht = weights.Tht, 
  #             ym = lapply(zm, function(x) x[id_Y]), 
  #             yv = lapply(zv, function(x) x[id_Y]), 
  #             nlambda = nlambda, lambdaratio = object$lambda.min.ratio, lambda = lambda, 
  #             nrho = nrho, rhoratio = object$rho.min.ratio, rho = rho, rho.x = rho.x, 
  #             alpha = object$alpha, lambda2 = lambda2, rho2 = rho2, connected = cnnctd, 
  #             maxit_em = maxit.em, thr_em = thr.em, maxit_bcd = maxit.bcd, thr_bcd = thr.bcd, 
  #             Xipt = if(!X.null) lapply(Zipt, function(x) x[, id_X, , , drop = FALSE]) else NULL, 
  #             Yipt = lapply(Zipt, function(x) x[, id_Y, , , drop = FALSE]), 
  #             B = B, mu = lapply(mu, function(x) x[, id_Y, , , drop = FALSE]), 
  #             R = R, S = lapply(S, function(x) x[id_Y, id_Y, , , drop = FALSE]), 
  #             Sgm = lapply(Sgm, function(x) x[id_Y, id_Y, , , drop = FALSE]), 
  #             Tht = lapply(Tht, function(x) x[id_Y, id_Y, , , drop = FALSE]), 
  #             Sxx = if(!X.null) lapply(S, function(x) x[id_X, id_X, , , drop = FALSE]) else NULL, 
  #             Sxy = if(!X.null) lapply(S, function(x) x[id_X, id_Y, , , drop = FALSE]) else NULL, 
  #             Sgmxx = if(!X.null) lapply(Sgm, function(x) x[id_X, id_X, , , drop = FALSE]) else NULL, 
  #             Sgmxy = if(!X.null) lapply(Sgm, function(x) x[id_X, id_Y, , , drop = FALSE]) else NULL, 
  #             Thtxx = if(!X.null) lapply(Tht, function(x) x[id_X, id_X, , , drop = FALSE]) else NULL, 
  #             Thtxy = if(!X.null) lapply(Tht, function(x) x[id_X, id_Y, , , drop = FALSE]) else NULL,
  #             Adj_yy = lapply(Adj, function(x) x[id_Y, id_Y, , , drop = FALSE]), 
  #             Adj_xy = if(!X.null) lapply(B, function(x) 1*(x[-1L, , , , drop = FALSE] != 0)) else NULL, 
  #             Adj_xx = if(!X.null) lapply(Adj, function(x) 1*(x[id_X, id_X, , , drop = FALSE] != 0)) else NULL, 
  #             dfB = dfB, dfTht = dfTht, ncomp = ncomp, Ck = Ck, pk = pk, 
  #             nit = nit.tot, conv = conv, offset = offset,
  #             subrout = subrout, trace = trace)
  out$conv <- switch(as.character(out$conv),
                     "-1" = "memory allocation error",
                     "0" = "Ok",
                     "1" = "maximum number of iterations has been exceeded",
                     "2" = "error in E-step",
                     "3" = "matrix inversion failed")
  return(out)
}

##### GoF functions #####
QFun2 <- function(object, mle, verbose = FALSE, ...) {
  # testing class attribute
  if (!inherits(object, c("jcggm", "jcglasso"))) stop(sQuote("QFun"), "function is not available for an object of class ", sQuote(class(object)))
  # Testing 'mle'
  if (missing(mle)) mle <- FALSE
  else {
    if (class(object)[1L] == "jcggm" & !mle) stop("For an object of class ", sQuote("jcggm"), " argument ", sQuote("mle"), " must be equal to TRUE")
  }
  if (!is.vector(mle)) stop(sQuote("mle"), " is not a vector")
  if (length(mle) > 1L) stop(sQuote("mle"), " is not an object of length ", sQuote(1))
  if (!is.logical(mle)) stop(sQuote("mle"), " is not a logical object")
  if (class(object)[1L] == "jcggm" & mle) mle <- FALSE
  n <- object$nobs
  K <- length(n)
  p <- object$nresp
  q <- object$npred
  id_X <- object$InfoStructure$id_X
  id_Y <- object$InfoStructure$id_Y
  Adj <- object$InfoStructure$Adj
  rho2 <- object$rho2
  rho1 <- object$rho
  nrho <- object$nrho
  lambda2 <- object$lambda2
  lambda1 <- object$lambda
  nlambda <- object$nlambda
  X.null <- is.null(object$Xipt)
  
  lambda_lab <- paste0("lambda_", seq_len(nlambda))
  rho_lab <- paste0("rho_", seq_len(nrho))
  value <- df <- dfB <- dfTht <- array(0, dim = c(nlambda, nrho, K), dimnames = list(lambda_lab, rho_lab, NULL))
  value.tot <- df.tot <- dfB.tot <- dfTht.tot <- matrix(0, nrow = nlambda, ncol = nrho, dimnames = list(lambda_lab, rho_lab))
  value.comp <- df.comp <- dfB.comp <- dfTht.comp <- vector(mode = "list", length = K)
  const1 <- (- 0.5 * n * p * log(2 * pi))*0
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
        test1 <- test2 <- vector(mode = "logical", length = K)
        i1 <- ii[hh, 1L]
        j1 <- ii[hh, 2L]
        for(k in seq_len(K)) {
          test1[k] <- all(Adj[[k]][id_Y, id_Y, i, j] == Adj[[k]][id_Y, id_Y, i1, j1])
          if (q == 0L) test2[[k]] <- TRUE
          else test2[k] <- all(Adj[[k]][id_X, id_Y, i, j] == Adj[[k]][id_X, id_Y, i1, j1])
        }
        if (all(test1 & test2)) {
          NEW <- FALSE
          break
        }
      }
      if (NEW) {
        if(!is.null(list(...)$trace)) cat(paste0("Model: ", h, " (lambda.id: ", i, " and rho.id: ", j, ")"))
        out.mle <- jcggm(object, lambda.id = i, rho.id = j, ...)
        for(k in seq_len(K)){
          S <- out.mle$S[[k]][, , 1L, 1L]
          Tht <- out.mle$Tht[[k]][, , 1L, 1L]
          # S <- if(!X.null) cbind(rbind(out.mle$S[[k]][, , 1L, 1L], out.mle$Sxy[[k]][, , 1L, 1L]), 
          #                        rbind(t(out.mle$Sxy[[k]][, , 1L, 1L]), out.mle$Sxx[[k]][, , 1L, 1L])) 
          # else out.mle$S[[k]][, , 1L, 1L]
          # Tht <- if(!X.null) cbind(rbind(out.mle$Tht[[k]][, , 1L, 1L], out.mle$Thtxy[[k]][, , 1L, 1L]), 
          #                          rbind(t(out.mle$Thtxy[[k]][, , 1L, 1L]), out.mle$Thtxx[[k]][, , 1L, 1L])) 
          # else out.mle$Tht[[k]][, , 1L, 1L]
          
          # Thtyy <- out.mle$Tht[[k]][, , 1L, 1L]
          # Syy <- out.mle$S[[k]][, , 1L, 1L]
          dfB[i, j, k] <- out.mle$dfB[[k]][p + 1L, 1L, 1L]
          dfTht[i, j, k] <- out.mle$dfTht[[k]][1L, 1L]
          df[i, j, k] <- dfB[i, j, k] + p + dfTht[i, j, k] + (p + q)
          # if(!is.null(out.mle$Xipt)) {
          #   Thtxx <- out.mle$Thtxx[[k]][, , 1L, 1L]
          #   Sxx <- out.mle$Sxx[[k]][, , 1L, 1L]
          #   Sxy <- out.mle$Sxy[[k]][, , 1L, 1L]
          #   B <- out.mle$B[[k]][-1, , 1L, 1L]
          #   uno <- determinant(Thtxx)$modulus - sum(Sxx * Thtxx)
          #   due <- determinant(Thtyy)$modulus - sum(Thtyy * (Syy - 2*crossprod(B, Sxy) + crossprod(B, Sxx) %*% B))
          #   value[i, j, k] <- const1[k] + const2[k] * (uno + due)
          #   df[i, j, k] <- df[i, j, k] + sum(Thtxx[upper.tri(Thtxx, diag = TRUE)] != 0)
          # } else value[i, j, k] <- const1[k] + const2[k] * (determinant(Thtyy)$modulus - sum(Syy * Thtyy))
          value[i, j, k] <- const1[k] + const2[k] * (determinant(Tht)$modulus - sum(S * Tht))
        }
      } 
      else {
        for(k in seq_len(K)) {
          value[i, j, k] <- value[i1, j1, k]
          dfB[i, j, k] <- dfB[i1, j1, k]
          dfTht[i, j, k] <- dfTht[i1, j1, k]
          df[i, j, k] <- df[i1, j1, k]
        }
      }
    } 
    else {
      for(k in seq_len(K)) {
        S <- object$S[[k]][, , i, j]
        Tht <- object$Tht[[k]][, , i, j]
        # S <- if(!X.null) cbind(rbind(object$S[[k]][, , i, j], object$Sxy[[k]][, , i, j]), 
        #                        rbind(t(object$Sxy[[k]][, , i, j]), object$Sxx[[k]][, , i, j])) 
        # else object$S[[k]][, , i, j]
        # Tht <- if(!X.null) cbind(rbind(object$Tht[[k]][, , i, j], object$Thtxy[[k]][, , i, j]), 
        #                          rbind(t(object$Thtxy[[k]][, , i, j]), object$Thtxx[[k]][, , i, j])) 
        # else object$Tht[[k]][, , i, j]
        
        # Thtyy <- object$Tht[[k]][, , i, j]
        # Syy <- object$S[[k]][, , i, j]
        dfB[i, j, k] <- object$dfB[[k]][p + 1L, i, j]
        dfTht[i, j, k] <- object$dfTht[[k]][i, j]
        df[i, j, k] <- dfB[i, j, k] + p + dfTht[i, j, k] + (p + q)
        # if(!is.null(object$Xipt)) {
        #   Thtxx <- object$Thtxx[[k]][, , i, j]
        #   Sxx <- object$Sxx[[k]][, , i, j]
        #   Sxy <- object$Sxy[[k]][, , i, j]
        #   B <- object$B[[k]][-1, , i, j]
        #   uno <- determinant(Thtxx)$modulus - sum(Sxx * Thtxx)
        #   due <- determinant(Thtyy)$modulus - sum(Thtyy * (Syy - 2*crossprod(B, Sxy) + crossprod(B, Sxx) %*% B))
        #   value[i, j, k] <- const1[k] + const2[k] * (uno + due)
        #   df[i, j, k] <- df[i, j, k] + sum(Thtxx[upper.tri(Thtxx, diag = TRUE)] != 0)
        # } else value[i, j, k] <- const1[k] + const2[k] * (determinant(Thtyy)$modulus - sum(Syy * Thtyy))
        value[i, j, k] <- const1[k] + const2[k] * (determinant(Tht)$modulus - sum(S * Tht))
      }
    }
  }
  for(k in seq_len(K)) {
    value.tot <- value.tot + value[, , k]
    value.comp[[k]] <- value[, , k]
    dfB.tot <- dfB.tot + dfB[, , k]
    dfB.comp[[k]] <- dfB[, , k]
    dfTht.tot <- dfTht.tot + dfTht[, , k]
    dfTht.comp[[k]] <- dfTht[, , k]
    df.tot <- df.tot + df[, , k]
    df.comp[[k]] <- df[, , k]
  }
  if(verbose) close(pb)
  out <- list(value = value.tot, df = df.tot, dfB = dfB.tot, dfTht = dfTht.tot,
              value.comp = value.comp, df.comp = df.comp, dfB.comp = dfB.comp, dfTht.comp = dfTht.comp, 
              n = n, p = p, q = q, rho2 = rho2, rho1 = rho1, nrho = nrho,
              lambda2 = lambda2, lambda1 = lambda1, nlambda = nlambda, 
              nu = object$nu, alpha = object$alpha, model = object$model)
  class(out) <- "QFun2"
  out
}

print.QFun2 <- function (x, digits = 3L, ...){
  nrho <- x$nrho
  rho1 <- x$rho1
  rho2 <- x$rho2
  nlambda <- x$nlambda
  lambda1 <- x$lambda1
  lambda2 <- x$lambda2
  rho.x <- x$nu
  alpha <- x$alpha
  K <- length(x$value)
  q <- x$q
  
  if (nrho == 1L | nlambda == 1L) {
    value <- drop(x$value)
    df <- drop(x$df)
  }
  else {
    value <- as.vector(t(x$value))
    df <- as.vector(t(x$df))
  }
  if (nrho > 1L & nlambda > 1L) {
    lambda1 <- rep(lambda1, each = nrho)
    lambda2 <- rep(lambda2, each = nrho)
    rho1 <- rep(rho1, times = nlambda)
    rho2 <- rep(rho2, times = nlambda)
  }
  
  dots <- list(...)
  if (is.null(dots$print.gap)) dots$print.gap <- 2L
  if (is.null(dots$quote)) dots$quote <- FALSE
  if (is.null(dots$row.names)) dots$row.names <- FALSE
  tbl <- data.frame(nu = if(is.null(rho.x)) rho1+rho2 else rho.x, 
                    rho = rho1+rho2, lambda = lambda1+lambda2, alpha = alpha, df = df, value = value)
  names(tbl)[6L] <- "Q-Values"
  if (q == 0L) tbl <- tbl[, -c(1L,3L), drop = FALSE]
  
  cat("\nQ-values of the fitted", sQuote(x$model[1L]), ifelse((nrho == 1L & nlambda == 1L), "model", "models"))
  cat("\n\nDetails:\n")
  if (nlambda <= nrho) f <- rep(seq_len(nlambda), each = nrho)
  else f <- rep(seq_len(nrho), times = nlambda)
  
  tbl.list <- split(tbl, f = f, drop = FALSE)
  do.call(function(...) print.listof(tbl.list, digits = digits, ...), dots)
  
  invisible(tbl)
}

AIC.jcglasso <- function(object, k = 2, mle, components = FALSE, Qfun = NULL, ...){
  K <- length(object$nobs)
  if(length(k) == 1) k <- rep(2, K)
  if(any(k <= 0)) stop(sQuote("k"), " is not a positive value")
  type <- ifelse(all(k == 2), "AIC", "GoF")
  
  out_QFun <- if(is.null(Qfun)) QFun2(object, mle, ...) else Qfun
  nrho <- out_QFun$nrho
  nlambda <- out_QFun$nlambda
  if(components) {
    value <-  out_QFun$value.comp
    df <- out_QFun$df.comp
    val <- vector(mode = "list", length = K)
    for(h in seq_len(K)) {
      val[[h]] <- -2 * value[[h]] + k[h] * df[[h]]
    }
    out <- list(value_gof = val, df = df, dfB = out_QFun$dfB.comp, dfTht = out_QFun$dfTht.comp, value = value,
                n = out_QFun$n, p = out_QFun$p, q = out_QFun$q, rho2 = out_QFun$rho2,
                rho1 = out_QFun$rho1, nrho = out_QFun$nrho, lambda2 = out_QFun$lambda2,
                lambda1 = out_QFun$lambda1, nlambda = out_QFun$nlambda, nu = out_QFun$nu, alpha = out_QFun$alpha,
                type = type, model = out_QFun$model, components = components)
  }
  else {
    value <-  out_QFun$value
    df <- sapply(seq_len(K), function(h) k[h] * out_QFun$df.comp[[h]])
    if(is.vector(df)) df <- t(df)
    df <- matrix(rowSums(df), nrow = nlambda, ncol = nrho, dimnames = dimnames(value))
    val <- -2 * value + df
    out <- list(value_gof = val, df = out_QFun$df, dfB = out_QFun$dfB, dfTht = out_QFun$dfTht, value = value,
                n = out_QFun$n, p = out_QFun$p, q = out_QFun$q, rho2 = out_QFun$rho2,
                rho1 = out_QFun$rho1, nrho = out_QFun$nrho, lambda2 = out_QFun$lambda2,
                lambda1 = out_QFun$lambda1, nlambda = out_QFun$nlambda, alpha = out_QFun$alpha,
                type = type, model = out_QFun$model, components = components)
  }
  
  class(out) <- "GoF2"
  out
}

BIC.jcglasso <- function(object, g = 0, type, mle, components = FALSE, Qfun = NULL, ...){
  # type == "FD"
  # ebic measure as proposed in Foygel and Drton (2010)
  #            val <- -2 * value + dfTht * (log(n) + 4 * g * log(p))
  # type == "CC"
  # generalization based on the ebic measure proposed in Chen and Chen (2008)
  # See also Chen and Chen (2012) Statistica Sinica 22, pg. 555-574
  
  # testing 'g'
  if (!is.vector(g)) stop(sQuote("g"), " is not a vector")
  if (length(g) != 1) stop(sQuote("g"), " is not a vector of length ", sQuote(1))
  if ((g < 0) | (g > 1)) stop(sQuote("g"), " does not belong to the closed interval [0, 1]")
  if (g == 0) {
    # classical BIC measure
    out <- AIC(object, k = log(object$nobs), mle, components = components, Qfun = Qfun, ...)
    out$type <- "BIC"
  } 
  else {
    out_QFun <- if(is.null(Qfun)) QFun2(object, mle, ...) else Qfun
    n <- out_QFun$n
    p <- out_QFun$p
    q <- out_QFun$q
    K <- length(n)
    k <- if (type == "FD") log(n) + 4 * g * log(p) else log(n) + 2 * g * log(q)
    nrho <- out_QFun$nrho
    nlambda <- out_QFun$nlambda
    if(components) {
      value <-  out_QFun$value.comp
      df <- out_QFun$df.comp
      val <- vector(mode = "list", length = K)
      for(h in seq_len(K)) {
        val[[h]] <- -2 * value[[h]] + k[h] * df[[h]]
      }
      out <- list(value_gof = val, df = df, dfB = out_QFun$dfB.comp, dfTht = out_QFun$dfTht.comp, value = value,
                  n = out_QFun$n, p = out_QFun$p, q = out_QFun$q, rho2 = out_QFun$rho2,
                  rho1 = out_QFun$rho1, nrho = out_QFun$nrho, lambda2 = out_QFun$lambda2,
                  lambda1 = out_QFun$lambda1, nlambda = out_QFun$nlambda, nu = out_QFun$nu, alpha = out_QFun$alpha,
                  type = type, model = out_QFun$model, components = components)
    }
    else {
      value <-  out_QFun$value
      df <- sapply(seq_len(K), function(h) k[h] * out_QFun$df.comp[[h]])
      if(is.vector(df)) df <- t(df)
      df <- matrix(rowSums(df), nrow = nlambda, ncol = nrho, dimnames = dimnames(value))
      if (missing(type)) type <- ifelse(q == 0L, "FD", "CC")
      else {
        if (!is.vector(type)) stop(sQuote("type"), "is not a vector")
        if (length(type) != 1L) stop(sQuote("type"), "is not a vector of length ", sQuote("1"))
        if (!is.character(type)) stop(sQuote("type"), "is not an object of type ", sQuote("character"))
        if (!is.element(type, c("FD", "CC")))
          stop(sQuote(type), " is not available. Please, set ", sQuote("type"), " argument equal to ", sQuote("FD"), " or ", sQuote("CC"))
        if (type == "CC" & q == 0L) stop("measure proposed in Chen and Chen (2008, 2012) can not be computed because q = 0")
      }
      val <- -2 * value + df
      type <- paste0("eBIC_", type)
      out <- list(value_gof = val, df = out_QFun$df, dfB = out_QFun$dfB, dfTht = out_QFun$dfTht, value = value,
                  n = out_QFun$n, p = out_QFun$p, q = out_QFun$q, rho2 = out_QFun$rho2,
                  rho1 = out_QFun$rho1, nrho = out_QFun$nrho, lambda2 = out_QFun$lambda2,
                  lambda1 = out_QFun$lambda1, nlambda = out_QFun$nlambda, nu = out_QFun$nu, alpha = out_QFun$alpha,
                  type = type, model = out_QFun$model, components = components)
    }
  }
  class(out) <- "GoF2"
  out
}

print.GoF2 <- function (x, digits = 3L, ...){
  type <- x$type
  nrho <- x$nrho
  rho1 <- x$rho1
  rho2 <- x$rho2
  nlambda <- x$nlambda
  lambda1 <- x$lambda1
  lambda2 <- x$lambda2
  rho.x <- x$nu
  alpha <- x$alpha
  q <- x$q
  
  if (nrho > 1L & nlambda > 1L) {
    lambda1 <- rep(lambda1, each = nrho)
    lambda2 <- rep(lambda2, each = nrho)
    rho1 <- rep(rho1, times = nlambda)
    rho2 <- rep(rho2, times = nlambda)
  }
  if(x$components){
    K <- length(x$value_gof)
    if (nrho == 1L | nlambda == 1L) {
      val <- sapply(x$value_gof, function(y) drop(y))
      df <- sapply(x$df, function(y) as.vector(y))
    }
    else {
      val <- sapply(x$value_gof, function(y) as.vector(t(y)))
      df <- sapply(x$df, function(y) as.vector(t(y)))
    }
    if(is.vector(df)) df <- t(df)
    if(is.vector(val)) val <- t(val)
    tbl <- data.frame(nu = if(is.null(rho.x)) rho1+rho2 else rho.x, 
                      rho = rho1+rho2, lambda = lambda1+lambda2, alpha = alpha, 
                      df = df, df.Tot = rowSums(df), val = val, val.Tot = rowSums(val))
    names(tbl)[5L:(5L + K - 1L)] <- paste0("df.Comp", seq_len(K))
    names(tbl)[-c(1L:(5L + K), ncol(tbl))] <- paste0(type, ".Comp", seq_len(K))
  }
  else{
    if (nrho == 1L | nlambda == 1L) {
      val <- drop(x$value_gof)
      df <- drop(x$df)
    }
    else {
      val <- as.vector(t(x$value_gof))
      df <- as.vector(t(x$df))
    }
    tbl <- data.frame(nu = if(is.null(rho.x)) rho1+rho2 else rho.x, 
                      rho = rho1+rho2, lambda = lambda1+lambda2, alpha = alpha, df = df, val = val)
    names(tbl)[6L] <- type
  }
  
  cat("\nSequence of", sQuote(type), "values of the fitted", sQuote(x$model[1L]), ifelse(nrho == 1L & nlambda == 1L, "model", "models"))
  cat("\n\nDetails:\n")
  dots <- list(...)
  if (is.null(dots$print.gap)) dots$print.gap <- 2L
  if (is.null(dots$quote)) dots$quote <- FALSE
  if (is.null(dots$row.names)) dots$row.names <- FALSE
  if (q == 0L) tbl <- tbl[, -c(1L,3L), drop = FALSE]
  
  if (nlambda <= nrho) f <- rep(seq_len(nlambda), each = nrho)
  else f <- rep(seq_len(nrho), times = nlambda)
  tbl.list <- split(tbl, f = f, drop = FALSE)
  do.call(function(...) print.listof(tbl.list, digits = digits, ...), dots)
  
  invisible(tbl)
}

select.jcglasso <- function(object, GoF = AIC, rho.id, lambda.id, ...){
  if (!is.element(class(GoF), c("function", "GoF2")))
    stop (sQuote("GoF"), " is not either a goodness-of-fit function (AIC or BIC) neither an object of class ", sQuote("GoF"))
  n <- object$nobs
  K <- length(n)
  p <- object$nresp
  q <- object$npred
  nrho <- object$nrho
  nlambda <- object$nlambda
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
    val <- which(GoF$value_gof == min(GoF$value_gof), arr.ind = TRUE)
    if (is.matrix(val)) val <- val[dim(val)[1L], ]
    lambda.id <- val[1L]
    rho.id <- val[2L]
  }
  else {
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
  }
  # reshape cglasso object
  for(k in seq_len(K)) {
    object$Zipt[[k]] <- object$Zipt[[k]][ , , lambda.id, rho.id, drop = FALSE]
    object$B[[k]] <- object$B[[k]][ , , lambda.id, rho.id, drop = FALSE]
    object$mu[[k]] <- object$mu[[k]][ , , lambda.id, rho.id, drop = FALSE]
    object$R[[k]] <- object$R[[k]][ , , lambda.id, rho.id, drop = FALSE]
    object$S[[k]] <- object$S[[k]][ , , lambda.id, rho.id, drop = FALSE]
    object$Sgm[[k]] <- object$Sgm[[k]][ , , lambda.id, rho.id, drop = FALSE]
    object$Tht[[k]] <- object$Tht[[k]][ , , lambda.id, rho.id, drop = FALSE]
    if(q > 0) object$Omega[[k]] <- object$Omega[[k]][ , , lambda.id, rho.id, drop = FALSE]
    object$dfB[[k]] <- object$dfB[[k]][ , lambda.id, rho.id, drop = FALSE]
    object$dfTht[[k]] <- object$dfTht[[k]][lambda.id, rho.id, drop = FALSE]
    object$InfoStructure$Adj[[k]] <- object$InfoStructure$Adj[[k]][ , , lambda.id, rho.id, drop = FALSE]
    object$InfoStructure$ncomp[[k]] <- object$InfoStructure$ncomp[[k]][lambda.id, rho.id, drop = FALSE]
    object$InfoStructure$Ck[[k]] <- object$InfoStructure$Ck[[k]][ , lambda.id, rho.id, drop = FALSE]
    object$InfoStructure$pk[[k]] <- object$InfoStructure$pk[[k]][ , lambda.id, rho.id, drop = FALSE]
  }
  object$nit <- object$nit[ , lambda.id, rho.id, drop = FALSE]
  object$connected <- object$connected[ , lambda.id, rho.id, drop = FALSE]
  
  object$nlambda <- 1L
  object$lambda.min.ratio <- 1
  object$lambda <- object$lambda[lambda.id]
  object$lambda2 <- object$lambda2[lambda.id]
  object$nrho <- 1L
  object$rho.min.ratio <- 1
  object$rho <- object$rho[rho.id]
  object$rho2 <- object$rho2[rho.id]
  
  object
}

plot.GoF2 <- function (x, add.line = TRUE, arg.line = list(lty = 2L, lwd = 2L, col = "red"), 
                       add.text = FALSE, arg.text = list(side = 3L), 
                       arg.points = list(pch = 2L), ...) {
  if (!is.vector(add.line)) 
    stop(sQuote("add.line"), " is not a vector")
  if (length(add.line) != 1L) 
    stop(sQuote("add.line"), " is not a vector of length ", sQuote("1"))
  if (!is.logical(add.line)) 
    stop(sQuote("add.line"), " is not an object of type ", sQuote("logical"))
  if (!is.list(arg.line)) 
    stop(sQuote("arg.line"), " is not a list")
  if (!is.vector(add.text)) 
    stop(sQuote("add.text"), " is not a vector")
  if (length(add.text) != 1L) 
    stop(sQuote("add.text"), " is not a vector of length ", sQuote("1"))
  if (!is.logical(add.text)) 
    stop(sQuote("add.text"), " is not an object of type ", sQuote("logical"))
  if (!is.list(arg.text)) 
    stop(sQuote("arg.text"), " is not a list")
  nrho <- x$nrho
  rho1 <- x$rho1
  rho2 <- x$rho2
  nlambda <- x$nlambda
  lambda1 <- x$lambda1
  lambda2 <- x$lambda2
  alpha <- x$alpha
  lambda <- lambda1 + lambda2
  rho <- rho1 + rho2
  if (nrho == 1L & nlambda == 1L) 
    stop("plot method is not available because ", sQuote("nlambda = 1"), 
         " and ", sQuote("nrho = 1"))
  dots <- list(...)
  if (is.null(arg.line$lty)) 
    arg.line$lty <- 2L
  if (is.null(arg.line$lwd)) 
    arg.line$lwd <- 2L
  if (is.null(arg.line$col)) 
    arg.line$col <- "red"
  if (is.null(arg.text$side)) 
    arg.text$side <- 3L
  if (is.null(arg.points$pch)) 
    arg.points$pch <- 2L
  if (nrho == 1L | nlambda == 1L) {
    val <- drop(x$value_gof)
    pen <- if (nlambda == 1L) rho else lambda
    df <- drop(x$df)
    if(is.vector(df)) df <- t(df)
    minval <- which.min(val)
    if (is.null(dots$main)) 
      dots$main <- "Tuning Parameter Selection"
    if (is.null(dots$sub)) {
      dots$sub <- switch(x$type, 
                         "AIC" = "Akaike Information Criterion", 
                         "BIC" = "Bayesian Information Criterion", 
                         "GoF" = "Measure of Goodness of Fit", 
                         "eBIC_FD" = "extended Bayesian Information Criterion",
                         "eBIC_CC" = "extended Bayesian Information Criterion")
    }
    if (is.null(dots$xlab)) 
      dots$xlab <- ifelse(nlambda == 1L, expression(log(rho)), expression(log(lambda)))
    if (is.null(dots$ylab)) 
      dots$ylab <- "Values"
    if (is.null(dots$type)) 
      dots$type <- "b"
    do.call(function(...) plot(x = log(pen), y = val, ...), dots)
    if (add.line) {
      do.call(function(...) abline(v = log(pen)[minval], ...), arg.line)
      if (add.text) {
        if (is.null(arg.text$text)) 
          arg.text$text <- paste0("df = ", df[minval])
        do.call(function(...) mtext(at = log(pen)[minval], ...), arg.text)
      }
    }
  }
  else {
    rho <- rev(rho)
    lambda <- rev(lambda)
    val <- matrix(x$value_gof, nrow = nlambda, ncol = nrho)
    val <- val[nlambda:1, nrho:1]
    minval <- min(val)
    rc.cord <- drop(which(val == minval, arr.ind = TRUE))
    if (is.matrix(rc.cord)) 
      rc.cord <- rc.cord[dim(rc.cord)[1L], ]
    if (is.null(dots$main)) 
      dots$main <- "Tuning Parameters Selection"
    if (is.null(dots$xlab)) 
      dots$xlab <- expression(log(lambda))
    if (is.null(dots$ylab)) 
      dots$ylab <- expression(log(rho))
    if (is.null(dots$sub)) {
      dots$sub <- switch(x$type,
                         "AIC" = "Akaike Information Criterion",
                         "BIC" = "Bayesian Information Criterion",
                         "GoF" = "Measure of Goodness of Fit",
                         "eBIC_FD" = "extended Bayesian Information Criterion",
                         "eBIC_CC" = "extended Bayesian Information Criterion")
    }
    do.call(function(...) contour(x = log(lambda), y = log(rho), z = val, ...), dots)
    if (add.line) {
      do.call(function(...) abline(h = log(rho)[rc.cord[2L]], ...), arg.line)
      do.call(function(...) abline(v = log(lambda)[rc.cord[1L]], ...), arg.line)
    }
    do.call(function(...) points(x = log(lambda)[rc.cord[1L]], y = log(rho)[rc.cord[2L]], ...), arg.points) 
    
  }
}

##### to_graph functions #####
to_graph2 <- function(object, GoF = AIC, rho.id, lambda.id, weighted = FALSE, simplify = TRUE, ...) {
  K <- length(object$nobs)
  q <- object$npred
  pY <- object$nresp
  if (!inherits(object, c("jcggm", "jcglasso"))) stop(sQuote("to_graph"), "function is not available for an object of class ", sQuote(class(object)))
  nrho <- object$nrho
  nlambda <- object$nlambda
  id_X <- object$InfoStructure$id_X
  id_Y <- object$InfoStructure$id_Y
  if (missing(lambda.id) & missing(rho.id)) {
    if (!is.element(class(GoF), c("function", "GoF2")))
      stop (sQuote("GoF"), " is not either a goodness-of-fit function (AIC or BIC) or an object of class ", sQuote("GoF"))
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
    nlambda <- object$nlambda
    lambda.id <- 1L
    nrho <- object$nrho
    rho.id <- 1L
  } 
  else {
    # testing 'pen2.id'
    if (missing(lambda.id)) {
      if (q == 0 | nlambda == 1L) lambda.id <- 1L
      else stop(sQuote("lambda.id"), " is missing")
    } 
    else {
      if (!is.vector(lambda.id)) stop(sQuote("lambda.id"), " is not a vector")
      if (length(lambda.id) != 1L) stop(sQuote("lambda.id"), " is not a vector of length ", sQuote("1"))
      if (any(abs(as.integer(lambda.id) - lambda.id) > 0)) stop(sQuote("lambda.id"), " is not an object of type ", dQuote("integer"))
      if (lambda.id <= 0) stop(sQuote("lambda.id"), " is not a positive integer")
      if (lambda.id > nlambda) stop("some entry in ", sQuote("lambda.id"), " is larger than ", sQuote(nlambda))
    }
    # testing 'rho.id'
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
  if (!is.logical(simplify)) stop(sQuote("simplify"), " is not an object of type ", sQuote("logical"))
  if (!is.logical(weighted)) stop(sQuote("weighted"), " is not an object of type ", dQuote("logical"))
  Adj <- object$InfoStructure$Adj
  nmsY <- colNames2(object$Z[[1]])$Y
  outTht <- lapply(Adj, function(x, pY) {
    # dimnames(x) <- list(nmsY, nmsY, NULL, NULL)
    out <- graph_from_adjacency_matrix(adjmatrix = x[id_Y, id_Y, lambda.id, rho.id], mode = "undirected", diag = FALSE)
    V(out)$type <- rep("response", pY)
    V(out)$color <- NA
    V(out)$frame.color <- NA
    V(out)$size <- 12
    V(out)$label.cex <- 1
    V(out)$label.font <- 2
    V(out)$label.dist <- 0
    V(out)$label.color <- rep("black", pY)
    E(out)$type <- "undirect"
    E(out)$color <- "gray50"
    E(out)$width <- 1
    E(out)$arrow.mode <- 0
    # rmvTht <- degree(out) == 0
    # if (any(rmvTht) & sum(rmvTht) < vcount(out))
    #   out <- delete.vertices(out, rmvTht)
    out
  }, pY = pY)
  
  outB <- outThtxx <- NULL
  if (q > 0L){
    # Adj_xx <- lapply(object$Thtxx, function(x, lambda.id, rho.id) 1*(x[, , lambda.id, rho.id] != 0), 
    #                  lambda.id = lambda.id, rho.id = rho.id)
    outThtxx <- lapply(Adj, function(x, qX) {
      out <- graph_from_adjacency_matrix(adjmatrix = x[id_X, id_X, lambda.id, rho.id], mode = "undirected", diag = FALSE)
      V(out)$type <- rep("response", qX)
      V(out)$color <- NA
      V(out)$frame.color <- NA
      V(out)$size <- 12
      V(out)$label.cex <- 1
      V(out)$label.font <- 2
      V(out)$label.dist <- 0
      V(out)$label.color <- "black"
      E(out)$type <- "undirect"
      E(out)$color <- "gray50"
      E(out)$width <- 1
      E(out)$arrow.mode <- 0
      out
    }, qX = q)
    
    # Adj_xy <- object$InfoStructure$Adj_xy
    nmsX <- colNames2(object$Z[[1]])$X
    for(k in seq_len(K)) {
      V(outTht[[k]])$label.color <- ifelse(apply(Adj[[k]][id_X, id_Y, lambda.id, rho.id, drop = FALSE], 2L, function(x) any(x == 1)), "black", "gray50")
      E(outTht[[k]])$color <- "gray85"
      V(outThtxx[[k]])$label.color <- ifelse(apply(Adj[[k]][id_X, id_Y, lambda.id, rho.id, drop = FALSE], 1L, function(x) any(x == 1)), "darkblue", "gray50")
      E(outThtxx[[k]])$color <- "gray85"
    }
    outB <- sapply(1:K, function(k, pX, pY, nmsX, nmsY, lambda.id, rho.id) {
      out <- Adj[[k]][id_X, id_Y, lambda.id, rho.id, drop = FALSE]
      eB <- data.frame(from = rep(nmsX, pY), to = rep(nmsY, each = pX))[c(out) != 0,]
      vB <- data.frame(id = c(nmsY, nmsX))
      outB <- graph_from_data_frame(d = eB, directed = TRUE, vertices = vB)
      V(outB)$type <- rep(c("response", "predictor"), c(pY, pX))
      V(outB)$color <- rep(c(NA, NA), c(pY, pX))
      V(outB)$frame.color <- NA
      V(outB)$size <- 12
      V(outB)$label.cex <- 1
      V(outB)$label.font <- 2
      V(outB)$label.dist <- 0
      V(outB)$label.color <- rep(c(NA, "darkblue"), c(pY, pX))
      E(outB)$type <- "direct"
      E(outB)$color <- "red"
      E(outB)$width <- .2
      E(outB)$arrow.mode <- 2
      outB
    }, pX = q, pY = pY, nmsX = nmsX, nmsY = nmsY, lambda.id = lambda.id, rho.id = rho.id, simplify = FALSE)
  }
  if (weighted) {
    for(k in seq_len(K)) {
      Tht <- cov2cor(object$Tht[[k]][id_Y, id_Y, lambda.id, rho.id])
      weight <- Tht[lower.tri(Tht)][Tht[lower.tri(Tht)] != 0]
      E(outTht[[k]])$lty <- ifelse(weight > 0, "solid", "dashed")
      E(outTht[[k]])$width <- as.numeric(as.character(cut(abs(weight), 
                                             breaks = c(0, .0001, .01, 1), 
                                             labels = c(.5, 1.5, 2.5))))
      if(q > 0L) {
        Tht <- cov2cor(object$Omega[[k]][, , lambda.id, rho.id])
        weight <- Tht[lower.tri(Tht)][Tht[lower.tri(Tht)] != 0]
        E(outThtxx[[k]])$lty <- ifelse(weight > 0, "solid", "dashed")
        E(outThtxx[[k]])$width <- as.numeric(as.character(cut(abs(weight), 
                                                 breaks = c(0, .0001, .01, 1), 
                                                 labels = c(.5, 1.5, 2.5))))
        
        B <- object$B[[k]][-1L, , lambda.id, rho.id]
        weight <- B[B != 0]
        E(outB[[k]])$lty <- ifelse(weight > 0, "solid", "dashed")
        E(outB[[k]])$width <- as.numeric(as.character(cut(abs(weight), 
                                             breaks = c(0, .0001, .01, 1), 
                                             labels = c(.5, 1.5, 2.5))))
      }
    }
  }
  if (simplify) {
    for(k in seq_len(K)) {
      rmvTht <- degree(outTht[[k]]) == 0
      if (any(rmvTht)& sum(rmvTht) < vcount(outTht[[k]])) outTht[[k]] <- delete.vertices(outTht[[k]], rmvTht)
      if (q != 0L) {
        rmvThtxx <- degree(outThtxx[[k]]) == 0
        if (any(rmvThtxx)& sum(rmvThtxx) < vcount(outThtxx[[k]])) outThtxx[[k]] <- delete.vertices(outThtxx[[k]], rmvThtxx)
        
        rmvB <- degree(outB[[k]]) == 0
        if(any(rmvB) & sum(rmvB) < vcount(outB[[k]])) outB[[k]] <- delete.vertices(outB[[k]], rmvB)
      }
    }
  }
  out <- list(Gyy = outTht, Gxy = outB, Gxx = outThtxx)
  class(out) <- "jcglasso2igraph"
  out
}

##### to_graph tools #####
getGraph2 <- function(x, type = c("Gyy", "Gxy", "Gxx", "conditional", "bipartite")){
  if (!is.jcglasso2igraph(x))
    stop(sQuote(deparse(substitute(x))), " is not an object of class", sQuote("jcglasso2igraph"))
  type <- match.arg(type)
  if (is.null(x$Gxy) & is.element(type, c("Gxy", "conditional", "bipartite")))
    stop(sQuote(type), " is not available. Please, use type = ", dQuote("Gyy"))
  if (type == "Gyy") gr <- x$Gyy
  if (type == "Gxy") gr <- x$Gxy
  if (type == "Gxx") gr <- x$Gxx
  if (type == "conditional") {
    K <- length(x$Gyy)
    gr <- vector(mode = "list", length = K)
    for(k in seq_len(K)) {
      gr.e <- rbind(as_data_frame(x$Gyy[[k]], what = "edges"), 
                    as_data_frame(x$Gxy[[k]], what = "edges"))
      v.yy <- as_data_frame(x$Gyy[[k]], what = "vertices")
      v.xy <- as_data_frame(x$Gxy[[k]], what = "vertices")
      id <- !is.element(v.xy$name, v.yy$name)
      gr.v <- rbind(v.yy, v.xy[id, ])
      gr[[k]] <- graph_from_data_frame(d = gr.e, directed = TRUE, vertices = gr.v)
    }
  }
  if (type == "bipartite") {
    K <- length(x$Gyy)
    gr <- vector(mode = "list", length = K)
    for(k in seq_len(K)) {
      gr.e <- rbind(as_data_frame(x$Gyy[[k]], what = "edges"), 
                    as_data_frame(x$Gxy[[k]], what = "edges"),
                    as_data_frame(x$Gxx[[k]], what = "edges"))
      v.yy <- as_data_frame(x$Gyy[[k]], what = "vertices")
      v.xy <- as_data_frame(x$Gxy[[k]], what = "vertices")
      v.xx <- as_data_frame(x$Gxx[[k]], what = "vertices")
      id <- !is.element(v.xy$name, v.yy$name)
      gr.v <- rbind(v.yy, v.xy[id, , drop = FALSE])
      id <- !is.element(v.xx$name, v.xy$name)
      gr.v <- rbind(gr.v, v.xx[id, , drop = FALSE])
      gr[[k]] <- graph_from_data_frame(d = gr.e, directed = TRUE, vertices = gr.v)
    }
  }
  gr
}

##### to_graph S3 methods #####
is.jcglasso2igraph <- function(x) inherits(x, 'jcglasso2igraph')

print.jcglasso2igraph <- function(x, ...) print.listof(x, ...)

plot.jcglasso2igraph <- function(x, type, which, highlight.connections = NULL, ...){
  if(missing(type)) type <- ifelse(is.null(x$Gxy), "Gyy", "conditional")
  gr <- getGraph2(x, type = type)
  opt <- list(...)
  K <- length(gr)
  if (missing(which)) which <- seq_len(K)
  else {
    if(!is.vector(which)) stop(sQuote("which"), " is not a vector")
    if(any(abs(as.integer(which) - which) > 0)) stop(sQuote("which"), " is not an object of type ", dQuote("integer"))
    if(min(which) <= 0) stop("some entry in ", sQuote("which"), " is not a positive integer")
    if(max(which) > K) stop("some entry in ", sQuote("which"), " is larger than ", sQuote(K))
  }
  if(!is.null(highlight.connections)){
    deg <- sapply(gr, degree, mode = "all", simplify = FALSE)
    if(is.logical(highlight.connections)) 
      if(is.logical(highlight.connections)) {
        if(!highlight.connections) stop(sQuote("highlight.connections"), " must be TRUE, or a value in (0, 1), or an integer > 1 or a character vector")
        prbs <- .95
        high <- "0" 
      }
    if(is.numeric(highlight.connections)) {
      if(highlight.connections > 0 & highlight.connections < 1){
        prbs <- highlight.connections
        high <- "1"
      } 
      else {
        prbs <- 0
        high <- "2"
      }
    }
    if(is.character(highlight.connections)) {
      prbs <- 0
      high <- "3"
    }
    connections <- quantile(unlist(deg), probs = prbs)
    id <- sapply(deg, function(x, high) {
      temp <- switch(high, 
                     "0" = which(x >= connections), 
                     "1" = which(x >= connections),
                     "2" = match(names(tail(sort(x), highlight.connections)), names(x)), 
                     "3" = match(highlight.connections, names(x)))
      if(all(is.na(temp))) stop(sQuote("highlight.connections"), " must be TRUE, or a value in (0, 1), or an integer > 1 or a character vector")
      names(temp) <- names(x)[temp]
      temp
    }, high = high, simplify = FALSE)
    id2 <- sort(unique(unlist(sapply(id, function(x) names(x), simplify = FALSE))))
    if(length(id2) != 0) {
      col <- rainbow(length(id2))
      names(col) <- id2
    }
    
    for(j in which) {
      # Set colors to plot the selected edges.
      ecol <- rep(gray(.8, .5), ecount(gr[[j]]))
      vcol <- rep(gray(.8, .5), vcount(gr[[j]]))
      
      if(length(id2) != 0 & length(id[[j]]) != 0) {
        for(i in 1:length(id[[j]])){
          inc.edges <- incident(gr[[j]],  V(gr[[j]])[id[[j]][i]], mode="all")
          ecol[inc.edges] <- col[names(id[[j]][i])]
          vcol[id[[j]][i]] <- col[names(id[[j]][i])]
        }
      }
      plot(gr[[j]], vertex.color=vcol, edge.color=ecol)
    }
  } 
  else {
    #    if(is.null(opt$layout)) opt$layout <- layout_with_kk
    for(k in which){ do.call(function(...) plot(gr[[k]], ...), opt) }  
  }
  invisible(NULL)
}


