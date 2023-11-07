## Wrapper functions to make the new class compatible with the 
## datacggm class of the R package cglasso

## Create a Dataset from a Conditional Gaussian Graphical Model 
## with Censored and/or Missing Values
## (Much of the documentation for these functions can be found in 
##  the manuals for the corresponding functions of the datacggm class.)

##### datajcggm class #####
datajcggm <- function (Y, lo = -Inf, up = +Inf, X = NULL, 
                       control = list(maxit = 10000, thr = 0.0001)) {
  big <- .Machine$double.xmax
  if (missing(Y)) stop(sQuote("Y"), " is missing")
  Yorig <- Y
  if(!is.list(Y)) Yorig <- list(Y)
  K <- length(Y)
  if(!is.null(X)) {
    Xorig <- X
    if(!is.list(X)) Xorig <- list(X)
    if(length(X) != K) stop(sQuote("X"), " and ", sQuote("Y"), " have a different number of layers")
  }
  Zout <- vector(mode = "list", length = K)
  for(k in seq_len(K)) {
    thr <- big / 2
    Y <- Yorig[[k]]
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
      X <- Xorig[[k]]
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
    loz <- if(q != 0) lox else loy
    upz <- if(q != 0) upx else upy
    if (!all(loz < upz)) stop(sQuote("lo"), " is not less than ", sQuote("up"))
    for (m in seq_len(p + q)) {
      lo.id <- which(Z[, m] < loz[m])
      up.id <- which(Z[, m] > upz[m])
      if (length(lo.id) != 0L) {
        message("message: some entries in column ", sQuote(znames[m]), 
                " are below the lower censoring value. These entries are replaced with ", 
                sQuote(loz[m]))
        Z[lo.id, m] <- loz[m]
      }
      if (length(up.id) != 0L) {
        message("message: some entries in column ", sQuote(znames[m]), 
                " are over the upper censoring value. These entries are replaced with ", 
                sQuote(upz[m]))
        Z[up.id, m] <- upz[m]
      }
    }
    Zna <- -9999
    Z[is.na(Z)] <- Zna
    storage.mode(n) <- "integer"
    storage.mode(p) <- "integer"
    storage.mode(Y) <- "double"
    storage.mode(loz) <- "double"
    storage.mode(upz) <- "double"
    storage.mode(Zna) <- "double"
    Info <- .Fortran(cglasso:::C_setup, n = n, p = p + q, Y = Z, lo = loz, up = upz, 
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
    Zout[[k]] <- Z
  }
  
  class(Zout) <- "datajcggm"
  Zout
}

##### datajcggm S3 methods #####
is.datajcggm <- function (x) inherits(x, "datajcggm")

dim.datajcggm <- function (x) {
  K <- 
  out <- list(Y = lapply(getMatrix2(x, "Y"), dim), X = lapply(getMatrix2(x, "X"), dim))
  out
}

nresp.datajcggm <- function (object, ...) unique(sapply(object, \(x) x$Info$p))

nobs.datajcggm <- function (object, ...) sapply(object, \(x) x$Info$n)

npred.datajcggm <- function (object, ...) unique(sapply(object, \(x) x$Info$q))

dimnames.datajcggm <- function (x) {
  Y <- getMatrix2(x, "Y")
  if (npred(x) != 0) X <- getMatrix2(x, "X")
  K <- length(x)
  lapply(seq_len(K), \(k) {
    out <- list(Y = dimnames(Y[[k]]), X = NULL)
    if (npred(x) != 0) out$X <- dimnames(X[[k]])
    out  
  })
}

as.character.datajcggm <- function (x, digits = 3L, ...) {
  p <- nresp(x); q <- npred(x); K <- length(x)
  Z <- lapply(getMatrix2(x, "both"), \(x) {
    format(cbind(x$Y, x$X), digits = digits, ...)
  })
  marker <- event2(x)
  lapply(seq_len(K), \(k) {
    Z[[k]][marker[[k]] == 1] <- paste0(Z[[k]][marker[[k]] == 1], "+")
    Z[[k]][marker[[k]] == -1] <- paste0(Z[[k]][marker[[k]] == -1], "-")
    if(q != 0) 
      list(Y = Z[[k]][, 1:p, drop = FALSE], X = Z[[k]][, (p+1):(p+q), drop = FALSE])
    else list(Y = Z[[k]][, 1:p, drop = FALSE], X = NULL)  
  })
}

print.datajcggm <- function (x, which, digits = 3L, n = 10L, width = getOption("width"), ...) {
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
  p <- nresp(x); q <- npred(x); K <- length(x)
  if (missing(which)) which <- seq_len(K)
  if (q != 0L) 
  cat("Printing", sQuote("datajcggm"), "object\n")
  for(k in which) {
    cat("\nY:", nobs(x), "x", p, "matrix\n\n")
    out[[k]]$Y <- fit_matrix(out[[k]]$Y, n = n, width = width)
    print(out[[k]]$Y, quote = FALSE)
    dn <- nobs(x)[k] - n
    dp <- p - dim(out[[k]]$Y)[2L]
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
      out[[k]]$X <- fit_matrix(out[[k]]$X, n = n, width = width)
      cat("\n\nX:", nobs(x), "x", npred(x), "matrix\n\n")
      print(out[[k]]$X, quote = FALSE)
      dq <- npred(x) - dim(out[[k]]$X)[2L]
      if (dn > 0 | dq > 0) {
        msg <- paste0("# with ", ifelse(dn > 0, paste0(dn, " more ", ifelse(dn > 1, "rows", "row")), ""), 
                      ifelse(dn > 0 & dq > 0, ", and ", ""), 
                      ifelse(dq > 0, paste0(dq, " more ", ifelse(dq > 1, "variables", "variable")), ""))
        cat(msg)
      }
    }
    cat("\n")
  }
}

summary.datajcggm <- function (object, which, n, quantile.type = 7L, digits = 3L, quote = FALSE, ...) {
  if (missing(n)) n <- Inf
  if (missing(which)) which <- seq_len(length(object))
  Y <- getMatrix2(object, "Y")
  X <- getMatrix2(object, "X")
  R <- event2(object)
  p <- nresp(object)
  q <- npred(object)
  lo <- lapply(lower2(object), \(x) c(x$Y, x$X))
  up <- lapply(upper2(object), \(x) c(x$Y, x$X))
  for(k in which) {
    out <- list(Y = NULL, X = NULL)
    lcs <- paste0(format(apply(R[[k]] == -1, 2L, function(x) mean(x) * 100), digits = digits), "%")
    rcs <- paste0(format(apply(R[[k]] == +1, 2L, function(x) mean(x) * 100), digits = digits), "%")
    nas <- paste0(format(apply(R[[k]] == +9, 2L, function(x) mean(x) * 100), digits = digits), "%")
    ynm <- colNames2(object)[[k]]$Y
    qq <- matrix(0, nrow = p, ncol = 6L)
    for (m in seq_len(p)) {
      obs <- Y[[k]][R[[k]][, m] == 0, m]
      qq[m, -4L] <- stats::quantile(obs, names = FALSE, type = quantile.type)
      qq[m, 4L] <- mean(obs)
    }
    tbl <- apply(cbind(lo[[k]][1:p], qq, up[[k]][1:p]), 2L, format, digits = digits)
    tbl <- cbind(tbl, nas[1:p], lcs[1:p], rcs[1:p])
    colnames(tbl) <- c("Lower", "Min.", "1st Qu.", "Median", 
                       "Mean", "3rd Qu.", "Max.", "Upper", 
                       "NA%", "LC%", "RC%")
    rownames(tbl) <- ynm
    attr(tbl, "class") <- "table"
    out$Y <- tbl
    if (q > 0L) {
      xnm <- colNames2(object)[[k]]$X
      qq <- matrix(0, nrow = q, ncol = 6L)
      for (m in seq_len(q)) {
        obs <- X[[k]][R[[k]][, p + m] == 0, m]
        qq[m, -4L] <- stats::quantile(obs, names = FALSE, type = quantile.type)
        qq[m, 4L] <- mean(obs)
      }
      tbl <- apply(cbind(lo[[k]][(p+1):(p+q)], qq, up[[k]][(p+1):(p+q)]), 2L, format, digits = digits)
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
  }
  
  invisible(out)
}

##### datajcggm auxiliary function methods #####

getMatrix2 <- function (x, name = c("Y", "X", "both"), ordered = FALSE) {
  if (!is.datajcggm(x)) 
    stop(sQuote(x), " is not an object of class ", sQuote("datajcggm"))
  if (!is.logical(ordered)) 
    stop(sQuote(ordered), " id not an object of type ", sQuote("logical"))
  name <- match.arg(name)
  K <- length(x)
  lapply(seq_len(K), \(k) {
    if (ordered) {
      row.order <- x[[k]]$Info$order
      Y <- x[[k]]$Y[row.order, , drop = FALSE]
      X <- x[[k]]$X[row.order, , drop = FALSE]
    }
    else {
      Y <- x[[k]]$Y
      X <- x[[k]]$X
    }
    out <- switch(name, Y = Y, X = X, both = list(Y = Y, X = X))
    out
  })
}

rowNames2 <- function (x) {
  if (!is.datajcggm(x)) 
    stop(sQuote(x), " is not an object of class ", sQuote("datajcggm"))
  K <- length(x)
  Y <- getMatrix2(x, "Y")
  X <- getMatrix2(x, "X")
  lapply(seq_len(K), \(k) {
    if (is.null(x[[k]]$X)) 
      out <- list(Y = rownames(Y[[k]]), X = NULL)
    else out <- list(Y = rownames(Y[[k]]), X = rownames(X[[k]]))
    out
  })
}

colNames2 <- function (x) {
  if (!is.datajcggm(x)) 
    stop(sQuote(x), " is not an object of class ", sQuote("datajcggm"))
  K <- length(x)
  Y <- getMatrix2(x, "Y")
  X <- getMatrix2(x, "X")
  lapply(seq_len(K), \(k) {
    if (is.null(x[[k]]$X)) 
      out <- list(Y = colnames(Y[[k]]), X = NULL)
    else out <- list(Y = colnames(Y[[k]]), X = colnames(X[[k]]))
    out
  })
}

event2 <- function (x, ordered = FALSE) {
  if (!is.datajcggm(x)) 
    stop("'x' is not an object of class 'datajcggm'")
  if (!is.logical(ordered)) 
    stop(sQuote(ordered), " id not an object of type ", sQuote("logical"))
  K <- length(x)
  lapply(seq_len(K), \(k) {
    R <- x[[k]]$Info$R
    if (ordered) 
      R <- R[x[[k]]$Info$order, , drop = FALSE]
    R
  })
}

lower2 <- function (x) {
  if (!is.datajcggm(x)) 
    stop(sQuote(x), " is not an object of class ", sQuote("datajcggm"))
  small <- -.Machine$double.xmax
  p <- nresp(x); q <- npred(x); K <- length(x)
  lapply(seq_len(K), \(k) {
    out <- x[[k]]$Info$lo
    out[out == small] <- -Inf
    if (q != 0L) {
      out <- list(Y = out[1:p], X = out[(p+1):(p+q)])
    } else out <- list(Y = out[1:p], X = NULL)
    out
  })
}

upper2 <- function (x) {
  if (!is.datajcggm(x)) 
    stop(sQuote(x), " is not an object of class ", sQuote("datajcggm"))
  big <- .Machine$double.xmax
  p <- nresp(x); q <- npred(x); K <- length(x)
  lapply(seq_len(K), \(k) {
    out <- x[[k]]$Info$up
    out[out == big] <- +Inf
    if (q != 0L) {
      out <- list(Y = out[1:p], X = out[(p+1):(p+q)])
    } else out <- list(Y = out[1:p], X = NULL)
    out
  })
}

ColMeans2 <- function (x) {
  if (!is.datajcggm(x)) 
    stop(sQuote("x"), " is not a ", sQuote("datajcggm"), " object")
  p <- nresp(x); q <- npred(x); K <- length(x)
  lapply(seq_len(K), \(k) {
    if (q != 0L) {
      out <- list(Y = x[[k]]$Info$zm[1:p], X = x[[k]]$Info$zm[(p+1):(p+q)])
    } else out <- list(Y = x[[k]]$Info$zm, X = NULL)
    out
  })
}

ColVars2 <- function (x) {
  if (!is.datajcggm(x)) 
    stop(sQuote("x"), " is not a ", sQuote("datajcggm"), " object")
  p <- nresp(x); q <- npred(x); K <- length(x)
  lapply(seq_len(K), \(k) {
    if (q != 0L) {
      out <- list(Y = x[[k]]$Info$zv[1:p], X = x[[k]]$Info$zv[(p+1):(p+q)])
    } else out <- list(Y = x[[k]]$Info$zv, X = NULL)
    out
  })
}
