##### GoF functions, used for model selection #####
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
                n = out_QFun$n, p = out_QFun$p, q = out_QFun$q,
                rho = out_QFun$rho, nrho = out_QFun$nrho, lambda = out_QFun$lambda, nlambda = out_QFun$nlambda, 
                nu = out_QFun$nu, alpha1 = out_QFun$alpha1, alpha2 = out_QFun$alpha2, alpha3 = out_QFun$alpha3,
                type = type, model = out_QFun$model, components = components)
  }
  else {
    value <-  out_QFun$value
    df <- sapply(seq_len(K), function(h) k[h] * out_QFun$df.comp[[h]])
    if(is.vector(df)) df <- t(df)
    df <- matrix(rowSums(df), nrow = nlambda, ncol = nrho, dimnames = dimnames(value))
    val <- -2 * value + df
    out <- list(value_gof = val, df = out_QFun$df, dfB = out_QFun$dfB, dfTht = out_QFun$dfTht, value = value,
                n = out_QFun$n, p = out_QFun$p, q = out_QFun$q,
                rho = out_QFun$rho, nrho = out_QFun$nrho, lambda = out_QFun$lambda, nlambda = out_QFun$nlambda, 
                alpha1 = out_QFun$alpha1, alpha2 = out_QFun$alpha2, alpha3 = out_QFun$alpha3, 
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
                  n = out_QFun$n, p = out_QFun$p, q = out_QFun$q, 
                  rho = out_QFun$rho, nrho = out_QFun$nrho, lambda = out_QFun$lambda, nlambda = out_QFun$nlambda, 
                  nu = out_QFun$nu, alpha1 = out_QFun$alpha1, alpha2 = out_QFun$alpha2, alpha3 = out_QFun$alpha3,
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
                  n = out_QFun$n, p = out_QFun$p, q = out_QFun$q, 
                  rho = out_QFun$rho, nrho = out_QFun$nrho, lambda = out_QFun$lambda, nlambda = out_QFun$nlambda, 
                  nu = out_QFun$nu, alpha1 = out_QFun$alpha1, alpha2 = out_QFun$alpha2, alpha3 = out_QFun$alpha3,
                  type = type, model = out_QFun$model, components = components)
    }
  }
  class(out) <- "GoF2"
  out
}

##### GoF S3 methods #####
print.GoF2 <- function (x, digits = 3L, ...){
  type <- x$type
  nrho <- x$nrho
  rho <- x$rho
  nlambda <- x$nlambda
  lambda <- x$lambda
  nu <- x$nu
  alpha1 <- x$alpha1
  alpha2 <- x$alpha2
  alpha3 <- x$alpha3
  q <- x$q
  
  if (nrho > 1L & nlambda > 1L) {
    lambda <- rep(lambda, each = nrho)
    rho <- rep(rho, times = nlambda)
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
    tbl <- data.frame(rho = rho, lambda = lambda, nu = if(is.null(nu)) rho else nu, 
                      alpha1 = alpha1, alpha2 = alpha2, alpha3 = alpha3, 
                      df = df, df.Tot = rowSums(df), 
                      val = val, val.Tot = rowSums(val))
    names(tbl)[7L:(7L + K - 1L)] <- paste0("df.Comp", seq_len(K))
    names(tbl)[-c(1L:(7L + K), ncol(tbl))] <- paste0(type, ".Comp", seq_len(K))
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
    tbl <- data.frame(rho = rho, lambda = lambda, nu = if(is.null(nu)) rho else nu, 
                      alpha1 = alpha1, alpha2 = alpha2, alpha3 = alpha3, 
                      df = df, val = val)
    names(tbl)[8L] <- type
  }
  
  cat("\nSequence of", sQuote(type), "values of the fitted", sQuote(x$model[1L]), ifelse(nrho == 1L & nlambda == 1L, "model", "models"))
  cat("\n\nDetails:\n")
  dots <- list(...)
  if (is.null(dots$print.gap)) dots$print.gap <- 2L
  if (is.null(dots$quote)) dots$quote <- FALSE
  if (is.null(dots$row.names)) dots$row.names <- FALSE
  if (q == 0L) tbl <- tbl[, -c(2L, 3L, 5L, 6L), drop = FALSE]
  
  if (nlambda <= nrho) f <- rep(seq_len(nlambda), each = nrho)
  else f <- rep(seq_len(nrho), times = nlambda)
  tbl.list <- split(tbl, f = f, drop = FALSE)
  do.call(function(...) print.listof(tbl.list, digits = digits, ...), dots)
  
  invisible(tbl)
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
  rho <- x$rho
  nlambda <- x$nlambda
  lambda <- x$lambda
  alpha1 <- x$alpha1
  alpha2 <- x$alpha2
  alpha3 <- x$alpha3
  
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

##### GoF auxiliary functions #####
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
  object$Zipt <- object$Zipt[ , , , lambda.id, rho.id, drop = FALSE]
  object$B <- object$B[ , , , lambda.id, rho.id, drop = FALSE]
  object$mu <- object$mu[ , , , lambda.id, rho.id, drop = FALSE]
  object$R <- object$R[ , , , lambda.id, rho.id, drop = FALSE]
  object$S <- object$S[ , , , lambda.id, rho.id, drop = FALSE]
  object$Sgm <- object$Sgm[ , , , lambda.id, rho.id, drop = FALSE]
  object$Tht <- object$Tht[ , , , lambda.id, rho.id, drop = FALSE]
  if(q > 0) object$Omega <- object$Omega[ , , , lambda.id, rho.id, drop = FALSE]
  object$dfB <- object$dfB[ , , lambda.id, rho.id, drop = FALSE]
  object$dfTht <- object$dfTht[, lambda.id, rho.id, drop = FALSE]
  object$dfOmg <- object$dfOmg[, lambda.id, rho.id, drop = FALSE]
  object$InfoStructure$Adj <- object$InfoStructure$Adj[ , , , lambda.id, rho.id, drop = FALSE]
  object$InfoStructure$ncomp <- object$InfoStructure$ncomp[lambda.id, rho.id, drop = FALSE]
  object$InfoStructure$Ck <- object$InfoStructure$Ck[ , lambda.id, rho.id, drop = FALSE]
  object$InfoStructure$pk <- object$InfoStructure$pk[ , lambda.id, rho.id, drop = FALSE]
  object$nit <- object$nit[ , lambda.id, rho.id, drop = FALSE]
  object$connected <- object$connected#[ , lambda.id, rho.id, drop = FALSE]
  
  object$nlambda <- 1L
  object$lambda.min.ratio <- 1
  object$lambda <- object$lambda[lambda.id]
  object$nrho <- 1L
  object$rho.min.ratio <- 1
  object$rho <- object$rho[rho.id]
  
  object
}

## QFun function, to extracts the values of the Q-function from an R object inheriting class ‘jcglasso’ ##
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
  rho <- object$rho
  nrho <- object$nrho
  lambda <- object$lambda
  nlambda <- object$nlambda
  alpha1 <- object$alpha1
  alpha2 <- object$alpha2
  alpha3 <- object$alpha3
  X.null <- (q == 0L)
  
  lambda_lab <- paste0("lambda_", seq_len(nlambda))
  rho_lab <- paste0("rho_", seq_len(nrho))
  value <- df <- dfB <- dfTht <- dfOmg <- array(0.0, dim = c(nlambda, nrho, K), 
                                                dimnames = list(lambda_lab, rho_lab, NULL))
  value.tot <- df.tot <- dfB.tot <- dfTht.tot <- dfOmg.tot <- matrix(0.0, nrow = nlambda, ncol = nrho, 
                                                                     dimnames = list(lambda_lab, rho_lab))
  value.comp <- df.comp <- dfB.comp <- dfTht.comp <- vector(mode = "list", length = K)
  const1 <- (- 0.5 * n * p * log(2.0 * pi))
  const2 <- 0.5 * n
  ntp <- nlambda * nrho
  ii <- cbind(i = rep(seq_len(nlambda), each = nrho), j = rep(seq_len(nrho), times = nlambda))
  if(verbose) {
    cat("\nComputing Q-values of the fitted", sQuote(object$model[1L]), 
        ifelse(nrho == 1L & nlambda == 1L, "model", "models"), "\n")
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
          test1[k] <- all(Adj[id_Y, id_Y, k, i, j] == Adj[id_Y, id_Y, k, i1, j1])
          if (q == 0L) test2[k] <- TRUE
          else test2[k] <- all(Adj[id_X, id_Y, k, i, j] == Adj[id_X, id_Y, k, i1, j1])
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
          S <- out.mle$S[id_Y, id_Y, k, 1L, 1L]
          Tht <- out.mle$Tht[id_Y, id_Y, k, 1L, 1L]
          
          dfB[i, j, k] <- out.mle$dfB[p + 1L, k, 1L, 1L]
          dfTht[i, j, k] <- out.mle$dfTht[k, 1L, 1L]
          if(!X.null) dfOmg[i, j, k] <- out.mle$dfOmg[k, 1L, 1L]
          df[i, j, k] <- dfB[i, j, k] + p + dfTht[i, j, k] + p
          # if(!X.null) df[i, j, k] <- df[i, j, k] + dfOmg[i, j, k] + q
          value[i, j, k] <- const1[k] + const2[k] * (determinant(Tht)$modulus - sum(S * Tht))
        }
      } 
      else {
        for(k in seq_len(K)) {
          value[i, j, k] <- value[i1, j1, k]
          dfB[i, j, k] <- dfB[i1, j1, k]
          dfTht[i, j, k] <- dfTht[i1, j1, k]
          if(!X.null) dfOmg[i, j, k] <- dfOmg[i1, j1, k]
          df[i, j, k] <- df[i1, j1, k]
        }
      }
    } 
    else {
      for(k in seq_len(K)) {
        S <- object$S[id_Y, id_Y, k, i, j]
        Tht <- object$Tht[id_Y, id_Y, k, i, j]
        
        dfB[i, j, k] <- object$dfB[p + 1L, k, i, j]
        dfTht[i, j, k] <- object$dfTht[k, i, j]
        if(!X.null) dfOmg[i, j, k] <- object$dfOmg[k, i, j]
        df[i, j, k] <- dfB[i, j, k] + p + dfTht[i, j, k] + p
        # if(!X.null) df[i, j, k] <- df[i, j, k] + dfOmg[i, j, k] + q
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
              n = n, p = p, q = q, rho = rho, nrho = nrho, lambda = lambda, nlambda = nlambda, 
              nu = object$nu, alpha1 = object$alpha1, alpha2 = object$alpha2, alpha3 = object$alpha3, model = object$model)
  class(out) <- "QFun2"
  out
}

## QFun S3 method ##
print.QFun2 <- function (x, digits = 3L, ...){
  nrho <- x$nrho
  rho <- x$rho
  nlambda <- x$nlambda
  lambda <- x$lambda
  nu <- x$nu
  alpha1 <- x$alpha1
  alpha2 <- x$alpha2
  alpha3 <- x$alpha3
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
    lambda <- rep(lambda, each = nrho)
    rho <- rep(rho, times = nlambda)
  }
  
  dots <- list(...)
  if (is.null(dots$print.gap)) dots$print.gap <- 2L
  if (is.null(dots$quote)) dots$quote <- FALSE
  if (is.null(dots$row.names)) dots$row.names <- FALSE
  tbl <- data.frame(rho = rho, lambda = lambda, nu = if(is.null(nu)) rho else nu, 
                    alpha1 = alpha1, alpha2 = alpha2, alpha3 = alpha3, 
                    df = df, value = value)
  names(tbl)[8L] <- "Q-Values"
  if (q == 0L) tbl <- tbl[, -c(2L, 3L, 5L, 6L), drop = FALSE]
  
  cat("\nQ-values of the fitted", sQuote(x$model[1L]), 
      ifelse((nrho == 1L & nlambda == 1L), "model", "models"))
  cat("\n\nDetails:\n")
  if (nlambda <= nrho) f <- rep(seq_len(nlambda), each = nrho)
  else f <- rep(seq_len(nrho), times = nlambda)
  
  tbl.list <- split(tbl, f = f, drop = FALSE)
  do.call(function(...) print.listof(tbl.list, digits = digits, ...), dots)
  
  invisible(tbl)
}