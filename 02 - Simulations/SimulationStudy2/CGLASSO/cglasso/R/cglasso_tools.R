impute <- function(object, type = c("mar", "censored", "both"), lambda.new, rho.new) {
    type <- match.arg(type)
    if (!inherits(object, c("cggm", "cglasso"))) stop(sQuote("imput"), "function is not available for an object of class ", sQuote(class(object)))
    if ( class(object)[1L] == "cggm") {
        if (!missing(lambda.new)) stop("Argument ", sQuote("lambda.new"), " is not available for an object of class ", sQuote("cggm"))
        if (!missing(rho.new)) stop("Argument ", sQuote("rho.new"), " is not available for an object of class ", sQuote("cggm"))
        p <- npred(object$Z)
        out <- drop(object$Yipt)
        Id <- event(object$Z)
        if (type == "mar") {
            lo <- lower(object$Z)
            up <- upper(object$Z)
            for (h in seq_len(p)) {
                out[Id[, h] == -1L, h] <- lo[h]
                out[Id[, h] == +1L, h] <- up[h]
            }
        }
        if (type == "censored")
            out[Id == +9L] <- NA
        return(out)
    }
    if (missing(lambda.new) & missing(rho.new)) return(object$Yipt)
    nrho <- object$nrho
    rho <- object$rho
    nlambda <- object$nlambda
    lambda <- object$lambda
    n <- nobs(object$Z)
    p <- nresp(object$Z)
    q <- npred(object$Z)
    if (missing(lambda.new)) {
        if (q > 0) stop(sQuote("lambda.new"), " is missing")
        else lambda.new <- lambda
    }
    else {
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
    Yin <- object$Yipt
    Yout <- matrix(0, nrow = n, ncol = p, dimnames = list(dimnames(Yin)[[1L]], dimnames(Yin)[[2L]]))
    Id <- event(object$Z)
    ########################
    # setting storage.mode #
    ########################
    storage.mode(rho.new) <- "double"
    storage.mode(lambda.new) <- "double"
    storage.mode(nrho) <- "integer"
    storage.mode(rho) <- "double"
    storage.mode(nlambda) <- "integer"
    storage.mode(lambda) <- "double"
    storage.mode(n) <- "integer"
    storage.mode(p) <- "integer"
    storage.mode(Yin) <- "double"
    storage.mode(Id) <- "integer"
    storage.mode(Yout) <- "double"
    out <- .Fortran(C_impute, newrho = rho.new, newlambda = lambda.new, nrho = nrho,
                      rho = rho, nlambda = nlambda, lambda = lambda, n = n, p = p,
                      Yin = Yin, Id = Id, Yout = Yout)
    if (type == "mar") {
        lo <- lower(object$Z)
        up <- upper(object$Z)
        for (h in seq_len(p)) {
            out$Yout[Id[, h] == -1L, h] <- lo[h]
            out$Yout[Id[, h] == +1L, h] <- up[h]
        }
    }
    if (type == "censored")
        out$Yout[Id == +9L] <- NA
    out$Yout
}

function2xy <- function(frml) {
    x <- match.arg(all.vars(frml)[2L], c("rho", "lambda"))
    y <- as.character(terms(frml)[[2L]])
    if (length(y) == 2L) {
        if(y[1L] != "diag" | y[2L] != "Theta")
            stop(sQuote(paste0(y[1L], "(", y[2L], ")")), " is not available as response. Please, use ", sQuote("diag(Theta)"))
        y <- "diag(Theta)"
    }
    given <- as.character(terms(frml)[[3L]])
    if (length(given) == 1L) given <- NULL
    else {
        if (given[3L] == "rho") {
            if (x == "rho") stop("You can not condition on ", sQuote(given[3L]), ". Please, use ", sQuote("lambda"))
            return(list(x = x, y = y, given = NULL))
        }
        if (given[3L] == "lambda") {
            if (x == "lambda") stop("You can not condition on ", sQuote(given[3L]), ". Please, use ", sQuote("rho"))
            return(list(x = x, y = y, given = NULL))
        }
        given.terms <- gsub("[[:blank:]]", "", given[3L])
        given.terms <- unlist(strsplit(given.terms, split = "\n", fixed = TRUE))
        if (length(given.terms) == 1L) given <- eval(parse(text = given.terms))
        else {
            given.terms <- given.terms[2L]
            given.name <- unlist(strsplit(given.terms, "=", fixed = TRUE))[1L]
            if (!is.element(given.name, c("lambda.id", "rho.id")))
                stop("Please, use ", sQuote(ifelse(x == "rho", "lambda.id", "rho.id")), " instead of ", sQuote(given.name))
                #stop("Please, use ", sQuote("lambda.id"), " or ", sQuote("rho.id"), " instead of ", sQuote(given.name))
            if ((x == "rho" &  given.name == "rho.id") | (x == "lambda" &  given.name == "lambda.id"))
                stop("You can not condition on ", sQuote(given.name), ". Please, use ", sQuote(ifelse(given.name == "rho.id", "lambda.id", "rho.id")))
            given <- eval(parse(text = given.terms))
        }
    }
    out <- list(penalty = x, what = y, given = given)
    out
}
