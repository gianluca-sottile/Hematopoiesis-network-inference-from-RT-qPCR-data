to_graph <- function(object, GoF = AIC, lambda.id, rho.id, weighted = FALSE, simplify = TRUE, ...) {
    if (!inherits(object, c("cggm", "cglasso"))) stop(sQuote("to_graph"), "function is not available for an object of class ", sQuote(class(object)))
    q <- npred(object$Z)
    nlambda <- object$nlambda
    nrho <- object$nrho
    if (missing(lambda.id) & missing(rho.id)) {
        if (!is.element(class(GoF), c("function", "GoF")))
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
        object <- select.cglasso(object, GoF = GoF)
        nlambda <- object$nlambda
        lambda.id <- 1L
        nrho <- object$nrho
        rho.id <- 1L
    } else {
        # testing 'lambda.id'
        nlambda <- object$nlambda
        if (missing(lambda.id)) {
            if (q == 0 | nlambda == 1L) lambda.id <- 1L
            else stop(sQuote("lambda.id"), " is missing")
        } else {
            if (!is.vector(lambda.id)) stop(sQuote("lambda.id"), " is not a vector")
            if (length(lambda.id) != 1L) stop(sQuote("lambda.id"), " is not a vector of length ", sQuote("1"))
            if (any(abs(as.integer(lambda.id) - lambda.id) > 0)) stop(sQuote("lambda.id"), " is not an object of type ", dQuote("integer"))
            if (lambda.id <= 0) stop(sQuote("lambda.id"), " is not a positive integer")
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
    }
    if (!is.logical(simplify)) stop(sQuote("simplify"), " is not an object of type ", sQuote("logical"))
    if (!is.logical(weighted)) stop(sQuote("weighted"), " is not an object of type ", dQuote("logical"))
    Adj_yy <- object$InfoStructure$Adj_yy[, , lambda.id, rho.id]
    pY <- nrow(Adj_yy)
    nmsY <- rownames(Adj_yy)
    outTht <- graph_from_adjacency_matrix(adjmatrix = Adj_yy, mode = "undirected", diag = FALSE)
    # setting vertex attribute
    V(outTht)$type <- rep("response", pY)
    V(outTht)$color <- NA
    V(outTht)$frame.color <- NA
    V(outTht)$size <- 12
    V(outTht)$label.cex <- 1
    V(outTht)$label.font <- 2
    V(outTht)$label.dist <- 0
    E(outTht)$type <- "undirect"
    E(outTht)$color <- "gray50"
    E(outTht)$width <- 1
    E(outTht)$arrow.mode <- 0
    outB <- NULL
    if (q != 0L){
        Adj_xy <- object$InfoStructure$Adj_xy[, , lambda.id, rho.id, drop = FALSE]
        V(outTht)$label.color <- ifelse(apply(Adj_xy, 2L, any), "black", "gray50")
        E(outTht)$color <- "gray85"
        pX <- nrow(Adj_xy)
        nmsX <- rownames(Adj_xy)
        eB <- data.frame(from = rep(nmsX, pY), to = rep(nmsY, each = pX))[c(Adj_xy) != 0,]
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
        E(outB)$width <- 1
        E(outB)$arrow.mode <- 2
    }
    if (weighted) {
        Tht <- object$Tht[, , lambda.id, rho.id]
        E(outTht)$weight <- Tht[lower.tri(Tht)][Tht[lower.tri(Tht)] != 0]
        E(outTht)$lty <- ifelse(E(outTht)$weight > 0, "solid", "dashed")
        if(!is.null(object$InfoStructure$Adj_xy)) {
            B <- object$B[-1, , lambda.id, rho.id]
            E(outB)$weight <- B[B != 0]
            E(outB)$lty <- ifelse(E(outB)$weight > 0, "solid", "dashed")
        }
    }
    if (simplify) {
        rmvTht <- degree(outTht) == 0
        if (any(rmvTht)& sum(rmvTht) < vcount(outTht)) outTht <- delete.vertices(outTht, rmvTht)
        if (q != 0L) {
            rmvB <- degree(outB) == 0
            if(any(rmvB) & sum(rmvB) < vcount(outB)) outB <- delete.vertices(outB, rmvB)
        }
    }
    out <- list(Gyy = outTht, Gxy = outB)
    class(out) <- "cglasso2igraph"
    out
}
