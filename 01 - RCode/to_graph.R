##### to_graph functions, to create Graphs from jcglasso or jcggm Objects #####
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
  nmsY <- colNames2(object$Z)[[1]]$Y
  outTht <- lapply(1:K, function(k, pY) {
    # dimnames(x) <- list(nmsY, nmsY, NULL, NULL)
    out <- graph_from_adjacency_matrix(adjmatrix = Adj[id_Y, id_Y, k, lambda.id, rho.id], 
                                       mode = "undirected", diag = FALSE)
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
    outThtxx <- lapply(1:K, function(k, qX) {
      out <- graph_from_adjacency_matrix(adjmatrix = Adj[id_X, id_X, k, lambda.id, rho.id], 
                                         mode = "undirected", diag = FALSE)
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
    nmsX <- colNames2(object$Z)[[1]]$X
    for(k in seq_len(K)) {
      V(outTht[[k]])$label.color <- ifelse(apply(Adj[id_X, id_Y, k, lambda.id, rho.id, drop = FALSE], 
                                                 2L, function(x) any(x == 1)), "black", "gray50")
      E(outTht[[k]])$color <- "gray85"
      V(outThtxx[[k]])$label.color <- ifelse(apply(Adj[id_X, id_Y, k, lambda.id, rho.id, drop = FALSE], 
                                                   1L, function(x) any(x == 1)), "darkblue", "gray50")
      E(outThtxx[[k]])$color <- "gray85"
    }
    outB <- sapply(1:K, function(k, pX, pY, nmsX, nmsY, lambda.id, rho.id) {
      out <- Adj[id_X, id_Y, k, lambda.id, rho.id, drop = FALSE]
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
      Tht <- -cov2cor(object$Tht[id_Y, id_Y, k, lambda.id, rho.id])
      weight <- Tht[lower.tri(Tht)][Tht[lower.tri(Tht)] != 0]
      E(outTht[[k]])$lty <- ifelse(weight > 0, "solid", "dashed")
      E(outTht[[k]])$width <- as.numeric(as.character(cut(abs(weight), 
                                                          breaks = c(0, .0001, .01, 1), 
                                                          labels = c(.5, 1.5, 2.5))))
      E(outTht[[k]])$weight <- weight
      if(q > 0L) {
        Tht <- -cov2cor(object$Omega[, , k, lambda.id, rho.id])
        weight <- Tht[lower.tri(Tht)][Tht[lower.tri(Tht)] != 0]
        E(outThtxx[[k]])$lty <- ifelse(weight > 0, "solid", "dashed")
        E(outThtxx[[k]])$width <- as.numeric(as.character(cut(abs(weight), 
                                                              breaks = c(0, .0001, .01, 1), 
                                                              labels = c(.5, 1.5, 2.5))))
        E(outThtxx[[k]])$weight <- weight
        
        B <- object$B[-1L, , k, lambda.id, rho.id]
        weight <- B[B != 0]
        E(outB[[k]])$lty <- ifelse(weight > 0, "solid", "dashed")
        E(outB[[k]])$width <- as.numeric(as.character(cut(abs(weight), 
                                                          breaks = c(0, .0001, .01, 1), 
                                                          labels = c(.5, 1.5, 2.5))))
        E(outB[[k]])$weight <- weight
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
      do.call(function(...) plot(gr[[j]], vertex.color=vcol, edge.color=ecol,...), opt)
    }
  } 
  else {
    #    if(is.null(opt$layout)) opt$layout <- layout_with_kk
    for(k in which){ do.call(function(...) plot(gr[[k]], ...), opt) }  
  }
  invisible(NULL)
}


##### to_graph auxiliary function #####
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
