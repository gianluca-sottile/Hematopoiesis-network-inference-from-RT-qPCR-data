is.cglasso2igraph <- function(x) inherits(x, 'cglasso2igraph')

print.cglasso2igraph <- function(x, ...) print.listof(x, ...)

plot.cglasso2igraph <- function(x, type, ...){
    if(missing(type)) type <- ifelse(is.null(x$Gxy), "Gyy", "both")
    gr <- getGraph(x, type = type)
    opt <- list(...)
#    if(is.null(opt$layout)) opt$layout <- layout_with_kk
    do.call(function(...) plot(gr, ...), opt)
    invisible(NULL)
}
