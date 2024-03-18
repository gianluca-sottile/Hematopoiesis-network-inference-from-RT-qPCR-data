getGraph <- function(x, type = c("Gyy", "Gxy", "both")){
    if (!is.cglasso2igraph(x))
        stop(sQuote(deparse(substitute(x))), " is not an object of class", sQuote("cglasso2igraph"))
    type <- match.arg(type)
    if (is.null(x$Gxy) & is.element(type, c("Gxy", "both")))
        stop(sQuote(type), " is not available. Please, use type = ", dQuote("Gyy"))
    if (type == "Gyy") gr <- x$Gyy
    if (type == "Gxy") gr <- x$Gxy
    if (type == "both") {
        gr.e <- rbind(as_data_frame(x$Gyy, what = "edges"), as_data_frame(x$Gxy, what = "edges"))
        v.yy <- as_data_frame(x$Gyy, what = "vertices")
        v.xy <- as_data_frame(x$Gxy, what = "vertices")
        id <- !is.element(v.xy$name, v.yy$name)
        gr.v <- rbind(v.yy, v.xy[id, ])
        gr <- graph_from_data_frame(d = gr.e, directed = TRUE, vertices = gr.v)
    }
    gr
}
