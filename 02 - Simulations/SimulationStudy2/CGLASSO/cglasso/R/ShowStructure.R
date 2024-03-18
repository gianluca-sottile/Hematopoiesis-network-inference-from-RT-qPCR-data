ShowStructure <- function(module = c("ALL", "DM", "MF", "MS", "NA"), description = TRUE, plot = TRUE) {
    if (!is.logical(description)) stop(sQuote(description), " is not a logical value")
    if (!is.logical(plot)) stop(sQuote(plot), " is not a logical value")
    if (!description & !plot) stop(sQuote("description"), " and ", sQuote("plot"), " are both equal to FALSE")
    module <- match.arg(module)
    cglasso.graph <- graph_from_literal("datacggm" -- "rcggm",
                                        "datacggm" -- "is.datacggm",
                                        "datacggm" -- "dim.datacggm",
                                        "datacggm" -- "nobs/nresp/npred",
                                        "datacggm" -- "rowNames/colNames",
                                        "datacggm" -- "dimnames.datacggm",
                                        "datacggm" -- "print.datacggm",
                                        "datacggm" -- "summary.datacggm",
                                        "datacggm" -- "getMatrix",
                                        "datacggm" -- "ColMeans/ColVars",
                                        "datacggm" -- "lower/upper",
                                        "datacggm" -- "event",
                                        "datacggm" -- "hist.datacggm",
                                        "datacggm" -- "qqcnorm",
                                        "datacggm" -- "cglasso",
                                        "cglasso" -- "print.cglasso",
                                        "cglasso" -- "plot.cglasso",
                                        "cglasso" -- "coef.cglasso",
                                        "cglasso" -- "fitted.cglasso",
                                        "cglasso" -- "residuals.cglasso",
                                        "cglasso" -- "QFun",
                                        "cglasso" -- "AIC.cglasso/BIC.cglasso",
                                        "cglasso" -- "summary.cglasso",
                                        "cglasso" -- "select.cglasso",
                                        "cglasso" -- "predict.cglasso",
                                        "cglasso" -- "impute",
                                        "cglasso" -- "cggm",
                                        "cglasso" -- "to_graph",
                                        "cggm" -- "print.cggm",
                                        "cggm" -- "coef.cglasso",
                                        "cggm" -- "fitted.cglasso",
                                        "cggm" -- "residuals.cglasso",
                                        "cggm" -- "predict.cggm",
                                        "cggm" -- "impute",
                                        "cggm" -- "QFun",
                                        "cggm" -- "plot.cggm",
                                        "cggm" -- "AIC.cglasso/BIC.cglasso",
                                        "cggm" -- "summary.cglasso",
                                        "cggm" -- "to_graph",
                                        "QFun" -- "print.QFun",
                                        "QFun" -- "AIC.cglasso/BIC.cglasso",
                                        "AIC.cglasso/BIC.cglasso" -- "print.GoF",
                                        "AIC.cglasso/BIC.cglasso" -- "plot.GoF",
                                        "AIC.cglasso/BIC.cglasso" -- "summary.cglasso",
                                        "AIC.cglasso/BIC.cglasso" -- "select.cglasso",
                                        "AIC.cglasso/BIC.cglasso" -- "plot.cglasso",
                                        "to_graph" -- "is.cglasso2igraph",
                                        "to_graph" -- "getGraph",
                                        "to_graph" -- "print.cglasso2igraph",
                                        "to_graph" -- "plot.cglasso2igraph")
                                        
    data.manipulation <- list(
        "datacggm" = "Create a Dataset from a Conditional Gaussian Graphical", "  Model with Censored and/or Missing Values",
        "is.datacggm" = "Is an Object of Class datacggm?",
        "print.datacggm" = "Print Method for a datacggm Object",
        "summary.datacggm" = "Summarizing Objects of Class datacggm",
        "dim.datacggm" = "Dimensions of a datacggm Object",
        "nobs/nresp/npred" = "Extract the Number of Observations/Responses/Predictors", "  from a datacggm Object",
        "dimnames.datacggm" = "Dimnames of a datacggm Object",
        "rowNames/colNames" = "Row and Column Names of a datacggm Object",
        "getMatrix" = "Retrieve Matrices Y and X from a datacggm Object",
        "event" = "Status Indicator Matrix from a datacggm Object",
        "lower/upper" = "Lower and Upper Limits from a datacggm Object",
        "ColMeans/ColVars" = "Form Column Means and Vars of a datacggm Object",
        "rcggm" = "Simulate from a Conditional Gaussian Graphical Model", "  with Censored and/or Missing Values",
        "hist.datacggm" = "Histogram for a datacggm Object",
        "qqcnorm" = "Quantile-Quantile Plots for a datacggm Object")
    class(data.manipulation) <- "simple.list"

    model.fitting <- list(
        "cglasso" = "Conditional Graphical Lasso Estimator",
        "print.cglasso" = "Print Method for a cglasso Object",
        "plot.cglasso" = "Plot Method for a cglasso Object",
        "coef.cglasso" = "Extract Model Coefficients",
        "fitted.cglasso" = "Extract Model Fitted Values",
        "residuals.cglasso" = "Extract Model Residuals",
        "predict.cglasso" = "Predict Method for cglasso Fits",
        "impute" = "Imputation of Missing and Censored Values",
        "cggm" = "Post-Hoc Maximum Likelihood Refitting", "  of a Conditional Graphical Lasso",
        "print.cggm" = "Print Method for a cggm Object",
        "plot.cggm" = "Plot Method for a cggm Object",
        "predict.cggm" = "Predict Method for cggm Fits")
    class(model.fitting) <- "simple.list"

    model.selection <- list(
        "QFun" = "Extract Q-Function",
        "print.QFun" = "Print Method for a QFun Object",
        "AIC.cglasso/BIC.cglasso" = "Goodness-of-fit Functions",
        "print.GoF" = "Print Method for a GoF Object",
        "summary.cglasso" = "Summarizing cglasso and cggm Fits",
        "plot.GoF" = "Plot Method for a GoF Object",
        "select.cglasso" = "Model Selection for Conditional", "  Graphical Lasso Estimator")
    class(model.selection) <- "simple.list"

    network.analysis <- list(
        "to_graph" = "Create Graphs from cglasso or cggm Objects",
        "is.cglasso2igraph" = "Is an Object of Class cglasso2igraph?",
        "print.cglasso2igraph" = "Print Method for a cglasso2igraph Object",
        "getGraph" = "Retrieve Graphs from a cglasso2igraph Object",
        "plot.cglasso2igraph" = "Plot Method for a cglasso2igraph Object")
    class(network.analysis) <- "simple.list"

    cglasso.description <- list(
      "Data Manipulation" = data.manipulation,
      "Model Fitting" = model.fitting,
      "Model Selection" = model.selection,
      "Network Analysis" = network.analysis)

      
#    n <- length(data.manipulation) + length(model.fitting) + length(model.selection) + length(network.analysis)
    data.manipulation.nm <- setdiff(names(data.manipulation), "")
    model.fitting.nm <- setdiff(names(model.fitting), "")
    model.selection.nm <- setdiff(names(model.selection), "")
    network.analysis.nm <- setdiff(names(network.analysis), "")
    n <- length(data.manipulation.nm) + length(model.fitting.nm) + length(model.selection.nm) + length(network.analysis.nm)
    
    V(cglasso.graph)$size <- 2
    V(cglasso.graph)$label.cex <- 0.8
    V(cglasso.graph)$label.dist <- 0.9
    
    V(cglasso.graph)[data.manipulation.nm]$label.color <-  "red4"
    V(cglasso.graph)[model.fitting.nm]$label.color <-  "blue4"
    V(cglasso.graph)[model.selection.nm]$label.color <- "darkgoldenrod4"
    V(cglasso.graph)[network.analysis.nm]$label.color <- "darkolivegreen"
    
    V(cglasso.graph)[data.manipulation.nm]$color <-  "red4"
    V(cglasso.graph)[model.fitting.nm]$color <-  "blue4"
    V(cglasso.graph)[model.selection.nm]$color <- "darkgoldenrod4"
    V(cglasso.graph)[network.analysis.nm]$color <- "darkolivegreen"

    V(cglasso.graph)[data.manipulation.nm]$frame.color <-  "red4"
    V(cglasso.graph)[model.fitting.nm]$frame.color <-  "blue4"
    V(cglasso.graph)[model.selection.nm]$frame.color <- "darkgoldenrod4"
    V(cglasso.graph)[network.analysis.nm]$frame.color <- "darkolivegreen"

    V(cglasso.graph)$label.degree <- rep(- pi / 4, n)
    # data manipulation module
    V(cglasso.graph)["datacggm"]$label.degree <- - 0.1
    V(cglasso.graph)["datacggm"]$label.dist <- + 2.1
    V(cglasso.graph)["is.datacggm"]$label.degree <- 0
    V(cglasso.graph)["is.datacggm"]$label.dist <- - 2.4
    V(cglasso.graph)["dimnames.datacggm"]$label.degree <- 0
    V(cglasso.graph)["dimnames.datacggm"]$label.dist <- 3.7
    V(cglasso.graph)["print.datacggm"]$label.degree <- 0
    V(cglasso.graph)["print.datacggm"]$label.dist <- 3
    V(cglasso.graph)["nobs/nresp/npred"]$label.degree <- 0
    V(cglasso.graph)["nobs/nresp/npred"]$label.dist <- 3.5
    V(cglasso.graph)["summary.datacggm"]$label.degree <- 0
    V(cglasso.graph)["summary.datacggm"]$label.dist <- - 3.5
    V(cglasso.graph)["dim.datacggm"]$label.degree <- 0
    V(cglasso.graph)["dim.datacggm"]$label.dist <- - 2.8
    V(cglasso.graph)["ColMeans/ColVars"]$label.degree <- 0
    V(cglasso.graph)["ColMeans/ColVars"]$label.dist <- - 3.5
    V(cglasso.graph)["rcggm"]$label.degree <- 0
    V(cglasso.graph)["rcggm"]$label.dist <- - 1.5
    V(cglasso.graph)["getMatrix"]$label.degree <- 0
    V(cglasso.graph)["getMatrix"]$label.dist <- - 2.0
    V(cglasso.graph)["event"]$label.degree <- 0
    V(cglasso.graph)["event"]$label.dist <- - 1.5
    V(cglasso.graph)["qqcnorm"]$label.degree <- 0.3
    V(cglasso.graph)["qqcnorm"]$label.dist <- + 1.7
    V(cglasso.graph)["hist.datacggm"]$label.degree <- 0
    V(cglasso.graph)["hist.datacggm"]$label.dist <- - 2.5

    # Model fitting
    V(cglasso.graph)["cglasso"]$label.degree <- 0
    V(cglasso.graph)["cglasso"]$label.dist <- - 1.7
    V(cglasso.graph)["plot.cglasso"]$label.degree <- 3.4
    V(cglasso.graph)["plot.cglasso"]$label.dist <- + 2
    V(cglasso.graph)["impute"]$label.degree <- + 2.5
    V(cglasso.graph)["impute"]$label.dist <- 1.0
    V(cglasso.graph)["fitted.cglasso"]$label.degree <- 0
    V(cglasso.graph)["fitted.cglasso"]$label.dist <- - 2.5
    V(cglasso.graph)["predict.cglasso"]$label.degree <- 0
    V(cglasso.graph)["predict.cglasso"]$label.dist <- - 2.8
    V(cglasso.graph)["print.cglasso"]$label.degree <- + pi / 4
    V(cglasso.graph)["plot.cggm"]$label.degree <- 0
    V(cglasso.graph)["plot.cggm"]$label.dist <- + 2.0
    V(cglasso.graph)["predict.cggm"]$label.degree <- 0
    V(cglasso.graph)["predict.cggm"]$label.dist <- + 2.5
    V(cglasso.graph)["print.cggm"]$label.degree <- 0
    V(cglasso.graph)["print.cggm"]$label.dist <- 2.2
    V(cglasso.graph)["cggm"]$label.degree <- - 0.2
    V(cglasso.graph)["cggm"]$label.dist <- + 1.5
    
    # Model selection
    V(cglasso.graph)["plot.GoF"]$label.degree <- 0
    V(cglasso.graph)["plot.GoF"]$label.dist <- + 2.0
    V(cglasso.graph)["print.GoF"]$label.degree <- 0
    V(cglasso.graph)["print.GoF"]$label.dist <- + 2.0
    V(cglasso.graph)["print.QFun"]$label.degree <- 0
    V(cglasso.graph)["print.QFun"]$label.dist <- + 2.0
    V(cglasso.graph)["AIC.cglasso/BIC.cglasso"]$label.degree <- 0.1
    V(cglasso.graph)["AIC.cglasso/BIC.cglasso"]$label.dist <- + 4.3
    V(cglasso.graph)["QFun"]$label.degree <- 0.4
    V(cglasso.graph)["QFun"]$label.dist <- + 1.3
    V(cglasso.graph)["summary.cglasso"]$label.degree <- 0
    V(cglasso.graph)["summary.cglasso"]$label.dist <- + 3.0
    
    # Network analysis
    V(cglasso.graph)["to_graph"]$label.degree <- 3.4
    V(cglasso.graph)["to_graph"]$label.dist <- + 1.8
    V(cglasso.graph)["is.cglasso2igraph"]$label.degree <- 0
    V(cglasso.graph)["is.cglasso2igraph"]$label.dist <- + 3.0
    V(cglasso.graph)["plot.cglasso2igraph"]$label.degree <- 0
    V(cglasso.graph)["plot.cglasso2igraph"]$label.dist <- - 3.5
    V(cglasso.graph)["print.cglasso2igraph"]$label.degree <- 0
    V(cglasso.graph)["print.cglasso2igraph"]$label.dist <- - 3.5
    V(cglasso.graph)["getGraph"]$label.degree <- + pi / 4
    V(cglasso.graph)["lower/upper"]$label.degree <- + pi / 4

    V(cglasso.graph)["datacggm"]$label.font <- 2
    V(cglasso.graph)["cglasso"]$label.font <- 2
    V(cglasso.graph)["cggm"]$label.font <- 2
    V(cglasso.graph)["QFun"]$label.font <- 2
    V(cglasso.graph)["AIC.cglasso/BIC.cglasso"]$label.font <- 2
    V(cglasso.graph)["to_graph"]$label.font <- 2

    legend.txt <- c("Data manipulation", "Model fitting", "Model selection", "Network analysis")
    legend.col <- c("red4", "blue4", "darkgoldenrod4", "darkolivegreen")
    
    if (module == "DM") {
        vids <- which(is.element(names(V(cglasso.graph)), data.manipulation.nm))
        legend.txt <- legend.txt[1L]
        legend.col <- legend.col[1L]
        cglasso.description <- cglasso.description["Data Manipulation"]
    }
    if (module == "MF") {
        vids <- which(is.element(names(V(cglasso.graph)), model.fitting.nm))
        legend.txt <- legend.txt[2L]
        legend.col <- legend.col[2L]
        cglasso.description <- cglasso.description["Model Fitting"]
    }
    if (module == "MS") {
        vids <- which(is.element(names(V(cglasso.graph)), model.selection.nm))
        legend.txt <- legend.txt[3L]
        legend.col <- legend.col[3L]
        cglasso.description <- cglasso.description["Model Selection"]
    }
    if (module == "NA") {
        vids <- which(is.element(names(V(cglasso.graph)), network.analysis.nm))
        legend.txt <- legend.txt[4L]
        legend.col <- legend.col[4L]
        cglasso.description <- cglasso.description["Network Analysis"]
    }
    if (module != "ALL")
        cglasso.graph <- induced_subgraph(cglasso.graph, vids = vids)
    
    if (plot) {
        cglasso.graph$layout <- layout_(cglasso.graph, with_kk())
        plot(cglasso.graph, main = "cglasso Package")
        legend(x = -1.6, y = -0.75,
               legend = legend.txt,
               text.col = legend.col,
               cex = 0.8,
               bty = "n",
               border = NULL)
    }
    if (description) {
        print.listof(cglasso.description)
        cat("NOTE: use", sQuote("?method.class"), "to get the documentation pages\n\n")
    }
    out <- list(description = cglasso.description, graph = cglasso.graph)
    invisible(out)
}
