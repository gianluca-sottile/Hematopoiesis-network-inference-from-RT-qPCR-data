\name{to_graph}
\alias{to_graph}

\title{Create Graphs from cglasso or cggm Objects}

\description{
\sQuote{\code{to_graph}} returns a named list of graphs using the results of an \R object of class \sQuote{\code{cglasso}} or \sQuote{\code{cggm}}. 
}

\usage{
to_graph(object, GoF = AIC, lambda.id, rho.id, weighted = FALSE, simplify = TRUE,
         ...)
}

\arguments{
\item{object}{an \R object inheriting class \sQuote{\code{cglasso}}, that is, the output of the model-fitting functions \code{\link{cglasso}} and \code{\link{cggm}}.}
\item{GoF}{a valid goodness-of-fit function, such as \code{\link{AIC.cglasso}} or \code{\link{BIC.cglasso}}, or an \R object of class \sQuote{\code{GoF}}.}
\item{lambda.id}{optional. If \code{object} has class \sQuote{\code{cglasso}}, this argument is an integer used to identify a specific fitted cglasso model, otherwise (i.e. if the \code{object} has class \sQuote{\code{cggm}}) this argument can be omitted. See section \sQuote{Details} for more details.}
\item{rho.id}{optional. If \code{object} has class \sQuote{\code{cglasso}}, this argument is an integer used to identify a specific fitted cglasso model, otherwise (i.e. if \code{object} has class \sQuote{\code{cggm}}) this argument can be omitted. See section \sQuote{Details} for more details.}
\item{weighted}{logical. Should weighted graphs be created? Default is \code{FALSE}.}
\item{simplify}{logical. Should isolated vertices be removed from the graph? Default is \code{TRUE}, i.e., isolated vertices are removed.}
\item{\ldots}{further arguments passed to the chosen goodness-of-fit function (argument \sQuote{\code{GoF}}).}
}

\details{
\sQuote{\code{to_graph}} returns a named list of graphs using the results of an \R object of class \sQuote{\code{cglasso}} or \sQuote{\code{cggm}}. 

If \code{object} has class \sQuote{\code{cglasso}}, then the goodness-of-fit function passed through the argument \code{GoF} is used to identify the adjacency matrix (\code{object$InfoStructure$Adj_yy}) describing the undirected edges among the \eqn{p}{p} response variables. If the model is fitted using \eqn{q}{q} predictors, then the matrix  describing the effects of the predictors onto the response variables (see object$InfoStructure$Adj_xy) is also returned. Finally, these matrices are used to return an undirected and directed graph. Opionally, the user can identify a specific fitted model using the arguments \code{lambda.id} and \code{rho.id}.

If \code{object} has class \sQuote{\code{cggm}}, then \code{GoF}, \code{lambda.id} and \code{rho.id} can be omitted.

If argument weighted is set equal to \sQuote{\code{TRUE}}, then the estimated precision matrix and, if available, the estimated regression coefficient matrix are used to return weighted graphs. In this case, edges associated with positive estimates are shown using a solid line. Otherwise, a dashed line is used.
}

\value{
\sQuote{\code{to_graph}} returns an \R object of S3 class \dQuote{\code{cglasso2igraph}}, i.e., a named list containing the
following components:
\item{Gyy}{an undirected graph representing the conditional dependence structure among the \eqn{p}{p} response variables.}
\item{Gxy}{a directed graph representing the effetcs of the \eqn{q}{q} predictors onto the \eqn{p}{p} response variables.}
Each component is an \R object of class \code{\link[igraph]{igraph}}.
}

\author{Luigi Augugliaro (\email{luigi.augugliaro@unipa.it})}

\seealso{
\code{\link{cglasso}}, \code{\link{cggm}} and \code{\link{plot.cglasso2igraph}}. For more details about the object of class \sQuote{\code{igraph}}, the interested reader is referred to the package \pkg{\link[igraph]{igraph}}.
}

\examples{
set.seed(123)
# Y ~ N(0, Sigma) and probability of left/right censored values equal to 0.05
n <- 100L
p <- 3L
rho <- 0.3
Sigma <- outer(1L:p, 1L:p, function(i, j) rho^abs(i - j))
Z <- rcggm(n = n, Sigma = Sigma, probl = 0.05, probr = 0.05)
out <- cglasso(. ~ ., data = Z)
out.graph <- to_graph(out)
out.graph

# Y ~ N(b0 + XB, Sigma)  and probability of left/right censored values equal to 0.05
n <- 100L
p <- 3L
q <- 2L
b0 <- runif(p)
B <- matrix(runif(q * p), nrow = q, ncol = p)
X <- matrix(rnorm(n * q), nrow = n, ncol = q)
rho <- 0.3
Sigma <- outer(1L:p, 1L:p, function(i, j) rho^abs(i - j))
Z <- rcggm(n = n, b0 = b0, X = X, B = B, Sigma = Sigma, probl = 0.05, probr = 0.05)
out <- cglasso(. ~ ., data = Z)
out.graph <- to_graph(out, lambda.id = 3, rho.id = 3, weighted = TRUE)
out.graph
}

\keyword{graphs}
