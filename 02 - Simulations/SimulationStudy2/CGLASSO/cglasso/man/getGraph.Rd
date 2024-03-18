\name{getGraph}
\alias{getGraph}

\title{Retrieve Graphs from a \sQuote{\code{cglasso2igraph}} Object}

\description{
\sQuote{\code{getGraph}} retrieves graphs from an \R object of class \sQuote{\code{cglasso2igraph}}. 
}

\usage{
getGraph(x, type = c("Gyy", "Gxy", "both"))
}

\arguments{
\item{x}{an object of class \sQuote{\code{cglasso2igraph}} (see also \code{\link{to_graph}}).}
\item{type}{a description of the required graph. Default is \sQuote{\code{Gyy}}.}
}

\value{
\code{getGraph} retrieves an \R object of class \sQuote{\code{igraph}} representing the graph required by the argument \code{type.}
}

\author{Luigi Augugliaro (\email{luigi.augugliaro@unipa.it})}


\seealso{
\code{\link{to_graph}}, \code{\link{is.cglasso2igraph}} and \code{\link{plot.cglasso2igraph}}.
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
getGraph(out.graph)


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
getGraph(out.graph, type = "Gyy")
getGraph(out.graph, type = "Gxy")
getGraph(out.graph, type = "both")
}

\keyword{classes}
\keyword{methods}
