\name{plot.cglasso2igraph}
\alias{plot.cglasso2igraph}

\title{Plot Method for a cglasso2igraph Object"}

\description{
\code{plot.cglasso2igraph} produces graphs from an \R object of class \sQuote{\code{cglasso2igraph}}.
}

\usage{
\method{plot}{cglasso2igraph}(x, type, \dots)
}

\arguments{
\item{x}{an \R object of class \sQuote{\code{cglasso2igraph}}, that is, the output of the function \code{\link{to_graph}}.}
\item{type}{a description of the required graph. Default is \sQuote{\code{both}} in a conditional glasso estimator or \sQuote{\code{Gyy}} in a glasso estimator (see \code{\link{getGraph}}).}
\item{\dots}{additional graphical arguments passed to the functions \code{\link[igraph]{plot.igraph}}.}
}

\author{Luigi Augugliaro (\email{luigi.augugliaro@unipa.it})}

\seealso{
\code{\link{to_graph}} and \code{\link{getGraph}}.
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
plot(out.graph, type = "Gyy")

out.graph <- to_graph(out, weighted = TRUE)
plot(out.graph,  type = "Gyy")


# Y ~ N(b0 +XB, Sigma)  and probability of left/right censored values equal to 0.05
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
plot(out.graph, type = "Gyy")
plot(out.graph, type = "Gxy")
plot(out.graph, type = "both")
}

\keyword{graphs}
