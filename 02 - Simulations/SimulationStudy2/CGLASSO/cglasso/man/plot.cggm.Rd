\name{plot.cggm}
\alias{plot.cggm}

\title{Plot Method for a \sQuote{\code{cggm}} Object}

\description{
The \code{plot.cggm} produces graphs from an \R object of class \sQuote{\code{cggm}}.
}

\usage{
\method{plot}{cggm}(x, type, weighted = FALSE, simplify = TRUE, \dots)
}

\arguments{
\item{x}{an \R object of class \sQuote{\code{cggm}}, that is, the output of the function \code{\link{cggm}}.}
\item{type}{a description of the required graph. Default depends on the type of fitted model, that is, \code{type = both} if a conditional censored graphical lasso estimator is fitted, otherwise \code{type} is equal to \code{Gyy} (see \code{\link{getGraph}}) for more details.}
\item{weighted}{logical. Should weighted graphs be created? Default is \code{FALSE}.}
\item{simplify}{logical. Should isolated vertices be removed from the graph? Default is \code{TRUE}, i.e., isolated vertices are removed.}
\item{\dots}{additional graphical arguments passed to the functions \code{\link[igraph]{plot.igraph}}.}
}

\author{Luigi Augugliaro (\email{luigi.augugliaro@unipa.it})}

\seealso{
\code{\link{cggm}}, \code{\link{to_graph}}, \code{\link{getGraph}} and \code{\link{plot.cglasso2igraph}}.
}

\examples{
set.seed(123)
# Y ~ N(XB, Sigma)  and
# 1. probability of left/right censored values equal to 0.05
# 2. probability of missing-at-random values equal to 0.05
n <- 100L
p <- 3L
q <- 2L
b0 <- runif(p)
B <- matrix(runif(q * p), nrow = q, ncol = p)
X <- matrix(rnorm(n * q), nrow = n, ncol = q)
rho <- 0.3
Sigma <- outer(1L:p, 1L:p, function(i, j) rho^abs(i - j))
Z <- rcggm(n = n, b0 = b0, X = X, B = B, Sigma = Sigma, probl = 0.05, probr = 0.5, 
           probna = 0.05)
out <- cglasso(. ~ ., data = Z)
out.mle <- cggm(out, lambda.id = 3L, rho.id = 3L)
plot(out.mle, type = "Gyy")
plot(out.mle, type = "Gxy")
plot(out.mle, type = "both")
}

\keyword{graphs}
