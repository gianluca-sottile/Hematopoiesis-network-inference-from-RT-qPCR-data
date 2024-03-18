\name{is.cglasso2igraph}
\alias{is.cglasso2igraph}

\title{Is an Object of Class \sQuote{cglasso2igraph}?}

\usage{
is.cglasso2igraph(x)
}

\description{
\code{is.cglasso2igraph} tests if its argument is an object of class \sQuote{cglasso2igraph}.
}

\arguments{
  \item{x}{object to be tested.}
}

\author{Luigi Augugliaro (\email{luigi.augugliaro@unipa.it})}

\seealso{
\code{\link{to_graph}}.
}

\examples{
set.seed(123)
n <- 100L
p <- 3L
rho <- 0.3
Sigma <- outer(1L:p, 1L:p, function(i, j) rho^abs(i - j))
Z <- rcggm(n = n, Sigma = Sigma, probl = 0.05, probr = 0.05)
out <- cglasso(. ~ ., data = Z)
out.graph <- to_graph(out)
is.cglasso2igraph(out.graph )
}

\keyword{methods}
\keyword{classes}
