\name{nobs + nresp + npred}

\alias{nobs}
\alias{nobs.datacggm}
\alias{nresp}
\alias{nresp.datacggm}
\alias{npred}
\alias{npred.datacggm}

\title{Extract the Number of Observations/Responses/Predictors from a datacggm Object}

\usage{
\method{nobs}{datacggm}(object, ...)
\method{nresp}{datacggm}(object, ...)
\method{npred}{datacggm}(object, ...)
}

\description{
Extract the number of observations, response variables and predictors from an object of class \code{datacggm}.
}

\arguments{
  \item{object}{an \R object of class \code{datacggm}.}
  \item{...}{further arguments to be passed to methods.}
}

\author{Luigi Augugliaro (\email{luigi.augugliaro@unipa.it})}

\seealso{
\code{\link{datacggm}}, \code{\link{rcggm}} and \code{\link{dim.datacggm}}.
}

\examples{
set.seed(123)
n <- 100
p <- 3
q <- 2
b0 <- rep(1, p)
X <- matrix(rnorm(n * q), n, q)
B <- matrix(rnorm(q * p), q, p)
Sigma <- outer(1:p, 1:p, function(i, j) 0.3^abs(i - j))
probl <- 0.05
probr <- 0.05
probna <- 0.05

Z <- rcggm(n = n, b0 = b0, Sigma = Sigma, probl = probl, probr = probr, 
           probna = probna)
nobs(Z)
nresp(Z)
npred(Z)

Z <- rcggm(b0 = b0, X = X, B = B, Sigma = Sigma, probl = probl, probr = probr, 
           probna = probna)
nobs(Z)
nresp(Z)
npred(Z)
}

\keyword{methods}
\keyword{classes}
