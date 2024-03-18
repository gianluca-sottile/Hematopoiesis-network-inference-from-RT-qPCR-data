\name{dim.datacggm}
\alias{dim}
\alias{dim.datacggm}

\title{Dimensions of a \dQuote{datacggm} Object}

\description{
Retrieve the dimension of an \R object of class \dQuote{datacggm}, that is, the sample size (\eqn{n}{n}), the number of response variables (\eqn{p}{p}) and the number of predictors (\eqn{q}{q}). 
}

\usage{
\method{dim}{datacggm}(x)
}

\arguments{
\item{x}{an \R object of class \sQuote{\code{datacggm}}.}
}

\details{
For an \R object of class \sQuote{\code{datacggm}}, \code{dim} retrieves a list with elements named \code{Y} and \code{X}. The first component is the \dQuote{\code{dim}} attribute of the matrix \code{Y} whereas the second component is the \dQuote{\code{dim}} attribute of the matrix \code{X}. If \code{X} is missing, the second component is \code{NULL}.
}

\value{
For an object of class \sQuote{\code{datacggm}}, \code{dim} retrieves a named list with components \code{Y} and \code{X} which are the \dQuote{\code{dim}} attribute of the two matrices, respectively. If \code{X} is missing then the last component is \code{NULL}.
}

\author{Luigi Augugliaro (\email{luigi.augugliaro@unipa.it})}

\seealso{
\code{\link{datacggm}}, \code{\link{rcggm}},\code{\link{nobs}}, \code{\link{nresp}} and \code{\link{npred}}.
}

\examples{
set.seed(123)
# a dataset from a censored Gaussian graphical model
n <- 100L
p <- 3L
b0 <- rep(0, p)
Z <- rcggm(n = n, b0 = b0, probl = 0.05, probr = 0.05)
dim(Z)

# a dataset from a  conditional censored Gaussian graphical model
n <- 100L
p <- 3L
q <- 2
b0 <- runif(p)
B <- matrix(runif(q * p), nrow = q, ncol = p)
X <- matrix(rnorm(n * q), nrow = n, ncol = q)
Z <- rcggm(n = n, b0 = b0, X = X, B = B, probl = 0.05, probr = 0.05)
dim(Z)
}

\keyword{array}

