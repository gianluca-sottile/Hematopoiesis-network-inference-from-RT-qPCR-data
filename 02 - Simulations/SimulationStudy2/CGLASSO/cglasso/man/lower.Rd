\name{lower + upper}
\alias{lower}
\alias{upper}

\title{Lower and Upper Limits from a \dQuote{datacggm} Object}

\description{
Functions \sQuote{\code{lower}} and \sQuote{\code{upper}} retrieve the lower and upper censoring values from an object of class \dQuote{datacggm}.
}

\usage{
lower(x)
upper(x)
}

\arguments{
\item{x}{an object of class \sQuote{\code{datacggm}}.}
}

\details{
For an \R object \code{x} of class \sQuote{\code{datacggm}}, \code{lower} (\code{upper}) retrieves the lower (upper) censoring values of the response variables.
}

\author{Luigi Augugliaro (\email{luigi.augugliaro@unipa.it})}

\seealso{
\code{\link{datacggm}}
}

\examples{
set.seed(123)
n <- 100L
p <- 3L
rho <- 0.3
Sigma <- outer(1L:p, 1L:p, function(i, j) rho^abs(i - j))

Z <- rcggm(n = n, Sigma = Sigma, probr = 0.05)
lower(Z)
upper(Z)

Z <- rcggm(n = n, Sigma = Sigma, probl = 0.05)
lower(Z)
upper(Z)

Z <- rcggm(n = n, Sigma = Sigma, probl = 0.05, probr = 0.05)
lower(Z)
upper(Z)
}

\keyword{array}
\keyword{manip}
