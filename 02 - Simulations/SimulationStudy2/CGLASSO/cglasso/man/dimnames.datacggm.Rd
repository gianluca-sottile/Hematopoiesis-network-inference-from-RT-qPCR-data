\name{dimnames.datacggm}
\alias{dimnames}
\alias{dimnames.datacggm}
\alias{dimnames<-.datacggm}

\title{Dimnames of a \dQuote{datacggm} Object}

\description{
Retrieve or set the dimnames of a \dQuote{datacggm} object.
}

\usage{
\method{dimnames}{datacggm}(x)

\method{dimnames}{datacggm}(x) <- value
}

\arguments{
\item{x}{an object of class \sQuote{\code{datacggm}}.}
\item{value}{a nested list with names \sQuote{\code{X}} and \sQuote{\code{Y}}. See below for mode details.}
}

\details{
For an \R object \sQuote{\code{x}} of class \sQuote{\code{datacggm}}, \code{dimnames} retrieves or sets the \code{dimnames} attribute for the matrices \code{Y} and \code{X} (see \sQuote{\code{\link{datacggm}}} for more details). When setting the \code{dimnames} attribute,  \code{value$Y} can be \code{NULL} (which is not stored) or a list of length two. In the last case, \code{value$Y} is passed to the method function  \code{\link[base]{dimnames}} for setting the \code{dimnames} attributes of the matrix retrieved by \code{getMatrix(x, "Y")}. In the same way, \code{value$X} can be used for setting the \code{dimnames} attributes of the matrix retrieved by \code{getMatrix(x, "X")}.
}

\author{Luigi Augugliaro (\email{luigi.augugliaro@unipa.it})}

\seealso{
\code{\link{datacggm}}, \code{\link{rcggm}}, \code{\link{getMatrix}}, \code{\link{rowNames}} and \code{\link{colNames}}. 
}

\examples{
set.seed(123)
# a dataset from a censored Gaussian graphical model
n <- 100L
p <- 3L
b0 <- rep(0, p)
Z <- rcggm(n = n, b0 = b0, probl = 0.05, probr = 0.05)
dimnames(Z)

dimnames(Z) <- list(Y = list(paste0("i", seq_len(n)), paste("Y", 1:p, sep = ":")))
dimnames(Z)

# the same as
# dimnames(Z)$Y <- list(paste0("i", seq_len(n)), paste("Y", 1:p, sep = ":"))
# dimnames(Z)


# a dataset from a  conditional censored Gaussian graphical model
n <- 100L
p <- 3L
q <- 2
b0 <- runif(p)
B <- matrix(runif(q * p), nrow = q, ncol = p)
X <- matrix(rnorm(n * q), nrow = n, ncol = q)
Z <- rcggm(n = n, b0 = b0, X = X, B = B, probl = 0.05, probr = 0.05)
dimnames(Z)

dimnames(Z) <- list(Y = list(NULL, paste("Y", 1:p, sep = ":")), 
                    X = list(NULL, paste("X", 1:q, sep = ":")))
dimnames(Z)

# the same as
# dimnames(Z)$Y <- list(NULL, paste("Y", 1:p, sep = ":"))
# dimnames(Z)$X <- list(NULL, paste("X", 1:q, sep = ":"))
# dimnames(Z)
}

\keyword{array}
\keyword{manip}


