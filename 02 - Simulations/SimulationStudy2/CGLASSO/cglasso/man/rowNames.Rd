\name{rowNames + colNames}
\alias{rowNames}
\alias{colNames}
\alias{rowNames<-}
\alias{colNames<-}

\title{Row and Column Names of a \dQuote{datacggm} Object}

\description{
Retrieve or set the row or column names of a \dQuote{datacggm} object.
}

\usage{
rowNames(x)
rowNames(x) <- value

colNames(x)
colNames(x) <- value
}

\arguments{
\item{x}{an \R object of class \sQuote{\code{datacggm}}.}
\item{value}{a named list with elements \sQuote{\code{X}} and \sQuote{\code{Y}}.}
}

\details{
For an \R object of class \sQuote{\code{datacggm}}, \code{rowNames} (\code{colNames}) retrieves or set \code{rownames} (\code{colnames}) for matrices \code{Y} and \code{X}. See below for some examples.
}

\author{Luigi Augugliaro (\email{luigi.augugliaro@unipa.it})}

\seealso{
\code{\link{datacggm}}, \code{\link{rcggm}} and the method function \code{\link{dimnames}}.
}

\examples{
set.seed(123)
# a dataset from a censored Gaussian graphical model
n <- 100L
p <- 3L
Z <- rcggm(n = n, p = p, probl = 0.05, probr = 0.05)

rowNames(Z)
rowNames(Z) <- list(Y = paste0("i", seq_len(n))) 
# the same as rowNames(Z)$Y <- paste0("i", seq_len(n))
rowNames(Z)

colNames(Z)
colNames(Z) <- list(Y = paste("Y", 1:p, sep = ":"))
# the same as colNames(Z)$Y <- paste("Y", 1:p, sep = ":")
colNames(Z)

# a dataset from a  conditional censored Gaussian graphical model
n <- 100L
p <- 3L
q <- 2
b0 <- runif(p)
B <- matrix(runif(q * p), nrow = q, ncol = p)
X <- matrix(rnorm(n * q), nrow = n, ncol = q)
Z <- rcggm(n = n, b0 = b0, X = X, B = B, probl = 0.05, probr = 0.05)

rowNames(Z)$Y <- paste0("i", seq_len(n))
rowNames(Z)$X <- paste0("i", seq_len(n))
dimnames(Z)

colNames(Z)$Y <- paste("Y", 1:p, sep = ":")
colNames(Z)$X <- paste("X", 1:q, sep = ":")
dimnames(Z)
}

\keyword{array}
\keyword{manip}
