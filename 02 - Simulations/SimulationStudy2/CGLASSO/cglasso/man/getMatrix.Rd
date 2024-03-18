\name{getMatrix}
\alias{getMatrix}

\title{Retrieve Matrices \sQuote{\code{Y}} and \sQuote{\code{X}} from a \sQuote{\code{datacggm}} Object}

\description{
\sQuote{\code{getMatrix}} retrieves matrices \sQuote{\code{Y}} and/or \sQuote{\code{X}} from an object of class \sQuote{\code{datacggm}}. 
}

\usage{
getMatrix(x, name = c("Y", "X", "both"), ordered = FALSE)
}

\arguments{
\item{x}{an object of class \sQuote{\code{\link{datacggm}}}.}
\item{name}{the name of the required matrix.}
\item{ordered}{logical value used to specify if the required matrix should be retrieved with rows ordered according to the patterns of censored values. See below for some example.}
}

\value{
\code{getMatrix} retrieves the matrix specified by \sQuote{\code{name}} and with row ordering specified by \sQuote{\code{ordered}}. A named list returned if \code{name} is \code{"both"}.
}

\author{Luigi Augugliaro (\email{luigi.augugliaro@unipa.it})}

\references{
Augugliaro, L., Sottile, G., and Vinciotti, V. (2020) <\doi{10.1007/s11222-020-09945-7}>.
The conditional censored graphical lasso estimator.
\emph{Statistics and Computing} \bold{30}, 1273--1289.

Augugliaro, L., Abbruzzo, A., and Vinciotti, V. (2020) <\doi{10.1093/biostatistics/kxy043}>.
\eqn{\ell_1}{l1}-Penalized censored Gaussian graphical model.
\emph{Biostatistics} \bold{21}, e1--e16.
}

\seealso{
\code{\link{datacggm}}, \code{\link{rcggm}}.
}

\examples{
set.seed(123)

# a dataset from a  conditional censored Gaussian graphical model
n <- 100L
p <- 3L
q <- 2
b0 <- runif(p)
B <- matrix(runif(q * p), nrow = q, ncol = p)
X <- matrix(rnorm(n * q), nrow = n, ncol = q)
Z <- rcggm(n = n, b0 = b0, X = X, B = B, probl = 0.05, probr = 0.05)
getMatrix(Z, name = "Y")

# in the following examples 'Y' and 'X' is returned with rows ordered 
# according to the  patterns of censored data
getMatrix(Z, name = "Y", ordered = TRUE)

getMatrix(Z, name = "X")
getMatrix(Z, name = "X", ordered = TRUE)

getMatrix(Z, name = "both")
getMatrix(Z, name = "both", ordered = TRUE)
}

\keyword{classes}
\keyword{methods}
