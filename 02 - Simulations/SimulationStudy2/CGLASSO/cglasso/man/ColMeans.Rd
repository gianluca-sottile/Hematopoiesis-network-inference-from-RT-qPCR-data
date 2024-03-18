\name{ColMeans + ColMeans}
\alias{ColMeans}
\alias{ColVars}

\title{Calculate Column Means and Vars of a \dQuote{datacggm} Object}

\description{
Retrieve column means and column variances of a \dQuote{datacggm} object.
}

\usage{
ColMeans(x)
ColVars(x)
}

\arguments{
\item{x}{an object of class \sQuote{\code{datacggm}}.}
}

\details{
For an \R object \code{x} of class \sQuote{\code{datacggm}}, \code{ColMeans} (\code{ColVars}) retrieves the column means (variances) of the matrices obtained by \code{getMatrix(x, "Y")} and \code{getMatrix(x, "X")}. For the response variables, marginal means and variances are estimated using a EM-algorithm under the assumption that the \eqn{p}{p} response variables are marginally normally distributed (see also Details section in \code{\link{datacggm}}). For the numeric predictor variables, marginal means and variances are computed by \code{\link{mean}} and \code{\link{var}}, whereas, for categorical data, \code{ColMeans} (\code{ColVars}) retrieves the statistical mode and the Gini-Simpson Index, respectively.
}

\value{
\code{ColMeans} (\code{ColVars}) returns a named list with the columns means (variances).
}

\author{Luigi Augugliaro (\email{luigi.augugliaro@unipa.it})}

\seealso{
\code{\link{datacggm}}, \code{\link{rcggm}}, \code{\link{qqcnorm}} and \code{\link{hist.datacggm}}.
}

\examples{
set.seed(123)
n <- 1000L
p <- 3L
b0 <- rep(0, p)
Z <- rcggm(n = n, b0 = b0, probl = 0.05, probr = 0.05)
ColMeans(Z)
ColVars(Z)

n <- 1000L
p <- 3L
q <- 2
b0 <- runif(p)
B <- matrix(runif(q * p), nrow = q, ncol = p)
X <- matrix(rnorm(n * q), nrow = n, ncol = q)
Z <- rcggm(n = n, b0 = b0, X = X, B = B, probl = 0.05, probr = 0.05)
ColMeans(Z)
ColVars(Z)
}
\keyword{array}
\keyword{manip}
