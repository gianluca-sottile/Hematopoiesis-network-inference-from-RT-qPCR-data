\name{event}
\alias{event}

\title{Status Indicator Matrix from a \sQuote{\code{datacggm}} Object}

\description{
The \sQuote{\code{event}} function retrieves the status indicator matrix from an object of class \sQuote{\code{datacggm}}.
}

\usage{
event(x, ordered = FALSE)
}

\arguments{
\item{x}{an object of class \sQuote{\code{datacggm}}.}
\item{ordered}{logical value used to specify if the rows of the status indicator matrix should be ordered according to the patterns of censored/missing values. Default \code{ordered = FALSE}.}
}

\details{
The \sQuote{\code{event}} function is used to retrieve the status indicator matrix, denoted by \code{R}, from an object of class \sQuote{\code{datacggm}}. The entries of the matrix are used to specify the status of the response variable:
\itemize{
\item \sQuote{\code{R[i, j] =  0}} means that \eqn{y_{ij}}{y_{ij}} is inside the open interval \code{(lo[j], up[j])};
\item \sQuote{\code{R[i, j] = -1}} means that \eqn{y_{ij}}{y_{ij}} is a left-censored value;
\item \sQuote{\code{R[i, j] = +1}} means that \eqn{y_{ij}}{y_{ij}} is a right-censored value;
\item \sQuote{\code{R[i, j] = +9}} means that \eqn{y_{ij}}{y_{ij}} is a missing value.
}
}

\value{
\code{event} returns a \eqn{(n\times p)}{(n x p)}-dimensional matrix.
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
\code{\link{datacggm}} and \code{\link{rcggm}}.
}

\examples{
set.seed(123)

# Y ~ N(b0 + XB, Sigma) and
# 1. probability of left/right censored values equal to 0.05
# 2. probability of missing-at-random euqals to 0.05
n <- 100L
p <- 3L
q <- 2
b0 <- runif(p)
B <- matrix(runif(q * p), nrow = q, ncol = p)
X <- matrix(rnorm(n * q), nrow = n, ncol = q)
rho <- 0.3
Sigma <- outer(1L:p, 1L:p, function(i, j) rho^abs(i - j))
Z <- rcggm(n = n, b0 = b0, X = X, B = B, Sigma = Sigma, probl = 0.05, 
           probr = 0.05, probna = 0.05)

# status indicator matrix
event(Z)

# in this case the status indicator matrix is returned with
# rows ordered according to the patterns of missing data
event(Z, ordered = TRUE)
}

\keyword{classes}
\keyword{methods}
