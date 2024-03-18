\name{fitted}
\alias{fitted}
\alias{fitted.cglasso}

\title{Extract Model Fitted Values}

\description{
The accessor function \code{fitted} extracts model fitted values from an \R object inheriting class \sQuote{\code{cglasso}}.
}

\usage{
\method{fitted}{cglasso}(object, lambda.id, rho.id, drop = TRUE, \dots)
}

\arguments{
\item{object}{an \R object inheriting class \sQuote{\code{cglasso}}, that is, the output of the model-fitting functions \code{\link{cglasso}} or \code{\link{cggm}}.}
\item{lambda.id}{a vector of integers used to specify the \eqn{\lambda}{lambda}-values.}
\item{rho.id}{a vector of integers used to specify the \eqn{\rho}{rho}-values.}
\item{drop}{logical. Dimensions can only be dropped if their extent is one.}
\item{\dots}{further arguments passed to or from other methods.}
}

\value{
Fitted values extracted from \sQuote{\code{object}} are returned.
}

\author{Luigi Augugliaro (\email{luigi.augugliaro@unipa.it})}

\seealso{
Model-fitting functions \code{\link{cglasso}}, \code{\link{cggm}} and the accessor functions \code{\link{coef.cglasso}}, \code{\link{residuals.cglasso}}, \code{\link{predict.cglasso}} and \code{\link{impute}}.
}

\examples{
set.seed(123)
# Y ~ N(0, Sigma) and probability of left/right censored values equal to 0.05
n <- 1000L
p <- 3L
rho <- 0.3
Sigma <- outer(1L:p, 1L:p, function(i, j) rho^abs(i - j))
Z <- rcggm(n = n, Sigma = Sigma, probl = 0.05, probr = 0.05)
out <- cglasso(. ~ ., data = Z)
fitted(out, rho.id = 3L, drop = TRUE)

# Y ~ N(b0 + XB, Sigma)  and probability of left/right censored values equal to 0.05
n <- 1000L
p <- 3L
q <- 2
b0 <- runif(p)
B <- matrix(runif(q * p), nrow = q, ncol = p)
X <- matrix(rnorm(n * q), nrow = n, ncol = q)
rho <- 0.3
Sigma <- outer(1L:p, 1L:p, function(i, j) rho^abs(i - j))
Z <- rcggm(n = n, b0 = b0, X = X, B = B, Sigma = Sigma, probl = 0.05, probr = 0.05)
out <- cglasso(. ~ ., data = Z)
fitted(out, lambda.id = 3L, rho.id = 3L, drop = TRUE)
}

\keyword{models}
\keyword{regression}

