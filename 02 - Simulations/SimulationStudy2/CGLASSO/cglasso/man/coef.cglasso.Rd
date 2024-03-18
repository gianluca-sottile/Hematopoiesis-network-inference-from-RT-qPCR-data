\name{coef}
\alias{coef}
\alias{coef.cglasso}

\title{Extract Model Coefficients}

\description{
The accessor function \code{coef} extracts model coefficients from an \R object inheriting class \sQuote{\code{cglasso}}.
}

\usage{
\method{coef}{cglasso}(object, type = c("all", "B", "Sigma", "Theta"), lambda.id, rho.id, 
     drop = TRUE, \dots)
}

\arguments{
\item{object}{an \R object inheriting class \sQuote{\code{cglasso}}, that is, the output of the model-fitting functions \code{\link{cglasso}} and \code{\link{cggm}}.}
\item{type}{a description of the desired estimates.}
\item{lambda.id}{an optional vector of integers used to specify the \eqn{\lambda}{lambda}-values.}
\item{rho.id}{an optional vector of integers used to specify the \eqn{\rho}{rho}-values.}
\item{drop}{logical. Dimensions of the required objects can only be dropped if their extent is one.}
\item{\dots}{further arguments passed to or from other methods.}
}

\value{
Coefficients extracted from \sQuote{\code{object}} are returned. By default, a named list storing all the estimated parameters is returned.
}

\author{Luigi Augugliaro (\email{luigi.augugliaro@unipa.it})}

\seealso{
Model-fitting functions \code{\link{cglasso}}, \code{\link{cggm}} and the accessor functions \code{\link{fitted.cglasso}}, \code{\link{residuals.cglasso}}, \code{\link{predict.cglasso}} and \code{\link{impute}}.
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
coef(out, type = "Theta", rho.id = 1:4)
coef(out, type = "Theta", rho.id = 3, drop = TRUE)

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
coef(out, type = "B", lambda.id = 3, rho.id = 1:4)
coef(out, type = "B", lambda.id = 3, rho.id = 3, drop = TRUE)
}

\keyword{regression}
\keyword{models}
