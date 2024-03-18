\name{residuals}
\alias{residuals}
\alias{residuals.cglasso}

\title{Extract Model Residuals}

\description{
Extracts model residuals from an \R object inheriting class \sQuote{\code{cglasso}}.
}

\usage{
\method{residuals}{cglasso}(object, type = c("observed", "working"), lambda.id, rho.id,
          drop = TRUE, \dots)
}

\arguments{
\item{object}{an \R object inheriting class \sQuote{\code{cglasso}}, that is, the output of the model-fitting functions \code{\link{cglasso}} and \code{\link{cggm}}.}
\item{type}{a description of the desired residuals (see section \sQuote{Details}).}
\item{lambda.id}{a vector of integers used to specify the \eqn{\lambda}{lambda}-values.}
\item{rho.id}{a vector of integers used to specify the \eqn{\rho}{rho}-values.}
\item{drop}{logical. Dimensions can only be dropped if their extent is one.}
\item{\dots}{further arguments passed to or from other methods.}
}

\details{
The accessor function \sQuote{\code{residuals}} returns an array storing the \sQuote{observed} or \sQuote{working} residuals, depending on the argument \code{type}.

The \sQuote{observed} residuals are defined as the difference between the observed response values and the corresponding fitted expected values. For missing and censored response values, the \sQuote{observed} residuals are set equal to \sQuote{\code{NA}}.

The \sQuote{working} residuals are obtained as a byproduct of the EM-algorithm used to fit the model. In the E-step, the algorithm computes the \sQuote{working} response matrix, that is, the response matrix with missing values replaced by the conditional expected values of the multivariate normal distribution and censored values replaced by the conditional expected values of the multivariate truncated normal distribution. The \sQuote{working} residuals are defined as the difference between \sQuote{working} responses and fitted expected values. 
}

\value{
The accessor function \sQuote{\code{residuals}} returns an array storing the desired residuals.
}

\author{Luigi Augugliaro (\email{luigi.augugliaro@unipa.it})}

\seealso{
Model-fitting functions \code{\link{cglasso}}, \code{\link{cggm}} and the accessor functions \code{\link{coef.cglasso}}, \code{\link{fitted.cglasso}}, \code{\link{predict.cglasso}} and \code{\link{impute}}.

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
residuals(out, type = "observed", rho.id = 3L, drop = TRUE)
residuals(out, type = "working", rho.id = 3L, drop = TRUE)

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
residuals(out, type = "observed", lambda.id = 3L, rho.id = 3L, drop = TRUE)
residuals(out, type = "working", lambda.id = 3L, rho.id = 3L, drop = TRUE)
}

\keyword{models}
\keyword{regression}

