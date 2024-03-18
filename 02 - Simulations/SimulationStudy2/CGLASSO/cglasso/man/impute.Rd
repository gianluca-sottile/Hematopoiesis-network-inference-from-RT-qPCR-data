\name{impute}
\alias{impute}

\title{Imputation of Missing and Censored Values}

\description{
Imputes multivariate missing and censored values.
}

\usage{
impute(object, type = c("mar", "censored", "both"), lambda.new, rho.new)
}

\arguments{
\item{object}{an \R object inheriting class \sQuote{\code{cglasso}}, that is, the output of the model-fitting functions \code{\link{cglasso}} and \code{\link{cggm}}.}
\item{type}{a description of the imputation required (see section \sQuote{Description}). By default only the missing-at-random values are imputed.}
\item{lambda.new}{value of the tuning parameter \eqn{\lambda}{lambda} at which the imputations are required.}
\item{rho.new}{value of the tuning parameter \eqn{\rho}{rho} at which the imputations are required.}
}

\details{
The \code{impute} function returns the response matrix with the imputed missing-at-random values (\sQuote{\code{type = "mar"}}), the imputed censored values (\sQuote{\code{type = "censored"}}), or both (\sQuote{\code{type = "both"}}).

If \sQuote{\code{type = "mar"}} then, for each row, the missing response values are replaced with the expected values of a multivariate normal distribution conditioned on the observed response values.

If \sQuote{\code{type = "censored"}} then, for each row, the censored response values are imputed using the expected values of a multivariate truncated normal distribution conditioned on the observed response values.

If \sQuote{\code{type = "both"}} then missing-at-random and censored values are imputed. In this case the returned matrix corresponds to the working response matrix computed during the E-Step.
}

\value{
The \code{impute} function returns a response matrix with the required imputed data.
}

\author{Luigi Augugliaro (\email{luigi.augugliaro@unipa.it})}

\seealso{
Model-fitting functions \code{\link{cglasso}}, \code{\link{cggm}} and the accessor functions \code{\link{coef.cglasso}}, \code{\link{fitted.cglasso}}, \code{\link{residuals.cglasso}} and \code{\link{predict.cglasso}}.
}

\examples{
set.seed(123)
# Y ~ N(0, Sigma) and 
# 1. probability of left/right censored values equal to 0.05
# 2. probability of missing-at-random values equal to 0.05
n <- 100L
p <- 3L
rho <- 0.3
Sigma <- outer(1L:p, 1L:p, function(i, j) rho^abs(i - j))
Z <- rcggm(n = n, Sigma = Sigma, probl = 0.05, probr = 0.05, probna = 0.05)

out <- cglasso(. ~ ., data = Z)
rho.new <- mean(out$rho)

# imputing missing values
Y.impute <- impute(out, type = "mar", rho.new = rho.new)

# imputing censored values
Y.impute <- impute(out, type = "censored", rho.new = rho.new)

# imputing missing and censored values
Y.impute <- impute(out, type = "both", rho.new = rho.new)


# Y ~ N(XB, Sigma) and 
# 1. probability of left/right censored values equal to 0.05
# 2. probability of missing-at-random values equal to 0.05
n <- 100L
p <- 3L
q <- 2
b0 <- runif(p)
B <- matrix(runif(q * p), nrow = q, ncol = p)
X <- matrix(rnorm(n * q), nrow = n, ncol = q)
rho <- 0.3
Sigma <- outer(1L:p, 1L:p, function(i, j) rho^abs(i - j))
Z <- rcggm(n = n, b0 = b0, X = X, B = B, Sigma = Sigma, probl = 0.05, probr = 0.05, 
           probna = 0.05)

out <- cglasso(. ~ ., data = Z)
lambda.new <- mean(out$lambda)
rho.new <- mean(out$rho)

# imputing missing values
Y.impute <- impute(out, type = "mar", lambda.new = lambda.new, rho.new = rho.new)

# imputing censored values
Y.impute <- impute(out, type = "censored", lambda.new = lambda.new, rho.new = rho.new)

# imputing missing and censored values
Y.impute <- impute(out, type = "both", lambda.new = lambda.new, rho.new = rho.new)
}

\keyword{models}
\keyword{regression}
