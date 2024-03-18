\name{predict}
\alias{predict}
\alias{predict.cglasso}
\alias{predict.cggm}

\title{Predict Method for cglasso and cggm Fits}

\description{
Obtains predictions from an \R object inheriting class \sQuote{\code{cglasso}}.
}

\usage{
\method{predict}{cglasso}(object, type = c("B", "mu", "Sigma", "Theta"), X.new, lambda.new, rho.new,
        \dots)
        
\method{predict}{cggm}(object, X.new, \dots)
}

\arguments{
\item{object}{an \R object inheriting class \sQuote{\code{cglasso}}, that is, the output of the model-fitting function \sQuote{\code{cglasso}} or \sQuote{\code{cggm}}.}
\item{type}{a description of prediction required.}
\item{X.new}{matrix of new values for \code{X} at which predictions are to be made. This argument is used only if \sQuote{\code{type = "mu"}}.}
\item{lambda.new}{value of the tuning parameter \eqn{\lambda}{lambda} at which predictions are required.}
\item{rho.new}{value of the tuning parameter \eqn{\rho}{rho} at which predictions are required.}
\item{\dots}{further arguments passed to or from other methods.}
}

\details{
If \code{object} has S3 class \sQuote{\code{cglasso}}, then for a new pair of the tuning parameters \eqn{\lambda}{lambda} and \eqn{\rho}{rho}, the \code{predict} function can be used to predict the estimate of the regression coefficient matrix (\sQuote{\code{type = "B"}}), the estimate of the covariance matrix (\sQuote{\code{type = "Sigma"}}) or the estimate of the precision matrix (\sQuote{\code{type = "Theta"}}). If \code{X.new} is missing and \sQuote{\code{type = "mu"}}, then the \code{predict} function returns the predicted values using the matrix of predictors \code{X}, otherwise the predicted fitted values are computed using the matrix \code{X.new}.

For a new pair of the tuning parameters \eqn{\lambda}{lambda} and \eqn{\rho}{rho}, the predicted values are computed using a bilinear interpolation.

If the \code{object} has S3 class \sQuote{\code{cggm}}, then the \code{predict} function returns only the predicted fitted values using the argument \code{X.new}. 
}

\value{
The matrix of predicted values.
}

\author{Luigi Augugliaro (\email{luigi.augugliaro@unipa.it})}

\seealso{
Model-fitting function \code{\link{cglasso}} and the other accessor functions \code{\link{coef.cglasso}}, \code{\link{fitted.cglasso}}, \code{\link{residuals.cglasso}} and \code{\link{impute}}.
}

\examples{
set.seed(123)
# Y ~ N(0, Sigma) and probability of left/right censored values equal to 0.05
n <- 100L
p <- 3L
rho <- 0.3
Sigma <- outer(1L:p, 1L:p, function(i, j) rho^abs(i - j))
Z <- rcggm(n = n, Sigma = Sigma, probl = 0.05, probr = 0.05)

out <- cglasso(. ~ ., data = Z)
rho.new <- mean(out$rho)
Theta.pred <- predict(out, rho.new = rho.new, type = "Theta")
Theta.pred

# Y ~ N(b0 + XB, Sigma)  and probability of left/right censored values equal to 0.05
n <- 100L
p <- 3L
q <- 2
b0 <- runif(p)
B <- matrix(runif(q * p), nrow = q, ncol = p)
X <- matrix(rnorm(n * q), nrow = n, ncol = q)
rho <- 0.3
Sigma <- outer(1L:p, 1L:p, function(i, j) rho^abs(i - j))
Z <- rcggm(n = n, b0 = b0, X = X, B = B, Sigma = Sigma, probl = 0.05, probr = 0.05)

out <- cglasso(. ~ ., data = Z)
rho.new <- mean(out$rho)
lambda.new <- mean(out$lambda)
Theta.pred <- predict(out, lambda.new = lambda.new, rho.new = rho.new, type = "Theta")
Theta.pred
}

\keyword{models}
\keyword{regression}
