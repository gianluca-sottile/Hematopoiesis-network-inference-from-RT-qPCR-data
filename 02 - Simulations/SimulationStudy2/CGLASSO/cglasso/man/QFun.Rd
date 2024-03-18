\name{QFun}
\alias{QFun}

\title{Extract Q-Function}

\description{
\sQuote{\code{QFun}} extracts the values of the Q-function from an \R object inheriting class \sQuote{\code{cglasso}}.
}

\usage{
QFun(object, mle, verbose = FALSE, ...)
}

\arguments{
\item{object}{an \R object inheriting class \sQuote{\code{cglasso}}, that is, the output of the model-fitting functions \code{\link{cglasso}} and \code{\link{cggm}}.}
\item{mle}{logical. Should Q-values be computed using the maximum likelihood estimates? Default depends on the class of the argument \code{object}: \code{mle = FALSE} for objects with class \code{cglasso} and \code{mle = TRUE} for objects with class \code{cggm}.}
\item{verbose}{logical for printing out a progress bar on the \R console. Default is \code{verbose = FALSE}.}
\item{\ldots}{further arguments passed to \code{\link{cggm}}.}
}

\details{
\sQuote{\code{QFun}} returns the value of the Q-function, i.e., the value of the function maximised in the M-step of the EM algorithm. The Q-function is defined as follows: \deqn{\frac{n}{2}\left\{\log det\Theta - tr(S\Theta) - p\log(2\pi)\right\},}{n/2 \{log det Tht - tr(S Tht) - p log(2 pi)\},} where \eqn{S}{S} is the \sQuote{working} empirical covariance matrix computed during the E-step.

\code{QFun} is used as a workhorse function to compute the measures of goodness-of-fit returned by the functions \code{\link{AIC.cglasso}} and \code{\link{BIC.cglasso}}.

The function \sQuote{\code{print.QFun}} is used the improve the readability of the results.
}

\value{
\sQuote{\code{QFun}} returns an \R object of S3 class \dQuote{\code{QFun}}, i.e., a named list containing the following components:
\item{value}{a matrix with the values of the Q-function.}
\item{df}{a matrix with the number of estimated non-zero parameters.}
\item{dfB}{a matrix with the number of estimated non-zero regression coefficients.}
\item{dfTht}{a matrix with the number of estimated non-zero partial correlation coefficients.}
\item{n}{the sample size.}
\item{p}{the number of response variables.}
\item{q}{the number of columns of the design matrix \code{X} used to fit the model.}
\item{lambda}{the \eqn{\lambda}{lambda}-values used to fit the model.}
\item{nlambda}{the number of \eqn{\lambda}{lambda}-values.}
\item{rho}{the \eqn{\rho}{rho}-values used to fit the model.}
\item{nrho}{the number of \eqn{\rho}{rho}-values.}
\item{model}{a description of the fitted model passed through the argument \code{object}.}
}

\author{Luigi Augugliaro (\email{luigi.augugliaro@unipa.it})}

\seealso{
\code{\link{AIC.cglasso}}, \code{\link{BIC.cglasso}}, \code{\link{cglasso}}, \code{\link{cggm}}, \code{\link{summary.cglasso}}, \code{\link{select.cglasso}} and \code{\link{to_graph}}.
}

\examples{
set.seed(123)
# Y ~ N(b0+ XB, Sigma)  and
# 1. probability of left/right censored values equal to 0.05
# 2. probability of missing-at-random values equal to 0.05
n <- 100L
p <- 3L
q <- 2L
b0 <- runif(p)
B <- matrix(runif(q * p), nrow = q, ncol = p)
X <- matrix(rnorm(n * q), nrow = n, ncol = q)
rho <- 0.3
Sigma <- outer(1:p, 1:p, function(i, j) rho^abs(i - j))
Z <- rcggm(n = n, b0 = b0, X = X, B = B, Sigma = Sigma, probl = 0.05, probr = 0.5, 
           probna = 0.05)
out <- cglasso(. ~ ., data = Z)
QFun(out)

out.mle <- cggm(out, lambda.id = 3L, rho.id = 3L)
QFun(out.mle)
}

\keyword{models}

