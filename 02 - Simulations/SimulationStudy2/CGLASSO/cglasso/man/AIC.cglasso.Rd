\name{AIC.cglasso}
\alias{AIC}
\alias{AIC.cglasso}

\title{Akaike Information Criterion}

\description{
\sQuote{\code{AIC}} computes the \sQuote{Akaike Information Criterion}.
}

\usage{
\method{AIC}{cglasso}(object, k = 2, mle, ...)
}

\arguments{
\item{object}{an \R object inheriting class \sQuote{\code{cglasso}}, that is, the output of the model-fitting functions \code{\link{cglasso}} and \code{\link{cggm}}.}
\item{k}{the \emph{penalty} parameter to be used; the default \code{k = 2} is the classical AIC.}
\item{mle}{logical. TRUE if the measure of goodness-of-fit should be computed using the maximum likelihood estimates. Default depends on the class of the argument \code{object}: \code{mle = FALSE} for objects of class \code{cglasso} and \code{mle = TRUE} for objects of class \code{cggm}.}
\item{\ldots}{further arguments passed to \code{\link{cggm}}.}
}

\details{
\sQuote{\code{AIC}} computes the following measure of goodness-of-fit (Ibrahim \emph{and other}, 2008): \deqn{-2\,\mbox{Q-function} + k\,\mbox{df},}{-2*Q-function + k*df,} where \eqn{k}{k} is the \emph{penalty} parameter and \eqn{\mbox{df}}{df} represents the number of unique non-zero parameters in the fitted model.

The values of the Q-function function are computed using \code{\link{QFun}}. By default, for an object of class \code{cglasso} these values are computed using the penalized estimates whereas maximum likelihood estimates are used if the \code{object} is of class \code{cggm} (see argument \sQuote{\code{mle}} in \code{\link{QFun}}).

The Akaike Information Criterion (AIC) is returned by letting \eqn{k = 2} (default value of the function \code{AIC}) whereas the \sQuote{Bayesian Information Criterion} (BIC) is returned by letting \eqn{k = \log(n)}{k = log(n)}, where \eqn{n} is the sample size.

Function \code{AIC} can be passed to the functions \code{\link{select.cglasso}} and \code{\link{summary.cglasso}} to select and print the optimal fitted model, respectively.

The function \code{\link{plot.GoF}} can be used to graphically evaluate the behaviour of the fitted models in terms of goodness-of-fit.
}

\value{
\sQuote{\code{AIC}} return an \R object of S3 class \dQuote{\code{GoF}}, i.e., a named list containing the following components:
\item{value_gof}{a matrix storing the values of the measure of goodness-of-fit used to evaluate the fitted models.}
\item{df}{a matrix storing the number of the estimated non-zero parameters.}
\item{dfB}{a matrix storing the number of estimated non-zero regression coefficients.}
\item{dfTht}{a matrix storing the number of estimated non-zero partial correlation coefficients.}
\item{value}{a matrix storing the values of the Q-function.}
\item{n}{the sample size.}
\item{p}{the number of response variables.}
\item{q}{the number of columns of the design matrix \code{X} used to fit the model.}
\item{lambda}{the \eqn{\lambda}{lambda}-values used to fit the model.}
\item{nlambda}{the number of \eqn{\lambda}{lambda}-values used.}
\item{rho}{the \eqn{\rho}{rho}-values used to fit the model.}
\item{nrho}{the number of \eqn{\rho}{rho}-values used.}
\item{type}{a description of the computed measure of goodness-of-fit.}
\item{model}{a description of the fitted model passed through the argument \code{object}.}
}

\author{Luigi Augugliaro (\email{luigi.augugliaro@unipa.it})}

\references{
Ibrahim, J.G., Zhu, H. and Tang, N. (2008) <\doi{10.1198/016214508000001057}>.
Model selection criteria for missing-data problems using the EM algorithm.
\emph{Journal of the American Statistical Association} \bold{103}, 1648--1658.

Sakamoto, Y., Ishiguro, M., and Kitagawa, G. (1986).
\emph{Akaike Information Criterion Statistics}.
D. Reidel Publishing Company.

Wit, E., Heuvel, E. V. D., & Romeijn, J. W. (2012) <\doi{10.1111/j.1467-9574.2012.00530.x}>.
All models are wrong...?: an introduction to model uncertainty.
\emph{Statistica Neerlandica}, 66(3), 217-236.
}

\seealso{
\code{\link{BIC.cglasso}}, \code{\link{cglasso}}, \code{\link{cggm}}, \code{\link{QFun}}, \code{\link{plot.GoF}} and \code{\link{summary.cglasso}}  
}

\examples{
set.seed(123)

# Y ~ N(b0 + XB, Sigma) and probability of left/right censored values equal to 0.05
n <- 100L
p <- 3L
q <- 2L
b0 <- runif(p)
B <- matrix(runif(q * p), nrow = q, ncol = p)
X <- matrix(rnorm(n * q), nrow = n, ncol = q)
rho <- 0.3
Sigma <- outer(1L:p, 1L:p, function(i, j) rho^abs(i - j))
Z <- rcggm(n = n, b0 = b0, X = X, B = B, Sigma = Sigma, probl = 0.05, probr = 0.05)
out <- cglasso(. ~ ., data = Z)
AIC(out)

out.mle <- cggm(out, lambda.id = 3L, rho.id = 3L)
AIC(out.mle)
}

\keyword{models}
