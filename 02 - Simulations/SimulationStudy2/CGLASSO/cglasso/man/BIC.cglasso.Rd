\name{BIC.cglasso}
\alias{BIC}
\alias{BIC.cglasso}

\title{Bayesian Information Criterion}

\description{
\sQuote{\code{BIC}} computes the Bayesian Information Criterion.
}

\usage{
\method{BIC}{cglasso}(object, g = 0, type, mle, ...)
}

\arguments{
\item{object}{an \R object inheriting class \sQuote{\code{cglasso}}, that is, the output of the model-fitting functions \code{\link{cglasso}} and \code{\link{cggm}}.}
\item{g}{a value belonging to the interval \eqn{[0, 1]}{[0, 1]}. Classical BIC is returned by letting \code{g = 0} (default value), whereas extended BIC corresponds to the case \code{g = 0.5}.}
\item{type}{character; if \code{g} is not zero then the measure proposed in Foygel \emph{and other} (2010) is returned by setting \code{type = "FD"}, otherwise (\code{type = "CC"}) returns the measure proposed in Chen \emph{and other} (2008, 2012). See section \sQuote{Details} for more details.}
\item{mle}{logical. TRUE if the measure of goodness-of-fit should be computed using the maximum likelihood estimates. Default depends on the class of the argument \code{object}: \code{mle = FALSE} for objects of class \code{cglasso} and \code{mle = TRUE} for objects of class \code{cggm}.}
\item{\dots}{further arguments passed to the model-fitting function \code{\link{cggm}}.}
}

\details{
\sQuote{\code{BIC}} computes the Bayesian Information Criterion (BIC) for models fitted by \code{\link{cglasso}} or \code{\link{cggm}}. As proposed in Ibrahim \emph{and other} (2008), \code{BIC} computes the measure of goodness-of-fit by replacing the log-likelihood function with the Q-function, that is, the function maximized in the M-Step of the EM-algorithm. The values of the Q-function are computed using \code{\link{QFun}}. By default, for an object of class \code{cglasso} these values are computed using the penalized estimates whereas, if the \code{object} has class \code{cggm}, maximum likelihood estimates are used (see argument \sQuote{\code{mle}} in \code{\link{QFun}}).

By default, \code{BIC} computes the standard BIC measure (\eqn{\gamma = 0}{g = 0}): \deqn{-2\,\mbox{Q-function} + \log(n)\,\mbox{df},}{-2*Q-function + log(n)*df,} where \eqn{n}{n} is the sample size and \eqn{\mbox{df}}{df} represents the number of unique non-zero parameters in the fitted model.

If \eqn{\gamma \ne 0}{g != 0}, the default depends on the number of predictors (\eqn{q}{q}).

If \eqn{q = 0}{q = 0}, \code{BIC} computes the measure of goodness-of-fit proposed in Foygel \emph{and other} (2010) (\code{type = "FD"}): %
\deqn{\mbox{eBIC} = -2\,\mbox{QFun} + (\log n + 4 \, \gamma \, log \, p)\,\mbox{df},}{eBIC = -2*QFun + (log(n) + 4*g*log(p))*df,} %
where \eqn{\gamma}{g} is a value belonging to the interval \eqn{[0, 1]}{(0, 1]} and indexing the measure of goodness-of-fit.

If \eqn{q \ne 0}{q != 0} , \code{BIC} computes the measure of goodness-of-fit proposed in Chen \emph{and other} (2008, 2012) (\code{type = "CC"}): %
\deqn{\mbox{eBIC} = -2\,\mbox{QFun} + (\log n + 2 \, \gamma \, log \, q)\,\mbox{df},}{eBIC = -2*QFun + (log(n) + 2*g*log(q))*df.} %

\code{BIC} can be passed to the functions  \code{\link{select.cglasso}} and \code{\link{summary.cglasso}} to select and print the optimal fitted model, respectively.

The function \code{\link{plot.GoF}} can be used to graphically evaluate the behaviour of the fitted models in terms of goodness-of-fit.
}

\value{
\sQuote{\code{BIC}} returns an \R object of S3 class \dQuote{\code{GoF}}, i.e. a named list containing the following components:
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
Foygel, R. and Drton, M. (2010).
Extended Bayesian Information Criteria for Gaussian Graphical Models.
In: Lafferty, J., Williams, C., Shawe-taylor, J., Zemel, R.s. and Culott, A. (editors),
\emph{Advances in Neural Information Processing Systems 23}. pp. 604--612.

Chen, J. and Chen, Z. (2008) <\doi{10.1093/biomet/asn034}>.
Extended Bayesian information criteria for model selection with large model spaces.
\emph{Biometrika}, Vol. 95(2), pp. 759--771.

Chen, J. and Chen, Z. (2012) <\doi{10.5705/ss.2010.216}>.
Extended BIC for small-n-large-p sparse GLM.
\emph{Statistica Sinica}, Vol. 22, pp. 555--574.

Wit, E., Heuvel, E. V. D., & Romeijn, J. W. (2012) <\doi{10.1111/j.1467-9574.2012.00530.x}>.
All models are wrong...?: an introduction to model uncertainty.
\emph{Statistica Neerlandica}, 66(3), 217-236.
}

\seealso{
\code{\link{cglasso}}, \code{\link{cggm}}, \code{\link{AIC.cglasso}}, \code{\link{QFun}}, \code{\link{plot.GoF}} and \code{\link{summary.cglasso}} 
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
BIC(out)                                    # standard BIC measure
BIC(out, mle = TRUE, g = 0.5, type = "FD")  # eBIC proposed in Foygel and other (2010)

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
BIC(out)                                    # standard BIC measure
BIC(out, mle = TRUE, g = 0.5, type = "FD")  # eBIC proposed in Foygel and other (2010)
BIC(out, mle = TRUE, g = 0.5, type = "CC")  # eBIC proposed in Chen and other (2008, 2010)
}

\keyword{models}
