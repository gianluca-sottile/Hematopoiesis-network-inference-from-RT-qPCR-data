\name{cggm}
\alias{cggm}

\title{Post-Hoc Maximum Likelihood Refitting of a Conditional Graphical Lasso}

\description{
\sQuote{\code{cggm}} is used to perform post-hoc maximum likelihood refitting of a selected conditional graphical lasso model with censored and/or missing values.
}

\usage{
cggm(object, GoF = AIC, lambda.id, rho.id, tp.min = 1.0E-6, ntp = 100L,
     maxit.em = 1.0E+4, thr.em = 1.0E-3, maxit.bcd = 1.0E+5, thr.bcd = 1.0E-4, 
     trace = 0L, \ldots)
}

\arguments{
\item{object}{an \R object of S3 class \sQuote{\code{cglasso}}, that is, the output of the function \code{\link{cglasso}}.}
\item{GoF}{a valid goodness-of-fit function, such as \code{\link{AIC.cglasso}} or \code{\link{BIC.cglasso}}, or an \R object of class \sQuote{\code{GoF}}.}
\item{lambda.id}{an optional integer used to identify a \eqn{\lambda}{lambda}-value stored in \code{object}. See section \sQuote{Details} for more details.}
\item{rho.id}{an optional integer used to identify a \eqn{\rho}{rho}-value stored in \code{object}. See section \sQuote{Details} for more details.}
\item{tp.min}{the smallest \eqn{\lambda}{lambda} and \eqn{\rho}{rho} value. See section \sQuote{Details} for more details.}
\item{ntp}{integer; number of \eqn{\lambda}{lambda}- and \eqn{\rho}{rho}-values used to compute the coefficients path. See section \sQuote{Details} for more details.}
\item{maxit.em}{maximum number of iterations of the EM algorithm. Default is \code{1.0E+4}.}
\item{thr.em}{threshold for the convergence of the EM algorithm. Default value is \code{1.0E-4}.}
\item{maxit.bcd}{maximum number of iterations of the glasso algorithm. Default is \code{1.0E+5}.}
\item{thr.bcd}{threshold for the convergence of the glasso algorithm. Default is \code{1.0E-4}.}
\item{trace}{integer for printing information out as iterations proceed: \code{trace = 0} no information is printed out on screen; \code{trace = 1} minimal information is printed; \code{trace = 2} detailed information is printed on screen.}
\item{\dots}{the \emph{penalty} parameter passed to the goodness-of-fit function (argument \sQuote{\code{GoF}}).}
}

\details{
The model-fitting function \code{cggm} is used to obtain the maximum likelihood estimates of a censored Gaussian graphical model whose structure was found by penalised cglasso inference. That is, given a fitted cglasso model, identified along the path via a goodness-of-fit function (passed by \code{GoF}) or via the optional arguments \code{lambda.id} and \code{rho.id}), \code{cggm} computes the maximum likelihood estimates of the parameters with zero constraints associated to the given structure.

Maximum likelihood estimates are computed using \code{cglasso} as workhorse function, that is, \code{cggm} fits a sequence of cglasso models, with length equal to \code{ntp}, reducing \eqn{\lambda}{lambda} and \eqn{\rho}{rho} until they are equal to the value \sQuote{\code{min(tp.min, lambda.min, rho.min)}}, where \code{lambda.min} and \code{rho.min} are the smallest \eqn{\lambda}{lambda}- and \eqn{\rho}{rho}-values stored in \code{object}. Maximum likelihood estimates are obtained from the last fitted cglasso model.

The model-fitting function \code{cggm} returns an \R object inheriting the S3 class \sQuote{\code{cglasso}}, for which all printing and plotting functions, designed for a \code{cglasso} object, can be used (see section \sQuote{\code{See Also}}). Function \code{\link{ShowStructure}} can be used to show the structure of the package.

NOTE: if the user tries a problem when the sample size is not large enough, then increase the number of fitted cglasso models (argument \code{ntp}).
}

\value{
\code{cggm} returns an object of S3 class \dQuote{\code{cggm}}, i.e., a list containing the
following components:
\item{call}{the call that produced this object.}
\item{Yipt}{an array storing the \sQuote{working response matrix}.}
\item{B}{the maximum likelihood estimate of the regression coefficient matrix.}
\item{mu}{the fitted expected values.}
\item{R}{the \sQuote{working residuals} matrix.}
\item{S}{the \sQuote{working empirical covariance matrix}.}
\item{Sgm}{the maximum likelihood estimate of the covariance matrix.}
\item{Tht}{the maximum likelihood estimate of the precision matrix.}
\item{dfB}{the number of estimated non-zero regression coefficients. Only for internal purpose.}
\item{dfTht}{the number of estimated non-zero (off-diagonal) partial correlation coefficients. Only for internal purpose.}
\item{InfoStructure}{a named list whose elements are used to store the information about the estimated networks. Only for internal purpose.}
\item{nit}{the number of EM steps.}
\item{Z}{the \sQuote{\code{datacggm}} object used to compute the censored graphical lasso estimator.}
\item{nlambda}{Only for internal purpose.}
\item{lambda}{the \eqn{\lambda}{lambda}-value of the selected cglasso model.}
\item{nrho}{Only for internal purpose.}
\item{rho}{the \eqn{\rho}{rho}-value of the selected cglasso model.}
\item{maxit.em}{maximum number of iterations of the EM algorithm.}
\item{thr.em}{threshold for the convergence of the EM algorithm.}
\item{maxit.bcd}{maximum number of iterations of the glasso algorithm.}
\item{thr.bcd}{threshold for the convergence of the glasso algorithm.}
\item{conv}{a description of the error that has occurred.}
\item{subrout}{the name of the Fortran subroutine where the error has occurred (for internal debug only).}
\item{trace}{the integer used for printing information on screen.}
\item{nobs}{the sample size}
\item{nresp}{the number of response variables used to fit the model.}
\item{npred}{the number of predictors used to fit the model.}
}

\author{Luigi Augugliaro (\email{luigi.augugliaro@unipa.it})}

\seealso{
\code{\link{cglasso}}, \code{\link{coef.cglasso}}, \code{\link{fitted.cglasso}}, \code{\link{residuals.cglasso}}, \code{\link{predict.cggm}}, \code{\link{impute}}, \code{\link{AIC.cglasso}}, \code{\link{BIC.cglasso}}, \code{\link{summary.cglasso}}, \code{\link{to_graph}}, \code{\link{plot.cglasso2igraph}} and \code{\link{ShowStructure}}.
}

\examples{
set.seed(123)
# Y ~ N(XB, Sigma)  and
# 1. probability of left/right censored values equal to 0.05
# 2. probability of missing-at-random values equal to 0.05
n <- 100L
p <- 3L
q <- 2L
b0 <- runif(p)
B <- matrix(runif(q * p), nrow = q, ncol = p)
X <- matrix(rnorm(n * q), nrow = n, ncol = q)
rho <- 0.3
Sigma <- outer(1L:p, 1L:p, function(i, j) rho^abs(i - j))
Z <- rcggm(n = n, b0 = b0, X = X, B = B, Sigma = Sigma, probl = 0.05, probr = 0.5, probna = 0.05)
out <- cglasso(. ~ ., data = Z)

# MLE of the censored Gaussian graphical model identified by 'BIC'
out.mle <- cggm(out, GoF = BIC)
out.mle

# accessor functions
coef(out.mle, drop = TRUE)
fitted(out.mle, drop = TRUE)
residuals(out.mle, type = "working", drop = TRUE)
impute(out.mle, type = "both")

# goodness-of-fit functions
AIC(out.mle)
BIC(out.mle)
summary(out.mle)

# network analysis
out.graph <- plot(out.mle)
out.graph
}

\keyword{regression}
