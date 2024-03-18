\name{cglasso}
\alias{cglasso}

\title{Conditional Graphical Lasso Estimator}

\description{
\sQuote{\code{cglasso}} fits the conditional graphical lasso model to datasets with censored and/or missing values.
}

\usage{
cglasso(formula, data, subset, contrasts = NULL, diagonal = FALSE,
        weights.B = NULL, weights.Tht = NULL, nlambda, lambda.min.ratio,
        lambda, nrho, rho.min.ratio, rho, maxit.em = 1.0E+4, thr.em = 1.0E-3,
        maxit.bcd = 1.0E+5, thr.bcd = 1.0E-4, trace = 0L,
        algorithm = c("glasso", "admm"))
}

\arguments{
\item{formula}{an object of class \sQuote{\code{formula}} (or one that can be coerced to that class): a symbolic description of the model to be fitted.}
\item{data}{an \R object of S3 class \sQuote{\code{datacggm}}, that is, the output of the function \code{\link{datacggm}}. See section \sQuote{Description} for more details.}
\item{subset}{an optional vector specifying a subset of observations to be used in the fitting process.}
\item{contrasts}{an optional list. See the \code{contrasts.arg} of \code{\link[stats]{model.matrix.default}}.}
\item{diagonal}{logical. Should diagonal entries of the concentration matrix be penalized? Default is \sQuote{\code{diagonal = FALSE}}.}
\item{weights.B}{an optional \eqn{q\times p}{q x p} dimensional matrix of non-negative weights used to penalize the regression coefficients (intercepts are unpenalized). This matrix can be also used to specify the unpenalized regression coefficients (\sQuote{\code{weights.B[i, j] = 0}}) or the structural zeros (\sQuote{\code{weights.B[i, j] = +Inf}}). By default, all weights are set equal to 1, meaning that the penalized regression coefficients are unweighted.}
\item{weights.Tht}{an optional symmetric matrix of non-negative weights used to penalize the partial regression coefficients. This matrix can be used to specify the unpenalized partial correlation coefficients (\sQuote{\code{weights.Tht[i, j] = 0}}) or the structural zeros in the precision matrix (\sQuote{\code{weights.Tht[i, j] = +Inf}}). By default, the off-diagonal entries of the matrix \sQuote{\code{weights.Tht}}  are set equal to 1, meaning that the partial correlation coefficients are unweighted, whereas the diagonal entries are equal to 0, that is the diagonal entries of the precision matrix are unpenalized.}
\item{nlambda}{integer. The number of \eqn{\lambda}{lambda}-values used to penalize the regression coefficients. By default \sQuote{\code{nlambda = 10}}.}
\item{lambda.min.ratio}{the smallest \eqn{\lambda}{lambda}-value is defined as a fraction of \sQuote{\code{lambda.max}} (i.e., the smallest \eqn{\lambda}{lambda}-value for which all the estimated regression coefficients are equal to zero). The default depends on the sample size \sQuote{\eqn{n}{n}} relative to the number of predictors \sQuote{\eqn{q}{q}}. If \sQuote{\eqn{q < n}{q < n}}, default is \sQuote{1.0E-6} otherwise the value \sQuote{1.0E-2} is used as default. A very small value of \sQuote{\code{lambda.min.ratio}} will lead to a saturated fitted model in the \sQuote{\eqn{q < n}{q < n}} case.}
\item{lambda}{an optional user-supplied decreasing sequence of \eqn{\lambda}{lambda}-values. By default \code{cglasso} computes a sequence of \eqn{\lambda}{lambda}-values using \code{nlambda} and \code{lambda.min.ratio}. If \code{lambda} is supplied then \code{nlambda} and \code{lambda.min.ratio} are overwritten. WARNING: use with care and avoid supplying a single \eqn{\lambda}{lambda}-value; supply instead a decreasing sequence.}
\item{nrho}{integer. The number of \eqn{\rho}{rho}-values used to penalize the partial correlation coefficients. By default \sQuote{\code{nrho = 10}}.}
\item{rho.min.ratio}{the smallest \eqn{\rho}{rho}-value is defined as a fraction of \sQuote{\code{rho.max}} (i.e., the smallest \eqn{\rho}{rho}-value for which all the estimated partial correlation coefficients are equal to zero). The default depends on the sample size \sQuote{\eqn{n}{n}} relative to the number of response variables \sQuote{\eqn{p}{p}}. If \sQuote{\eqn{p < n}{p < n}}, the default is \sQuote{1.0E-6} otherwise the value \sQuote{1.0E-2} is used as default. A very small value of \sQuote{\code{rho.min.ratio}} will lead to a saturated fitted model in the \sQuote{\eqn{p < n}{p < n}} case.}
\item{rho}{an optional user supplied decreasing sequence of \eqn{\rho}{rho}-values. By default \code{cglasso} computes a sequence of \eqn{\rho}{rho}-values using \code{nrho} and \code{rho.min.ratio}. If \code{rho} is supplied then \code{nrho} and \code{rho.min.ratio} are overwritten. WARNING: use with care and avoid supplying a single \eqn{\rho}{rho}-value; supply instead a decreasing sequence.}
\item{maxit.em}{maximum number of iterations of the EM algorithm. Default is \code{1.0E+4}.}
\item{thr.em}{threshold for the convergence of the EM algorithm. Default value is \code{1.0E-4}.}
\item{maxit.bcd}{maximum number of iterations of the glasso algorithm. Default is \code{1.0E+5}.}
\item{thr.bcd}{threshold for the convergence of the glasso algorithm. Default is \code{1.0E-4}.}
\item{trace}{integer for printing information out as iterations proceed: \code{trace = 0} no information is printed out; \code{trace = 1} minimal information is printed on screen; \code{trace = 2} detailed information is printed on screen.}
\item{algorithm}{the type of algorithm used to estimate the graphical lasso model.}
}
\details{
\code{cglasso} is the main model-fitting function and can be used to fit a broad range of extensions of the glasso estimator (Friedman \emph{and other}, 2008). It is specifically proposed to study datasets with censored and/or missing response values. To help the user, the \code{cglasso} function has been designed to automatically select the most suitable extension by using the information stored in the \sQuote{\code{\link{datacggm}}} object passed through \code{Z}. 

Below we sum up the available extenions:
\itemize{
\item if only left/right-censored are observed in the response matrix \code{Y} (without missing values) and no predictor matrix is stored in \code{Z}, then \code{cglasso} computes the censored glasso estimator (Augugliaro \emph{and other}, 2020a);
\item if only left/right-censored are observed in the response matrix \code{Y} (without missing values) and and a predictor matrix \code{X} is stored in \code{Z}, then \code{cglasso} computes the conditional censored glasso estimator (Augugliaro \emph{and other}, 2020b);
\item if only missing values are stored in the response matrix \code{Y} (without censored values), then \code{cglasso} computes the missglasso estimator (Stadler \emph{and other}, 2012);
\item starting with version 2.0.0, \code{cglasso} can also handle datasets with both missing and censored response values.
}
See section \sQuote{Examples} for some example.

The model-fitting function \code{cglasso} returns an \R object of S3 class \sQuote{\code{cglasso}} for which there are available a set of accessor functions, a set of functions designed to evaluate the goodness-of-fit of the fitted models and, finally, a set of functions developed to analyze the selected network. The function \code{\link{ShowStructure}} can be used to show the structure of the package.

The accessor functions \code{\link{coef.cglasso}}, \code{\link{fitted.cglasso}}, \code{\link{residuals.cglasso}}, \code{\link{predict.cglasso}} and \code{\link{impute}} can be used to extract various useful features of the object fitted by \code{cglasso}.

For an \R object returned by \code{cglasso}, the functions \code{\link{AIC.cglasso}} and \code{\link{BIC.cglasso}} can be used to evaluate the goodness-of-fit of the fitted models. Usually, these functions are used together with the function \code{\link{summary.cglasso}}, which gives more information about the sequence of fitted models. The plotting function \code{\link{plot.GoF}} can be used to graphically identify the optimal pair of the tuning parameters.

Given a pair of the tuning paremeters, the functions \code{\link{to_graph}} and \code{\link{plot.cglasso2igraph}} can be used to analyze and show the selected network. Finally, the function \code{\link{cggm}} can be used to produce post-hoc maximum likelihood refitting of the selected graphical model.  
}

\value{
\code{cglasso} returns an object of S3 class \dQuote{\code{cglasso}}, i.e., a named list containing the
following components:
\item{call}{the call that produced this object.}
\item{Yipt}{an array of dimension \sQuote{\code{n x p x nlambda x nrho}} storing the \sQuote{working response matrices}, that is, \code{Yipt[, , i, j]} is used as response matrix to compute the multilasso estimator (Augugliaro \emph{and other}, 2020b). The accessor function \sQuote{\code{\link{impute}}} can be used to extract the desidered imputed matrix.}
\item{B}{an array of dimension \sQuote{\code{(q + 1) x p x nlambda x nrho}} storing the penalized estimate of the regression coefficient matrices. The accessor function \sQuote{\code{\link{coef.cglasso}}} can be used to extract the desired estimates.}
\item{mu}{an array of dimension \sQuote{\code{n x p x nlambda x nrho}} storing the fitted values. The accessor function \sQuote{\code{\link{fitted.cglasso}}} can be used to extract the desired fitted values.}
\item{R}{an array of dimension \sQuote{\code{n x p x nlambda x nrho}} storing the \sQuote{working residuals}, that is, \code{R[, , i, j]} is defined as the difference between \code{Yipt[, , i, j]} and \code{mu[, , i, j]}. The accessor function \sQuote{\code{\link{residuals.cglasso}}} can be used to extract the desidered residual matrix.}
\item{S}{an array of dimension \sQuote{\code{p x p x nlambda x nrho}} storing the \sQuote{working empirical covariance matrices}, that is, \sQuote{\code{S[, , i, j]}} is used to compute the glasso estimator.}
\item{Sgm}{an array of dimension \sQuote{\code{p x p x nlambda x nrho}} storing the estimated covariance matrices. The accessor function  \sQuote{\code{\link{coef}}} can be used to extract the desired estimates.}
\item{Tht}{an array of dimension \sQuote{\code{p x p x nlambda x nrho}} storing the estimated precision matrices. The accessor funtion \sQuote{\code{\link{coef}}} can be used to extract the desired estimates.}
\item{dfB}{a matrix of dimension \sQuote{\code{(p + 1) x nlambda x nrho}} storing the number of estimated non-zero regression coefficients. Only for internal purpose.}
\item{dfTht}{a matrix of dimension \sQuote{\code{nlambda x nrho}} storing the number of estimated non-zero (off-diagonal) partial correlation coefficients. Only for internal purpose.}
\item{InfoStructure}{a named list whose elements contain information about the estimated networks. Only for internal purpose.}
\item{nit}{an array of dimension \sQuote{\code{2 x nlambda x nrho}} storing the number of EM steps.}
\item{Z}{the \sQuote{\code{datacggm}} object used to compute the censored graphical lasso estimator.}
\item{diagonal}{the flag used to specify if the diagonal entries of the precision matrix are penalized.}
\item{weights.B}{the matrix of non-negative weights used for the regression coefficients.}
\item{weights.Tht}{the matrix of non-negative weights used for the precision matrix.}
\item{nlambda}{the number of \eqn{\lambda}{\lambda}-values used.}
\item{lambda.min.ratio}{the value used to compute the smallest \eqn{\lambda}{lambda}-value.}
\item{lambda}{the sequence of \eqn{\lambda}{lambda}-values used to fit the model.}
\item{nrho}{the number of \eqn{\rho}{\rho}-values used.}
\item{rho.min.ratio}{the value used to compute the smallest \eqn{\rho}{rho}-value.}
\item{rho}{the sequence of \eqn{\rho}{rho}-values used to fit the model.}
\item{model}{a description of the fitted model.}
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

\references{
Augugliaro, L., Abbruzzo, A., and Vinciotti, V. (2020a) <\doi{10.1093/biostatistics/kxy043}>.
\eqn{\ell_1}{l1}-Penalized censored Gaussian graphical model.
\emph{Biostatistics} \bold{21}, e1--e16.

Augugliaro, L., Sottile, G., and Vinciotti, V. (2020b) <\doi{10.1007/s11222-020-09945-7}>.
The conditional censored graphical lasso estimator.
\emph{Statistics and Computing} \bold{30}, 1273--1289.

Friedman, J.H., Hastie, T., and Tibshirani, R. (2008) <\doi{10.1093/biostatistics/kxm045}>.
Sparse inverse covariance estimation with the graphical lasso.
\emph{Biostatistics} \bold{9}, 432--441.

Stadler, N. and Buhlmann, P. (2012) <\doi{10.1007/s11222-010-9219-7}>.
Missing values: sparse inverse covariance estimation and an extension to sparse regression.
\emph{Statistics and Computing} \bold{22}, 219--235.
}

\seealso{
\code{\link{datacggm}}, \code{\link{coef.cglasso}}, \code{\link{fitted.cglasso}}, \code{\link{residuals.cglasso}}, \code{\link{predict.cglasso}}, \code{\link{impute}}, \code{\link{AIC.cglasso}}, \code{\link{BIC.cglasso}}, \code{\link{summary.cglasso}}, \code{\link{select.cglasso}}, \code{\link{plot.GoF}}, \code{\link{to_graph}}, \code{\link{plot.cglasso2igraph}}, \code{\link{cggm}} and \code{\link{ShowStructure}}.
}

\examples{
set.seed(123)

# Model 1: censored glasso estimator (Augugliaro \emph{and other}, 2020a)
# Y ~ N(0, Sigma) and probability of left/right censored values equal to 0.05
n <- 1000L
p <- 3L
rho <- 0.3
Sigma <- outer(1L:p, 1L:p, function(i, j) rho^abs(i - j))
Z <- rcggm(n = n, Sigma = Sigma, probl = 0.05, probr = 0.05)
out <- cglasso(. ~ ., data = Z)
out

# Model 2: conditional censored glasso estimator (Augugliaro \emph{and other}, 2020b)
# Y ~ N(b0 + XB, Sigma)  and probability of left/right censored values equal to 0.05
n <- 1000L
p <- 3L
q <- 2L
b0 <- runif(p)
B <- matrix(runif(q * p), nrow = q, ncol = p)
X <- matrix(rnorm(n * q), nrow = n, ncol = q)
rho <- 0.3
Sigma <- outer(1L:p, 1L:p, function(i, j) rho^abs(i - j))
Z <- rcggm(n = n, b0 = b0, X = X, B = B, Sigma = Sigma, probl = 0.05, probr = 0.05)
out <- cglasso(. ~ ., data = Z)
out

# Model 3: missglasso estimator (Stadler \emph{and other}, 2012)
# Y ~ N(b0 + XB, Sigma)  and probability of missing-at-random values equal to 0.05
n <- 1000L
p <- 3L
q <- 2L
b0 <- runif(p)
B <- matrix(runif(q * p), nrow = q, ncol = p)
X <- matrix(rnorm(n * q), nrow = n, ncol = q)
rho <- 0.3
Sigma <- outer(1L:p, 1L:p, function(i, j) rho^abs(i - j))
Z <- rcggm(n = n, b0 = b0, X = X, B = B, Sigma = Sigma, probna = 0.05)
out <- cglasso(. ~ ., data = Z)
out

# Model 4: mixed estimator
# Y ~ N(b0 + XB, Sigma)  and
# 1. probability of left/right censored values equal to 0.05
# 2. probability of missing-at-random values equal to 0.05
n <- 1000L
p <- 3L
q <- 2L
b0 <- runif(p)
B <- matrix(runif(q * p), nrow = q, ncol = p)
X <- matrix(rnorm(n * q), nrow = n, ncol = q)
rho <- 0.3
Sigma <- outer(1L:p, 1L:p, function(i, j) rho^abs(i - j))
Z <- rcggm(n = n, X = X, b0 = b0, B = B, Sigma = Sigma, probl = 0.05, probr = 0.05,
           probna = 0.05)
out <- cglasso(. ~ ., data = Z)
out
}

\keyword{regression}
