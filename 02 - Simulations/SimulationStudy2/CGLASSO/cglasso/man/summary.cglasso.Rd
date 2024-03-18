\name{summary.cglasso}
\alias{summary.cglasso}

\title{Summarizing cglasso and cggm Fits}

\description{
\sQuote{\code{summary}} produces a summary of the models fitted using \code{cglasso} or \code{cggm}.
}

\usage{
\method{summary}{cglasso}(object, GoF = AIC, print.all = TRUE, digits = 3L, \dots)
}

\arguments{
\item{object}{an \R object inheriting class \code{cglasso}, that is, the output of the model-fitting functions \code{\link{cglasso}} and \code{\link{cggm}}.}
\item{GoF}{a valid goodness-of-fit function, such as \code{\link{AIC.cglasso}} or \code{\link{BIC.cglasso}}, or an \R object of class \sQuote{\code{GoF}}.}
\item{print.all}{logical. Should all summary statistics be printed?}
\item{digits}{the minimum number of significant digits to be used: see \sQuote{\link[base]{print.default}}.}
\item{\dots}{the \emph{penalty} parameter passed to the goodness-of-fit function (argument \sQuote{\code{GoF}}) or further arguments passed to \sQuote{\link[base]{print.data.frame}} and \sQuote{\link[base]{print.listof}}.}
}

\details{
The function \code{summary.cglasso} computes the summary statistics needed to evaluate the goodness-of-fit of the models fitted by \code{cglasso} or \code{cggm}. 

Model evaluation can be made in two ways. The easiest way is to use a valid goodness-of-fit function, such as \code{\link{AIC.cglasso}} or \code{\link{BIC.cglasso}}. In this case, further arguments can be passed to these functions by the argument \sQuote{\dots}. The second way consists in passing the output of a goodness-of-fit function, that is, an \R object of class \sQuote{\code{GoF}}. Usually, this approach is preferable when the computation of  the chosen goodness-of-fit measure is time-consuming, such as when the sample size is small relative to the number of parameters and the \code{AIC} or \code{BIC} functions are used to evaluate a sequence of fitted models. In these cases, we suggest the computation of several measures of goodness-of-fit in a preliminary step and then the use of \code{summary.cglasso} in a subsequent step to evaluate the sequence of fitted models.

To improve the readability of the results, the output is divided into two sections. 

The first section is structured in a sequence of tables showing, for each combination of the two tuning parameters \sQuote{\code{lambda}} and \sQuote{\code{rho}}, the number of estimated non-zero regression coefficients (\sQuote{\code{df.B}}), the number of estimated non-zero partial correlation coefficients (\sQuote{\code{df.Tht}}), the degrees-of-freedom (\sQuote{\code{df}}) and the value of the goodness-of-fit measure used to evaluate the fitted models. To help the user with the identification of the optimal fitted model, the last column of each table reports the ranking of the models (where the optimal model is marked with the symbol \sQuote{\code{<-}}).

The second section reports the summary statistics of the selected optimal model.
}

\value{
A named list with the following elements is silently returned:
\item{table}{data.frame containing the summary statistics used to evaluate the sequence of fitted models.}
\item{lambda.id}{position of the optimal \eqn{\lambda}{lambda}-value identified by the chosen goodness-of-fit function.}
\item{rho.id}{position of the optimal \eqn{\rho}{rho}-value identified by the chosen goodness-of-fit function.}
}

\author{Luigi Augugliaro (\email{luigi.augugliaro@unipa.it})}

\seealso{
\code{\link{cglasso}}, \code{\link{cggm}}, \code{\link{AIC.cglasso}} and \code{\link{BIC.cglasso}}.
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
summary(out)                                   # models evaluation by AIC
summary(out, GoF = BIC)                        # models evaluation by BIC
summary(out, GoF = BIC, mle = TRUE, g = 0.5)   # models evaluation by eBIC

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

out <- cglasso(. ~ ., data = Z, lambda = 0.01)
summary(out)                                   # models evaluation by AIC
summary(out, GoF = BIC)                        # models evaluation by BIC
summary(out, GoF = BIC, mle = TRUE, g = 0.5)   # models evaluation by eBIC

out <- cglasso(. ~ ., data = Z, rho = 0.01)
summary(out)                                   # models evaluation by AIC
summary(out, GoF = BIC)                        # models evaluation by BIC
summary(out, GoF = BIC, mle = TRUE, g = 0.5)   # models evaluation by eBIC

out <- cglasso(. ~ ., data = Z)
summary(out)                                   # models evaluation by AIC
summary(out, GoF = BIC)                        # models evaluation by BIC
summary(out, GoF = BIC, mle = TRUE, g = 0.5)   # models evaluation by eBIC
}

\keyword{regression}
\keyword{models}
