\name{select.cglasso}
\alias{select.cglasso}

\title{Model Selection for the Conditional Graphical Lasso Estimator}

\description{
\sQuote{\code{select.cglasso}} returns the optimal fitted model selected by a chosen measure of goodness-of-fit.
}

\usage{
select.cglasso(object, GoF = AIC, \dots)
}

\arguments{
\item{object}{an \R object of class \code{cglasso}.}
\item{GoF}{a valid goodness-of-fit function, such as \code{\link{AIC.cglasso}} or \code{\link{BIC.cglasso}}, or an \R object of class \sQuote{\code{GoF}}.}
\item{\dots}{further arguments passed to the chosen goodness-of-fit function (argument \sQuote{\code{GoF}}).}
}

\details{
The function \code{select.cglasso} evaluates the goodness-of-fit of the models fitted by \code{cglasso} and extracts the selected model.

Model evaluation can be made in two ways. The easiest way is to use a valid goodness-of-fit function, such as \code{\link{AIC.cglasso}} or \code{\link{BIC.cglasso}}. In this case, further arguments are passed to these functions by \sQuote{\dots}. The second way consists on passing the output of a goodness-of-fit function, that is, an \R object of class \sQuote{\code{GoF}}. Usually, this approach is preferable when the computation of the  chosen goodness-of-fit measure is time-consuming, such as when the sample size is small relative to the number of parameters and the \code{AIC} or \code{BIC} functions are used to evaluate a long sequence of fitted models. In these cases, we suggest the computation of several measures of goodness-of-fit in a preliminary step and then the use of the \code{select.cglasso} function in a subsequent step to select the optimal fitted model.
}

\value{
\sQuote{\code{select.cglasso}} returns the optimal fitted model.
}

\author{Gianluca Sottile (\email{gianluca.sottile@unipa.it})}

\seealso{
\code{\link{cglasso}}, \code{\link{AIC.cglasso}} and \code{\link{BIC.cglasso}}.
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
select.cglasso(out)                                    # models selection by AIC
select.cglasso(out, GoF = BIC)                         # models selection by BIC
select.cglasso(out, GoF = BIC, mle = TRUE, g = 0.5)    # models selection by eBIC

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
select.cglasso(out)                                    # models selection by AIC
select.cglasso(out, GoF = BIC)                         # models selection by BIC
select.cglasso(out, GoF = BIC, mle = TRUE, g = 0.5)    # models selection by eBIC

out <- cglasso(. ~ ., data = Z, rho = 0.01)
select.cglasso(out)                                    # models selection by AIC
select.cglasso(out, GoF = BIC)                         # models selection by BIC
select.cglasso(out, GoF = BIC, mle = TRUE, g = 0.5)    # models selection by eBIC

out <- cglasso(. ~ ., data = Z)
select.cglasso(out)                                    # models selection by AIC
select.cglasso(out, GoF = BIC)                         # models selection by BIC
select.cglasso(out, GoF = BIC, mle = TRUE, g = 0.5)    # models selection by eBIC
}

\keyword{regression}
\keyword{models}
