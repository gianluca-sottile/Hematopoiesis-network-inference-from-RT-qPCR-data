\name{plot.GoF}
\alias{plot.GoF}

\title{Plot for \sQuote{\code{GoF}} Object}

\description{
\sQuote{The \code{plot.GoF}} function produces plots to study the sequence of fitted models.
}

\usage{
\method{plot}{GoF}(x, add.line = TRUE, arg.line = list(lty = 2L, lwd = 2L, col = "red"), 
      add.text = FALSE, arg.text = list(side = 3L), arg.points = list(pch = 2L),
      \dots)
}

\arguments{
\item{x}{an \R object of class \sQuote{\code{GoF}}, that is, the output of a goodness-of-fit function such as \code{\link{AIC.cglasso}} or \code{\link{BIC.cglasso}}.}
\item{add.line}{logical; if \sQuote{\code{add.line = TRUE}} then a line is added to identify the optimal value of the tuning parameter.}
\item{arg.line}{a named list of graphical parameters passed to the function \code{\link[graphics]{abline}} (see also \code{\link[graphics]{par}}).}
\item{add.text}{logical; if \sQuote{\code{add.text = TRUE}} then a text is added to the line used to identify the optimal value of the tuning parameter.}
\item{arg.text}{a list of further parameters passed to the function \code{\link[graphics]{mtext}} (only if \sQuote{\code{add.text = TRUE}}).}
\item{arg.points}{a named list of graphical parameters passed to the function \code{\link[graphics]{points}}.}
\item{\dots}{additional graphical arguments passed to the functions \code{\link[graphics]{plot}},  \code{\link[graphics]{contour}} or \code{\link[graphics]{filled.contour}}.}
}

\details{
\code{plot.GoF} is the plotting method function of an \R object of class \sQuote{\code{GoF}}, that is, the output of a goodness-of-fit function (see \code{\link{AIC.cglasso}}, or \code{\link{BIC.cglasso}}). This function produces a plot aimed both to evaluate the sequence of fitted models in terms of goodness-of-fit and to identify the optimal values of the tuning parameters.

If a tuning parameter is held fixed, then \code{plot.GoF} produces a plot showing the chosen measure of goodness-of-fit as a function of the remaining tuning parameter. In this case, the optimal value is identified by a vertical dashed line. The degrees-of-freedom of the selected fitted model are also shown.

If the cglasso model is fitted using both a sequence of \eqn{\rho}{rho} and \eqn{\lambda}{lambda} values, then \code{plot.GoF} produces a contour plot and a triangle is used to identify the optimal pair of the two tuning parameters.
}

\author{Luigi Augugliaro (\email{luigi.augugliaro@unipa.it})}

\seealso{
\code{\link{cglasso}}, \code{\link{AIC.cglasso}}, \code{\link{BIC.cglasso}}, \code{\link{summary.cglasso}} and \code{\link{select.cglasso}}.
}

\examples{
set.seed(123)
n <- 1000L
p <- 3L
q <- 2
b0 <- runif(p)
B <- matrix(runif(q * p), nrow = q, ncol = p)
X <- matrix(rnorm(n * q), nrow = n, ncol = q)
rho <- 0.3
Sigma <- outer(1L:p, 1L:p, function(i, j) rho^abs(i - j))
Z <- rcggm(n = n, b0 = b0, X = X, B = B, Sigma = Sigma, probl = 0.05, probr = 0.05)

out <- cglasso(. ~ ., data = Z, nlambda = 1L)
plot(AIC(out))
plot(BIC(out))

out <- cglasso(. ~ ., data = Z, nrho = 1L)
plot(AIC(out))
plot(BIC(out))

out <- cglasso(. ~ ., data = Z)
plot(AIC(out))
plot(BIC(out))
}

\keyword{graphs}
