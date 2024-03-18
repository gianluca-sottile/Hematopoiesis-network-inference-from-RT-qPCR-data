\name{summary.datacggm}
\alias{summary.datacggm}

\title{Summarizing Objects of Class \sQuote{\code{datacggm}}}

\description{
The function \sQuote{\code{summary.datacggm}} produces summaries of an object of class \sQuote{\code{datacggm}}.
}

\usage{
\method{summary}{datacggm}(object, n, quantile.type = 7L, digits = 3L, quote = FALSE, \ldots)
}

\arguments{
\item{object}{an object of class \sQuote{\code{datacggm}} for which a summary is desired.}
\item{n}{maximum number of rows printed on video. By default all rows are printed.}
\item{quantile.type}{an integer between 1 and 9 selecting one of the nine quantile algorithms (see also \code{\link[stats]{quantile}} function).}
\item{digits}{integer, used for number formatting with \sQuote{\code{\link[base]{format}()}}.}
\item{quote}{logical, indicating whether or not strings should be printed with surrounding quotes.}
\item{\dots}{further arguments passed to \code{\link[base]{print.listof}}.}
}

\details{
The function \sQuote{\code{summary.datacggm}} extends the results given by \sQuote{\code{\link[base]{summary.matrix}()}} taking into account the specific structure of an object of class \sQuote{\code{datacggm}}. Summary statistics are printed out in two sections: the first section reports the summary statistics for the response matrix \sQuote{\code{Y}}, whereas the second section is devoted to summarising the matrix of predictors \sQuote{\code{X}}. 

For each response variable, the first section reports the standard summary statistics computed using only the observed values, that is the entries of the response matrix whose status indicators are equals to zero (see \sQuote{\code{\link{event}}} for more details), together with the information about the lower and upper censoring values (columns named \code{Lower} and \code{Upper}) and the percentage of missing-at-random, left-censored and right-censored response values (columns named \code{NA\%}, \code{LC\%} and \code{RC\%}, respectively).

The second section is divided into two subsections reporting the main summary statistics computed for numeric and categorical predictors, respectively.

The two sections are printed out in such a way that the readability of the summary statistics is improved in a high-dimensional setting, that is, when the number of predictors and response variables exceeds the sample size. 
}

\value{
\sQuote{\code{summary.dataccgm}} returns a named list with components
\item{Y}{a matrix with class \sQuote{\code{table}} obtained by computing the summary measures on the response variables.}
\item{X.numeric}{a matrix with class \sQuote{\code{table}} obtained by computing the summary measures on the numeric predictors.}
\item{X.categorical}{a list containing the summary measures on the categorical predictors.}
}

\author{Luigi Augugliaro (\email{luigi.augugliaro@unipa.it})}

\references{
Augugliaro, L., Sottile, G., and Vinciotti, V. (2020a) <\doi{10.1007/s11222-020-09945-7}>.
The conditional censored graphical lasso estimator.
\emph{Statistics and Computing} \bold{30}, 1273--1289.

Augugliaro, L., Abbruzzo, A., and Vinciotti, V. (2020b) <\doi{10.1093/biostatistics/kxy043}>.
\eqn{\ell_1}{l1}-Penalized censored Gaussian graphical model.
\emph{Biostatistics} \bold{21}, e1--e16.
}

\seealso{
\code{\link{datacggm}}, \code{\link{rcggm}} and \code{\link{event}}.
}

\examples{
set.seed(123)

# case 1: Y ~ N(0, Sigma) 
# 1. probability of left/right censored values equal to 0.05
# 2. probability of missing-at-random equals to 0.05
n <- 50L
p <- 100L
rho <- 0.3
Sigma <- outer(1L:p, 1L:p, function(i, j) rho^abs(i - j))
Z <- rcggm(n = n, Sigma = Sigma, probl = 0.05, probr = 0.05, probna = 0.05)
summary(Z, max = n * p)

# case 2: Y ~ N(b0 + XB, Sigma) and
# 1. probability of left/right censored values equal to 0.05
# 2. probability of missing-at-random equals to 0.05
n <- 50L
p <- 100L
q <- 100L
b0 <- runif(p)
B <- matrix(runif(q * p), nrow = q, ncol = p)
X <- matrix(rnorm(n * q), nrow = n, ncol = q)
rho <- 0.3
Sigma <- outer(1L:p, 1L:p, function(i, j) rho^abs(i - j))
Z <- rcggm(n = n, b0 = b0, X = X, B = B, Sigma = Sigma, probl = 0.05, probr = 0.05, 
           probna = 0.05)
summary(Z)
summary(Z, max = n * p)
}

\keyword{methods}

