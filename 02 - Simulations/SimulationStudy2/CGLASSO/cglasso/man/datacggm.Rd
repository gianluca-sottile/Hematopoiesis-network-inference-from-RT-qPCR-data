\name{datacggm}
\alias{datacggm}

\title{Create a Dataset from a Conditional Gaussian Graphical Model with Censored and/or Missing Values}

\description{
\sQuote{The \code{datacggm}} function is used to create a dataset from a conditional Gaussian graphical model with censored and/or missing values.
}

\usage{
datacggm(Y, lo = -Inf, up =  +Inf, X = NULL, control = list(maxit = 1.0E+4,
         thr = 1.0E-4))
}

\arguments{
\item{Y}{a \eqn{(n\times p)}{(n x p)}-dimensional matrix; each row is an observation from a conditional Gaussian graphical model with censoring vectors \code{lo} and \code{up}. Missing-at-random values are recorded as \sQuote{\code{NA}}.}
\item{lo}{the lower censoring vector; \code{lo[j]} is used to specify the lower censoring value for the random variable \eqn{Y_j}{Y_j}.}
\item{up}{the upper censoring vector; \code{up[j]} is used to specify the upper censoring value for the random variable \eqn{Y_j}{Y_j}.}
\item{X}{an optional \eqn{(n\times q)}{(n x q)}-dimensional data frame of predictors. If missing (default), a dataset from a Gaussian graphical model is returned otherwise a dataset from a conditional Gaussian graphical model is returned.}
\item{control}{a named list used to pass the arguments to the EM algorithm (see below for more details). The components are:
\itemize{
\item \code{maxit}: maximum number of iterations. Default is \code{1.0E+4}.
\item \code{thr}: threshold for the convergence. Default value is \code{1.0E-4}.}
}
}

\details{
The function \sQuote{\code{datacggm}} returns an \R object of class \sQuote{\code{datacggm}}, that is a named list containing the elements needed to fit a conditional graphical LASSO (cglasso) model to datasets with censored and/or missing values.

A set of specific method functions are developed to decsribe data with censored/missing values. For example, the method function \sQuote{\code{print.datacggm}}  prints out the left and right-censored values using the following rules: a right-censored value is labeled adding the symbol \sQuote{\code{+}} at the end of the value, whereas the symbol \sQuote{\code{-}} is used for the left-censored values (see examples below). The summary statistics can be obtained using the method function \sQuote{\code{\link{summary.datacggm}}}. The matrices \code{Y} and \code{X} are extracted from a \code{datacggm} object using the function \sQuote{\code{\link{getMatrix}}}.

For each column of the matrix \sQuote{\code{Y}}, mean and variance are estimated using a standard EM-algorithm based on the assumption of a Gaussian distribution. \sQuote{\code{maxit}} and \sQuote{\code{thr}} are used to set the number of iterations and the threshold for convergence, respectively. Marginal means and variances can be extracted using the accessor functions \sQuote{\code{\link{ColMeans}}} and \sQuote{\code{\link{ColVars}}}, respectively. Furthermore, the plotting functions \sQuote{\code{\link{hist.datacggm}}} and \sQuote{\code{\link{qqcnorm}}} can be used to inspect the marginal distribution of each column of the matrix \sQuote{\code{Y}}.

The status indicator matrix, denoted by \code{R}, can be extracted by using the function \code{\link{event}}. The entries of this matrix specify the status of an observation using the following code:
\itemize{
\item \sQuote{\code{R[i, j] =  0}} means that the \eqn{y_{ij}}{y_{ij}} is inside the open interval \code{(lo[j], up[j])};
\item \sQuote{\code{R[i, j] = -1}} means that the \eqn{y_{ij}}{y_{ij}} is a left-censored value;
\item \sQuote{\code{R[i, j] = +1}} means that the \eqn{y_{ij}}{y_{ij}} is a right-censored value;
\item \sQuote{\code{R[i, j] = +9}} means that the \eqn{y_{ij}}{y_{ij}} is a missing value.
}

See below for the other functions related to an object of class \sQuote{\code{datacggm}}.
}

\value{
\sQuote{\code{datacggm}} returns an \R object of S3 class \dQuote{\code{datacggm}}, that is, a nested named list containing the
following components:
\item{Y}{the \eqn{(n\times p)}{(n x p)}-dimensional matrix \code{Y}.}
\item{X}{the \eqn{(n\times q)}{(n x q)}-dimensional data frame \code{X}.}
\item{Info}{
	\itemize{
		\item \code{lo:} the lower censoring vector;
		\item \code{up:} the upper censoring vector;
		\item \code{R:} the status indicator matrix encoding the censored/missing values (mainly for internal purposes);
		\item \code{order:} an integer vector used for the ordering of the matrices \code{Y} and \code{X} (for internal purposes only);
		\item \code{Pattern:} a matrix encoding the information about the the patterns of censored/missing values (for internal purposes only);
		\item \code{ym:} the estimated marginal means of the random variables \eqn{Y_j}{Y_j};
		\item \code{yv:} the estimated marginal variances of the random variables \eqn{Y_j}{Y_j};
		\item \code{n:} the sample size;
		\item \code{p:} the number of response variables;
		\item \code{q:} the number of columns of the data frame \code{X}.
	}
}
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
Related to the \R objects of class \dQuote{\code{datacggm}} there are the accessor functions, \code{\link{rowNames}}, \code{\link{colNames}}, \code{\link{getMatrix}}, \code{\link{ColMeans}}, \code{\link{ColVars}}, \code{\link{upper}}, \code{\link{lower}}, \code{\link{event}}, \code{\link{qqcnorm}} and the method functions \code{\link{is.datacggm}}, \code{\link{dim.datacggm}}, \code{\link{summary.datacggm}} and \code{\link{hist.datacggm}}. The function \code{\link{rcggm}} can be used to simulate a dataset from a conditional Gaussian graphical model whereas \code{\link{cglasso}} is the model fitting function devoted to the l1-penalized censored Gaussian graphical model.
}

\examples{
set.seed(123)

# a dataset from a right-censored Gaussian graphical model
n <- 100L
p <- 3L
Y <- matrix(rnorm(n * p), n, p)
up <- 1
Y[Y >= up] <- up
Z <- datacggm(Y = Y, up = up)
Z

# a dataset from a  conditional censored Gaussian graphical model
n <- 100L
p <- 3L
q <- 2
Y <- matrix(rnorm(n * p), n, p)
up <- 1
lo <- -1
Y[Y >= up] <- up
Y[Y <= lo] <- lo
X <- matrix(rnorm(n * q), n, q)
Z <- datacggm(Y = Y, lo = lo, up = up, X = X)
Z

# a dataset from a  conditional censored Gaussian graphical model 
# and with missing-at-random values
n <- 100L
p <- 3L
q <- 2
Y <- matrix(rnorm(n * p), n, p)
NA.id <- matrix(rbinom(n * p, 1L, 0.01), n, p)
Y[NA.id == 1L] <- NA
up <- 1
lo <- -1
Y[Y >= up] <- up
Y[Y <= lo] <- lo
X <- matrix(rnorm(n * q), n, q)
Z <- datacggm(Y = Y, lo = lo, up = up, X = X)
Z
}

\keyword{classes}
\keyword{methods}
