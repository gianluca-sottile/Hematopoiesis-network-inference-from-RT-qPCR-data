\name{rcggm}
\alias{rcggm}

\title{Simulate Data from a Conditional Gaussian Graphical Model with Censored and/or Missing Values}

\description{
\sQuote{\code{rcggm}} function is used to produce one or more samples from a conditional Gaussian graphical model with censored and/or missing values.
}

\usage{
rcggm(n, p, b0, X, B, Sigma, probl, probr, probna, \dots)
}

\arguments{
\item{n}{the number of samples required (optional, see below for a description).}
\item{p}{the number of response variables (optional, see below for a description).}
\item{b0}{a vector of length \code{p} used to specify the intercepts. Default is a zero vector of length \code{p} (optional, see below for a description).}
\item{X}{a matrix of dimension \eqn{n\times q}{n x q} used to model the expected values of the response variables (optional, see below for a description).}
\item{B}{a matrix of dimension \eqn{q\times p}{q x p} used to specify the regression coefficients. If \code{X} is missing then \code{B} is set equal to \code{NULL} (optional, see below for a description).}
\item{Sigma}{a positive-definite symmetric matrix specifying the covariance matrix of the response variables. Default is the identity matrix (optional, see below for a description).}
\item{probl}{a vector giving the probabilities that the response variables are left-censored.}
\item{probr}{a vector giving the probabilities that the response variables are right-censored.}
\item{probna}{the probability that a response value is missing-at-random. By default \sQuote{\code{probna}} is set equal to zero.}
\item{\dots}{further arguments passed to the function \sQuote{\code{mvrnorm}}.}
}

\details{
\sQuote{The \code{rcggm}} function simulates a dataset from a conditional Gaussian graphical model with censored or missing values and returns an object of class \sQuote{\code{datacggm}}. Censoring values are implicitly specified using arguments \code{probl} and \code{probr}, that is, \code{lo} and \code{up} are computed in such a way that the average probabilities of left and right censoring are equal to \code{probl} and \code{probr}, respectively. Missing-at-random values are simulated using a Bernoulli distribution with probability \code{probna}.

The dataset is simulated through the following steps:
\enumerate{
\item lower and upper censoring values (\code{lo} and \code{up}) are computed according to the arguments \code{probl} and \code{probr};
\item The function \code{\link[MASS]{mvrnorm}} is used to simulate one or more samples from the multivariate Gaussian distribution specified by the arguments \code{b0}, \code{X}, \code{B} and \code{Sigma};
\item The response values that are outside of the interval \code{[lo, up]} are replaced with the corresponding censoring values;
\item if \code{probna} is greater than zero, then missing-at-random values are simulated using a Bernoulli distribution with probability \code{probna}.
} 

%The default setting of the multivariate Gaussian distribution used in step 2 depends on what arguments are specified. Table below sums up default models (specified arguments are marked by the symbol \sQuote{\code{x}}).

\tabular{cccccccl}{
Model \tab \code{n} \tab \code{p} \tab \code{b0} \tab \code{X} \tab \code{B} \tab \code{Sigma} \tab Gaussian distribution \cr
1 \tab \code{x} \tab \code{x} \tab           \tab          \tab          \tab              \tab \eqn{Y\sim N(0, I)}{Y ~ N(0, I)} \cr
2 \tab \code{x} \tab          \tab           \tab          \tab          \tab \code{x}     \tab \eqn{Y\sim N(0, \Sigma)}{Y ~ N(0, Sigma)} \cr
3 \tab \code{x} \tab          \tab \code{x}  \tab          \tab          \tab              \tab \eqn{Y\sim N(b0, I)}{Y ~ N(b0, I)} \cr
4 \tab \code{x} \tab          \tab \code{x}  \tab          \tab          \tab \code{x}     \tab \eqn{Y\sim N(b0, \Sigma)}{Y ~ N(b0, Sigma)} \cr
5 \tab          \tab          \tab           \tab\code{x}  \tab \code{x} \tab              \tab \eqn{Y\sim N(XB, I)}{Y ~ N(XB, I)} \cr
6 \tab          \tab          \tab           \tab\code{x}  \tab \code{x} \tab \code{x}     \tab \eqn{Y\sim N(XB, \Sigma)}{Y ~ N(XB, Sigma)} \cr
7 \tab          \tab          \tab \code{x}  \tab\code{x}  \tab \code{x} \tab              \tab \eqn{Y\sim N(b0 + XB, I)}{Y ~ N(b0 + XB, I)} \cr
8 \tab          \tab          \tab \code{x}  \tab\code{x}  \tab \code{x} \tab \code{x}     \tab \eqn{Y\sim N(b0 + XB, \Sigma)}{Y ~ N(b0 + XB, Sigma)}
}

The previous table sums up the default setting of the multivariate Gaussian distribution used in step 2 (specified arguments are marked by the symbol \sQuote{\code{x}}). See also below for some examples.
}

\value{
\code{rcggm} returns an object of class \sQuote{\code{datacggm}}. See \code{\link{datacggm}} for further details.
}

\author{Luigi Augugliaro (\email{luigi.augugliaro@unipa.it})}

\references{
Augugliaro, L., Sottile, G., and Vinciotti, V. (2020a) <\doi{10.1007/s11222-020-09945-7}>.
The conditional censored graphical lasso estimator.
\emph{Statistics and Computing} \bold{30}, 1273--1289.

Augugliaro, L., Abbruzzo, A., and Vinciotti, V. (2020b) <\doi{10.1093/biostatistics/kxy043}>.
\eqn{\ell_1}{l1}-Penalized censored Gaussian graphical model.
\emph{Biostatistics} \bold{21}, e1--e16.}

\seealso{
\code{\link{datacggm}}.
}

\examples{
set.seed(123)

n <- 100
p <- 3
q <- 2
b0 <- rep(1, p)
X <- matrix(rnorm(n * q), n, q)
B <- matrix(rnorm(q * p), q, p)
Sigma <- outer(1:p, 1:p, function(i, j) 0.3^abs(i - j))
probl <- 0.05
probr <- 0.05
probna <- 0.05

# Model 1: Y ~ N(0, I)
Z <- rcggm(n = n, p = p, probl = probl, probr = probr, probna = probna)
summary(Z)

# Model 2: Y ~ N(0, Sigma)
Z <- rcggm(n = n, Sigma = Sigma, probl = probl, probr = probr, probna = probna)
summary(Z)

# Model 3: Y ~ N(b0, I)
Z <- rcggm(n = n, b0 = b0, probl = probl, probr = probr, probna = probna)
summary(Z)

# Model 4: Y ~ N(b0, Sigma)
Z <- rcggm(n = n, b0 = b0, Sigma = Sigma, probl = probl, probr = probr, 
           probna = probna)
summary(Z)

# Model 5: Y ~ N(XB, I)
Z <- rcggm(X = X, B = B, probl = probl, probr = probr, probna = probna)
summary(Z)

# Model 6: Y ~ N(XB, Sigma)
Z <- rcggm(X = X, B = B, Sigma = Sigma, probl = probl, probr = probr, 
           probna = probna)
summary(Z)

# Model 7: Y ~ N(b0 + XB, I)
Z <- rcggm(b0 = b0, X = X, B = B, probl = probl, probr = probr, probna = probna)
summary(Z)

# Model 8: Y ~ N(b0 + XB, Sigma)
Z <- rcggm(b0 = b0, X = X, B = B, Sigma = Sigma, probl = probl, probr = probr, 
           probna = probna)
summary(Z)
}

\keyword{distribution}
\keyword{datagen}
\keyword{multivariate}

