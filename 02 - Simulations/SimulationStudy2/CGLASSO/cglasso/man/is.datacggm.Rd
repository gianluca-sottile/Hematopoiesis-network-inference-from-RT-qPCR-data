\name{is.datacggm}
\alias{is.datacggm}

\title{Is an Object of Class \sQuote{datacggm}?}

\usage{
is.datacggm(x)
}

\description{
\sQuote{\code{is.datacggm}} tests if its argument is an object of class \sQuote{datacggm}.
}

\arguments{
  \item{x}{object to be tested.}
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
\code{\link{datacggm}} and \code{\link{rcggm}}
}

\examples{
set.seed(123)
n <- 100L
p <- 3L
b0 <- rep(0, p)
Z <- rcggm(n = n, b0 = b0, probl = 0.05, probr = 0.05)
is.datacggm(Z) # TRUE
}

\keyword{methods}
\keyword{classes}
