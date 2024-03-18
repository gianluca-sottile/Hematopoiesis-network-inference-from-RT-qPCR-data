\name{qqcnorm}
\alias{qqcnorm}

\title{Quantile-Quantile Plots for a \code{datacggm} Object}

\description{
Creates a quantile-quantile plot for a censored Gaussian distribution.
}

\usage{
qqcnorm(x, which, max.plot = 1L, save.plot = FALSE, grdev = pdf, grdev.arg,
        main = "Censored Normal Q-Q Plot", xlab = "Theoretical Quantiles",
        ylab = "Sample Quantiles", plot.it = TRUE, plot.pch = c(2L, 20L),
        plot.col = c(2L, 1L), plot.cex = c(2L, 1L), abline = FALSE,
        line.col = "gray50", line.lwd = 1L, line.lty = 2L, \dots)
}

\arguments{
\item{x}{an object of class \sQuote{\code{datacggm}}.}
\item{which}{a vector of integers used to specify the response variables for which the histogram is required.}
\item{max.plot}{the maximum number of plots drawn in a single figure.}
\item{save.plot}{a logical variable or a string specifying the path of the directory where plots will be saved. Letting \sQuote{\code{save.plot = TRUE}}, the required plots will be saved as external files inside the current working directory. User can save these files on a specific working directory passing the absolute path through the argument \code{save.plot}. On exit, working directory will be always set the previous one.}
\item{grdev}{the graphics device used to save the required histograms on external files. See \sQuote{\code{\link[grDevices]{device}}} for more details.}
\item{grdev.arg}{additional parameters passed to the graphics device specified by \sQuote{\code{grdev}}.}
\item{main, xlab, ylab}{plot labels.}
\item{plot.it}{logical. Should the result be plotted?}
\item{plot.pch}{a pair of graphical parameters. The first entry specifies the symbol used for plotting the points associated to the censoring values \eqn{lo}{lo} and \eqn{up}{up} (by default a triangle), whereas the second entry specifies the symbol used for the remaining points (by default a black point).}
\item{plot.col}{a pair of graphical parameters. The first entry specifies the colour used for plotting the points associated to the censoring values \eqn{lo}{lo} and \eqn{up}{up} (by default a red triangle), whereas the second entry specifies the colour used for the remaining points (by default a black point).}
\item{plot.cex}{a pair of graphical parameters. The first entry specifies the size of the symbol used for plotting the points associated to the censoring values \eqn{lo}{lo} and \eqn{up}{up}, whereas the second entry specifies the size of the symbol used for the remaining points.}
\item{abline}{logical. Should the line \eqn{y = x}{y = x} be plotted?}
\item{line.col}{graphical parameter. If \sQuote{\code{abline = TRUE}}, then this argument specifies the colour of the line \eqn{y = x}{y = x}.}
\item{line.lwd}{graphical parameter. If \sQuote{\code{abline = TRUE}}, then this argument specifies the line width of the line \eqn{y = x}{y = x}.}
\item{line.lty}{graphical parameter. If \sQuote{\code{abline = TRUE}}, then this argument specifies the line type of the line \eqn{y = x}{y = x}.}
\item{\dots}{additional graphical parameter passed to \sQuote{\code{plot}}.}
}

\details{
\sQuote{\code{qqcnorm}} produces a censored normal QQ plot, that is, a graphical method for comparing the empirical distribution of a given response variable (specified by the argument \code{which}) to the censored Gaussian distribution, which is defined as:
\tabular{lcl}{
\eqn{\Phi((lo - \mu)/\sigma)}{\Phi((lo - \mu)/\sigma)} \tab if \tab \eqn{y\le lo}{y <= lo} \cr
\eqn{\phi((y - \mu)/\sigma) / \sigma}{\phi((lo - \mu)/\sigma) / \sigma} \tab if \tab \eqn{lo < y < up}{lo < y < up} \cr
\eqn{1 - \Phi((up - \mu)/\sigma)}{1 - \Phi((up - \mu)/\sigma)} \tab if \tab \eqn{y\ge up}{y >= up}
}
where \eqn{\phi}{\phi} and \eqn{\Phi}{\Phi} are the probability density function and the cumulative ditribution of the standard normal distribution, respectively, whereas \eqn{lo}{lo} and \eqn{up}{up} are the lower and upper censoring values, respectively.

The comparison is done by plotting the empirical quantiles (\eqn{y}{y}-coordinate) against the theoretical quantiles (\eqn{x}{x}-coordinate) of the censored Gaussian distribution, which are defined as follows:
\tabular{lcl}{
\eqn{lo}{lo} \tab if \tab \eqn{p \le \Phi((lo - \mu)/\sigma)}{p <= \Phi((lo - \mu)/\sigma)} \cr
\eqn{\mu + \sigma \Phi^{-1}(p)}{\mu + \sigma \Phi^{-1}(p)} \tab if \tab \eqn{\Phi((lo - \mu)/\sigma) < p < 1 - \Phi((up - \mu)/\sigma)}{\Phi((lo - \mu)/\sigma) < p < 1 - \Phi((up - \mu)/\sigma)} \cr
\eqn{up}{up} \tab if \tab \eqn{p \ge 1 - \Phi((up - \mu)/\sigma)}{p >= 1 - \Phi((up - \mu)/\sigma)} 
}
where \eqn{p\in(0, 1)}{0 < p < 1}. If the two distributions are similar, the points will approximately lie on the line \eqn{y = x}{y = x}. If the distributions are linearly related, the points will approximately lie on a line, but not necessarily on the line \eqn{y = x}{y = x}. In order to evaluate if the proportions of left/right-censored values are similar to the Gaussian tail probabilities, points corresponding to the censored values are plotted using a  specific symbol (see argument \sQuote{\code{plot.pch}}), colour (see argument \sQuote{\code{plot.col}}) and size (see argument \sQuote{\code{plot.cex}}).

Finally, maximum likelihood estimates of the marginal parameters \eqn{\mu}{\mu} and \eqn{\sigma}{\sigma} are computed as described in \sQuote{\code{\link{datacggm}}} and can be extracted from an \R of class \sQuote{\code{datacggm}} by using the functions \sQuote{\code{\link{ColMeans}}} and \sQuote{\code{\link{ColVars}}}, respectively.
}

\value{
A named list is silently returned. Each element of the list contains a two-columns matrix; first columns (named \sQuote{\code{x}}) contains the theoretical quantiles whereas second columns (named \sQuote{\code{y}}) contains the empirical quantiles.
}

\author{Gianluca Sottile (\email{gianluca.sottile@unipa.it})}

\references{
Augugliaro, L., Sottile, G., and Vinciotti, V. (2020) <\doi{10.1007/s11222-020-09945-7}>.
The conditional censored graphical lasso estimator.
\emph{Statistics and Computing} \bold{30}, 1273--1289.

Augugliaro, L., Abbruzzo, A., and Vinciotti, V. (2020) <\doi{10.1093/biostatistics/kxy043}>.
\eqn{\ell_1}{l1}-Penalized censored Gaussian graphical model.
\emph{Biostatistics} \bold{21}, e1--e16.
}

\seealso{
\code{\link{datacggm}}, \code{\link{rcggm}}, \code{\link{ColMeans}}, \code{\link{ColVars}} and \code{\link{hist.datacggm}}.
}

\examples{
set.seed(123)

# a dataset from a right-censored Gaussian graphical model
n <- 1000L
p <- 10L
Y <- matrix(rnorm(n * p), n, p)
up <- 1
Y[Y >= up] <- up
Z <- datacggm(Y = Y, up = up)
qqcnorm(Z, max.plot = 4L)

# a dataset from a  conditional censored Gaussian graphical model
n <- 1000L
p <- 10L
q <- 2
Y <- matrix(rnorm(n * p), n, p)
up <- 1
lo <- -1
Y[Y >= up] <- up
Y[Y <= lo] <- lo
X <- matrix(rnorm(n * q), n, q)
Z <- datacggm(Y = Y, lo = lo, up = up, X = X)
qqcnorm(Z, max.plot = 4L)
}

\keyword{hplot}

