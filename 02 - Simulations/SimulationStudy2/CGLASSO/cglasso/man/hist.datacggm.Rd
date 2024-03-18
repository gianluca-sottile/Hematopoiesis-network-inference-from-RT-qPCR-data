\name{hist.datacggm}
\alias{hist}
\alias{hist.datacggm}

\title{Histogram for a \code{datacggm} Object}

\description{
Creates histograms with a normal density as background and two segments corresponding to left and/or right censored values (if any) in a conditional censored Gaussian graphical model. 
}

\usage{
\method{hist}{datacggm}(x, breaks = "Sturges", include.lowest = TRUE, right = TRUE, nclass = NULL,
     which, max.hist = 1L, save.hist = FALSE, grdev = pdf, grdev.arg, 
     polygon.col = adjustcolor("grey", alpha.f = 0.25), polygon.border = NA, 
     segments.lwd = 4L, segments.lty = 2L, segments.col = "gray40", 
     points.pch = c(4L, 1L), points.cex = 1.8, points.col = rep("black", 2L), 
     legend = TRUE, \dots)
}

\arguments{
\item{x}{an object of class \sQuote{\code{datacggm}}.}
\item{breaks}{one of:
\itemize{
\item{a vector giving the breakpoints between histogram cells,}
\item{a function to compute the vector of breakpoints,}
\item{a single number giving the number of cells for the histogram,}
\item{a character string naming an algorithm to compute the number of cells,}
\item{a function to compute the number of cells.}
}
See \code{\link{hist}} for more details.
}
\item{include.lowest}{logical; if \code{TRUE}, an \code{x[i]} equal to the breaks value will be included in the first (or last, for \code{right = FALSE}) bar. This will be ignored (with a warning) unless \code{breaks} is a vector.}
\item{right}{logical; if \code{TRUE}, the histogram cells are right-closed (left open) intervals.}
\item{nclass}{numeric (integer). For S(-PLUS) compatibility only, \code{nclass} is equivalent to \code{breaks} for a scalar or character argument.}
\item{which}{a vector of integers used to specify the response variables for which the histogram is required.}
\item{max.hist}{the maximum number of histograms drawn in a single figure.}
\item{save.hist}{a logical variable or a string specifying the path of the directory where plots will be saved. Letting \sQuote{\code{save.plot = TRUE}}, the required plots will be saved as external files inside the current working directory. User can save these files on a specific working directory passing the absolute path through the argument \code{save.plot}. On exit, working directory will be always set the previous one.}
\item{grdev}{the graphics device used to save the required histograms on external files. See \sQuote{\code{\link[grDevices]{device}}} for more details.}
\item{grdev.arg}{additional parameters passed to the graphics device specified by \sQuote{\code{grdev}}.}
\item{polygon.col}{graphical parameter; the colour for filling the area underlying the Gaussian distribution.}
\item{polygon.border}{graphical parameter; the colour of the border of the Gaussian distribution.}
\item{segments.lwd}{graphical parameter; the line width used in plotting the segments associated to the symbols specified by \sQuote{\code{points.pch}}.}
\item{segments.lty}{graphical parameter; the line type used in plotting the segments associated to the symbols specified by \sQuote{\code{points.pch}}.}
\item{segments.col}{graphical parameter; the colour used in plotting the segments associated to the symbols specified by \sQuote{\code{points.pch}}.}
\item{points.pch}{graphical parameter; a pair of integers specifying the symbols used to plot the proportion of censored response variables (default \sQuote{\code{points.pch = 4L}}) and the probability to observe a censored value (default \sQuote{\code{points.pch = 1L}}), respectively.}
\item{points.cex}{graphical parameter; a numerical value giving the amount by which the symbols, specified by \sQuote{\code{points.pch}}, are magnified relative to the default.}
\item{points.col}{graphical parameter; the colours used to plot the symbols specified by \sQuote{\code{points.pch}}.}
\item{legend}{logical; if \sQuote{\code{TRUE}} legend is added to the plot.}
\item{\dots}{additional graphical parameters passed to \sQuote{\code{\link[graphics]{plot.histogram}}}.}
}

\details{
For each response variable, the method function \sQuote{\code{hist.datacggm}} plots a histogram using only the observed values. Densities, plotted on the \code{y}-axis, are computed in such a way that the sum of all densities plus the proportion of censored values is equal to one. 

To evaluate the distributional assumption underlying the censored Gaussian distribution, on each histogram the area underlying the Gaussian density function is also shown in the background (marginal parameters are estimated as described in \sQuote{\code{\link{datacggm}}} and \sQuote{\code{\link{ColMeans}}}). Furthermore, on each plot, the proportions of left/right censored values are graphically compared with the corresponding Gaussian tail probabilities. If the assumption about the censoring mechanism is satisfied, then the proportion of censored values and the tail probability are approximately equals to each other.
}

\author{Gianluca Sottile (\email{gianluca.sottile@unipa.it})}

\seealso{
\code{\link{datacggm}}, \code{\link{rcggm}}, \code{\link{ColMeans}}, \code{\link{ColVars}} and \code{\link{qqcnorm}}.
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
hist(Z, max.hist = 4L)
}
\keyword{hplot}

