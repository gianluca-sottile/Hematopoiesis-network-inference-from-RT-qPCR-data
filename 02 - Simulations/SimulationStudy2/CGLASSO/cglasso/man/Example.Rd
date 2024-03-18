\name{Example}
\docType{data}
\alias{Example}
\alias{MKMEP.Sim}
\alias{MM.Sim}

\title{
Simulated data for the cglasso vignette
}

\description{
Simulated data, used to demonstrate the features of datacggm class. It contains two simulated datasets that mimic the structural characteristics of the two real datasets MKMEP and MM. 
\enumerate{
\item \code{MKMEP.Sim}: is a reponse matrix containing 50 samples of 10 genes;
\item \code{MM.Sim}: is a list containing two components;
\itemize{
\item \code{Y}: is the response matrix containing 50 measurements of 10 miRNAs
\item \code{X}: is the predictor matrix containing 50 measurements of 5 variables (four numerical and one categorical)
}
}
In both datasets, the response matrices are right-censored with an upper limit of detection fixed to 40.
}

\usage{data("Example")}

\seealso{
\code{\link{cglasso}}, \code{\link{to_graph}}, and the method functions \code{\link{summary}}, \code{\link{coef}}, \code{\link{plot}}, \code{\link{AIC.cglasso}}, \code{\link{BIC.cglasso}}, \code{\link{MKMEP}} and \code{\link{MM}}.
}

\examples{
data("Example")

MM.Sim <- datacggm(Y = MM.Sim$Y, up = 40, X = MM.Sim$X)
ColMeans(MM.Sim)
ColVars(MM.Sim)
summary(MM.Sim)

MKMEP.Sim <- datacggm(Y = MKMEP.Sim, up = 40)
ColMeans(MKMEP.Sim)
ColVars(MKMEP.Sim)
summary(MKMEP.Sim)
}

\keyword{datasets}
