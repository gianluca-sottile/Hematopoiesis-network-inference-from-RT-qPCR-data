\name{ShowStructure}
\alias{ShowStructure}

\title{Show Package Structure}

\description{
\code{ShowStructure} gives information about the package structure.
}

\usage{
ShowStructure(module = c("ALL", "DM", "MF", "MS", "NA"), description = TRUE,
              plot = TRUE)
}

\arguments{
\item{module}{a description of the module for which the description is required.}
\item{plot}{logical. If \code{TRUE} a graph showing the package structure is plotted.}
\item{description}{logical. If \code{TRUE} a description of the available function is printed on screen.}
}

\value{
\code{ShowStructure} silently returns a named list with components:
\item{description}{a named list containing the description of the available functions.}
\item{graph}{an \R object of class \sQuote{\code{igraph}} describing the package structure.}
}

\author{Luigi Augugliaro (\email{luigi.augugliaro@unipa.it})}

\keyword{package}
