\name{MKMEP}
\docType{data}
\alias{MKMEP}

\title{
Megakaryocyte-Erythroid Progenitors
}

\description{
In a study about the formation of blood cells, Psaila \emph{and others} (2016) have recently identified three distinct sub-populations of cells, which are all derived from hematopoietic stem cells through cell differentiation. One of these sub-populations, denoted by MK-MEP, is a previously unknown, rare population of cells that are bipotent but primarily generate megakaryocytic progeny.

Multiplex RT-qPCR has been used to profile 63 genes and 48 single human MK-MEP cells. RT-qPCR data are typically right-censored with a limit of detection fixed by the manufacturer to 40. Raw data have been mean normalized using the method proposed in Pipelers \emph{and others} (2017). See Section 5 in Augugliaro \emph{and others} (2018) for more details.

\sQuote{\code{MKMEP}} is a matrix containing a subset of the data available from Psaila \emph{and others} (2016).
}

\usage{data("MKMEP")}

\references{
Augugliaro, L., Sottile, G., and Vinciotti, V. (2020) <\doi{10.1007/s11222-020-09945-7}>.
The conditional censored graphical lasso estimator.
\emph{Statistics and Computing} \bold{30}, 1273--1289.

Augugliaro, L., Abbruzzo, A., and Vinciotti, V. (2020) <\doi{10.1093/biostatistics/kxy043}>.
\eqn{\ell_1}{l1}-Penalized censored Gaussian graphical model.
\emph{Biostatistics} \bold{21}, e1--e16.

Pipelers, P., Clement, L., Vynck, M., Hellemans, J., Vandesompele, J. and Thas, O. (2017) <doi: 10.1371/journal.pone.0182832>.
A unified censored normal regression model for qPCR differential gene expression analysis.
\emph{PLoS One} \bold{12}, e0182832.

Psaila, B., Barkas, N., Iskander, D., Roy, A., Anderson, S., Ashley, N., Caputo, V. S., Lichtenberg, J., Loaiza, S., Bodine, D. M. \emph{and others}. (2016) <doi: 10.1186/s13059-016-0939-7>
Single-cell profiling of human megakaryocyte-erythroid progenitors identifies distinct megakaryocyte and erythroid differentiation pathways.
\emph{Genome Biology} \bold{17}, 83--102.
}

\seealso{
\code{\link{cglasso}}, \code{\link{to_graph}}, and the method functions \code{\link{summary}}, \code{\link{coef}}, \code{\link{plot}}, \code{\link{AIC.cglasso}} and \code{\link{BIC.cglasso}}.
}

\examples{
data("MKMEP")
is.matrix(MKMEP)
Z <- datacggm(Y = MKMEP, up = 40)
str(Z)
is.datacggm(Z)
print(Z, width = 80L)
}

\keyword{datasets}
