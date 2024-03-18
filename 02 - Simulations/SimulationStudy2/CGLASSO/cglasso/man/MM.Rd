\name{MM}
\docType{data}
\alias{MM}

\title{
The Rule of miRNA in Multiple Myeloma
}

\description{
MicroRNAs (miRNAs) are endogenous small non-coding RNAs, approximately 22 nucleotides in length, that play a regulatory role in gene expression by mediating mRNA cleavage or translation expression. Several studies have shown that a deregulation of the miRNAs can cause a disruption in the gene regulation mechanisms of the cell and that this might even lead to cancerous phenotypes.

Gutierrez \emph{and others} (2010) investigated the expression level of a set of 141 miRNAs, measured by RT-qPCR technology on a sample of 64 patients with multiple myeloma (MM). Patients were selected to represent the most relevant recurrent genetic abnormalities in MM. RT-qPCR data are typically right-censored and, in this study, the upper limit of detection was fixed to 40 cycles.

Formally, \code{MM} is a list with containing a subset of the date studied in Gutierrez \emph{and others} (2010). Its components are:
\itemize{
\item \code{Y}: the response matrix containing the cycle-threshold of the measured 49 miRNAs;
\item \code{X}: the predictor matrix containing the expression levels of 15 endogenous internal reference genes, called housekeeping, and a factor encoding the cytogenetic abnormalities.
}
}

\usage{data("MM")}

\references{
Augugliaro, L., Sottile, G., and Vinciotti, V. (2020) <\doi{10.1007/s11222-020-09945-7}>.
The conditional censored graphical lasso estimator.
\emph{Statistics and Computing} \bold{30}, 1273--1289.

Augugliaro, L., Abbruzzo, A., and Vinciotti, V. (2020) <\doi{10.1093/biostatistics/kxy043}>.
\eqn{\ell_1}{l1}-Penalized censored Gaussian graphical model.
\emph{Biostatistics} \bold{21}, e1--e16.

Gutierrez, N. Sarasquete, M., Misiewicz-Krzeminska, I., Delgado, M. De Las Rivas, J., Ticona, F.V., Ferminan, E., Martin-Jimenez, P., Chillon, C., Risueno, A., Hernandez, J.M., Garcia-Sanz, R., Gonzalez, M. and San Miguel, J.F. (2010) <\doi{10.1038/leu.2009.274}>.
Deregulation of microRNA expression in the different genetic subtypes of multiple myeloma and correlation with gene expression profiling.
\emph{Leukemia} \bold{24}(3), 629--637.
}

\seealso{
\code{\link{cglasso}}, \code{\link{to_graph}}, and the method functions \code{\link{summary}}, \code{\link{coef}}, \code{\link{plot}}, \code{\link{AIC.cglasso}} and \code{\link{BIC.cglasso}}.
}

\examples{
data("MM")
Z <- datacggm(Y = MM$Y, X = MM$X, up = 40)
print(Z, width = 80L)
}

\keyword{datasets}
