\name{cglasso-package}
\alias{cglasso-package}

\docType{package}

\title{Conditional Graphical LASSO for Gaussian Graphical Models with Censored and Missing Values}

\description{Conditional graphical lasso (cglasso) estimator (Yin \emph{and other}, 2011) is an extension of the graphical lasso (Yuan \emph{and other}, 2007) proposed to estimate the conditional dependence structure of a set of p response variables given q predictors. This package provides suitable extensions developed to study datasets with censored and/or missing values (Augugliaro \emph{and other}, 2020a and Augugliaro \emph{and other}, 2020b). Standard conditional graphical lasso is available as a special case. Furthermore, the cglasso package provides an integrated set of core routines for visualization, analysis, and simulation of datasets with censored and/or missing values drawn from a Gaussian graphical model.}

\details{
\tabular{ll}{
Package: \tab cglasso\cr
Type: \tab Package\cr
Version: \tab 2.0.5\cr
Date: \tab 2021-06-08\cr
License: \tab GPL (>=2)
}
}

\author{
Luigi Augugliaro [aut, cre],\cr
Gianluca Sottile [aut]\cr
Ernst C. Wit [aut]\cr
Veronica Vinciotti [aut]\cr
Maintainer: Luigi Augugliaro <\email{luigi.augugliaro@unipa.it}>
}

\references{
Augugliaro, L., Sottile, G., and Vinciotti, V. (2020a) <\doi{10.1007/s11222-020-09945-7}>.
The conditional censored graphical lasso estimator.
\emph{Statistics and Computing} \bold{30}, 1273--1289.

Augugliaro, L., Abbruzzo, A., and Vinciotti, V. (2020b) <\doi{10.1093/biostatistics/kxy043}>.
\eqn{\ell_1}{l1}-Penalized censored Gaussian graphical model.
\emph{Biostatistics} \bold{21}, e1--e16.

Friedman, J.H., Hastie, T., and Tibshirani, R. (2008) <\doi{10.1093/biostatistics/kxm045}>.
Sparse inverse covariance estimation with the graphical lasso.
\emph{Biostatistics} \bold{9}, 432--441.

Yin, J. and Li, H. (2001) <\doi{10.1214/11-AOAS494}>.
A sparse conditional Gaussian graphical model for analysis of genetical genomics data.
\emph{The Annals of Applied Statistics} \bold{5}(4), 2630--2650.

Yuan, M., and Lin, Y. (2007) <\doi{10.1093/biomet/asm018}>.
Model selection and estimation in the Gaussian graphical model.
\emph{Biometrika} \bold{94}, 19--35.
}

\keyword{multivariate}
\keyword{package}
\keyword{models}
