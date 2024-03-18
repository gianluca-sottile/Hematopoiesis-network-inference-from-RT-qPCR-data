\name{plot.cglasso}
\alias{plot.cglasso}

\title{Plot Method for \sQuote{\code{cglasso}} Object}

\description{
\sQuote{The \code{plot.cglasso}} function produces plots to study the coefficient paths of fitted cglasso models.
}

\usage{
\method{plot}{cglasso}(x, what = c("Theta", "diag(Theta)", "b0", "B"), 
     penalty = ifelse(x$nrho >= x$nlambda, "rho", "lambda"), given = NULL,
     GoF = AIC, add.labels, matplot.arg1, matplot.arg2, labels.arg, abline.arg,
     mtext.arg, save.plot, grdev = pdf, grdev.arg, digits = 4L, \dots)
}

\arguments{
\item{x}{an \R object of class \sQuote{\code{cglasso}}, that is, the output of the fitting function \code{\link{cglasso}}.}
\item{what}{a character or a formula specifying the required plot. If \code{what} is a character, then it is used to specify the conditional coefficient path. Allowed descriptors are: \sQuote{\code{Theta}}, the off-diagonal entries of the precision matrix, \sQuote{\code{diag(Theta)}}, the diagonal entries of the precision matrix, \sQuote{\code{b0}}, the intercepts of the conditional models, or equivalently the expected values of the response variables, and \sQuote{\code{B}}, the regression coefficients. Optionally, \sQuote{\code{what}} can be a model formula (see sections Description and Examples for more details).}
\item{penalty}{optional character argument used to specify the tuning parameter needed to plot the conditional coefficient path. Allowed descriptors are \sQuote{\code{rho}} and \sQuote{\code{lambda}}. This argument can be omitted if \code{what} is a model formula.}
\item{given}{an optional vector of integers identifying the conditioning values of the second tuning parameter to be used for the  coefficient profile plots of the parameters identified by the input \sQuote{\code{what}} across their corresponding tuning parameter. This argument is required only if \sQuote{\code{what}} is a character, otherwise it can be omitted.}
\item{GoF}{a valid goodness-of-fit function, such as \code{\link{AIC.cglasso}} and \code{\link{BIC.cglasso}}, or an \R object of class \sQuote{\code{GoF}}.}
\item{add.labels}{logical value. Should labels be added to the plot?}
\item{matplot.arg1}{a named list with the graphical parameters used to plot the paths of the estimates identified by \sQuote{\code{GoF}}.}
\item{matplot.arg2}{a named list with the graphical parameters used to plot the paths of the remaining estimates.}
\item{labels.arg}{a named list with the graphical parameters used to plot the labels.}
\item{abline.arg}{a named list with the graphical parameters used to plot the line identifying the optimal tuning parameter value.}
\item{mtext.arg}{a named list with the graphical parameters used to plot the text reported on the third axis.}
\item{save.plot}{a logical variable or a string specifying the path of the directory where plots will be saved. Letting \sQuote{\code{save.plot = TRUE}}, the required plots will be saved as external files inside the current working directory. User can save these files on a specific working directory passing the absolute path through the argument \code{save.plot}. On exit, working directory will be always set the previous one.}
\item{grdev}{the graphics device used to save the required histograms on external files. See \sQuote{\code{\link[grDevices]{device}}} for more details.}
\item{grdev.arg}{additional parameters passed to the graphics device specified by \sQuote{\code{grdev}}.}
\item{digits}{number of digits used to print the value of the second tuning parameter identified by \sQuote{\code{given}}.}
\item{\dots}{further arguments passed to the chosen goodness-of-fit function, such as \sQuote{\code{AIC}} or \sQuote{\code{BIC}}.}
}

\details{
The function \code{plot.cglasso} produces the coefficient profile plot. The output depends both on the type of fitted model and on the setting of the three main arguments, \code{what}, \code{penalty} and \code{given}. Below we give more details.

If the model fitting function \code{\link{cglasso}} is used to fit an l1-penalized censored Gaussian graphical model (see first part in Section Example), then the user can specify only the main argument \sQuote{\code{what}}, whereas \sQuote{\code{penalty}} and \sQuote{\code{given}} can be omitted. If \sQuote{\code{penalty}} is specified, then it must be equal to \dQuote{\code{rho}}. The main argument \sQuote{\code{what}} is used to specify the estimator that is plotted on the \eqn{y}{y}-axis. In this case, it must be equal to one of the following descriptors:
\itemize{
\item \dQuote{\code{Theta}}: the path of the estimated partial correlation coefficients is returned;
\item \dQuote{\code{diag(Theta)}}: the path of the diagonal entries of the estimated  precision matrix is returned;
\item \dQuote{\code{b0}}: the path of the estimated  expected values of the response variables is returned;
}

If a l1-penalized conditional censored Gaussian graphical model is fitted (see second part in Section Example), then all the main arguments can be used. In this case:
\enumerate{
\item \sQuote{\code{what}} must be equal to one of the following descriptors:
\itemize{
\item \dQuote{\code{Theta}}: the path of the off-diagonal entries of the estimated precision matrix is returned;
\item \dQuote{\code{diag(Theta)}}: the path of the estimated partial correlation coefficients is returned;
\item \dQuote{\code{b0}}: the path of the estimated intercepts is returned;
\item \dQuote{\code{B}}: the path of the estimated regression coefficients is returned for each response variable;
}
\item \sQuote{\code{penalty}} is used to specify the tuning parameter that is reported on the \eqn{x}{x}-axis and it must be equal to one of the following descriptors:
\itemize{
\item \dQuote{\code{rho}}: the tuning parameter used to penalize the precision matrix;
\item \dQuote{\code{lambda}}: the tuning parameter used to penalize the regression coefficient matrix;
}
\item \sQuote{\code{given}} is an optional vector of integers used to specify the values of the second tuning parameter. For instance, letting \code{what = "Theta"}, \code{penalty = "rho"} and \code{given = 1} then \code{plot.cglasso} returns the plot of the estimated partial correlation coefficients versus the \eqn{\rho}{rho}-values and given the first value of the used \eqn{\lambda}{lambda}-values. If \code{given} is left unspecified then a sequence of conditional profile plot is returned, a plot for each value of the second tuning parameter. 
}
Optionally, the user can specify a conditional profile plot using a model formula with the following template:

\sQuote{\code{what ~ penalty | given}}

For instance, the previous plot can be specified using the following model formula

\sQuote{\code{Theta ~ rho | 1}}

In this case, the arguments \code{penalty} and \code{given} must be left unspecified whereas \code{what} is used for the model formula (see also examples below).

The optional argument \code{GoF} can be used to identify the non-zero estimates by the output of a goodness-of-fit function, such as \code{\link{AIC.cglasso}} and \code{\link{BIC.cglasso}}. In this case, a vertical red dashed line is used the identify the optimal value of the tuning parameter reported on the \eqn{x}{x}-axis, whereas the path of the identified non-zero estimates is drawn using a solid black line; the remaining paths are drawn using gray dashed lines. This option is disabled if we let \sQuote{\code{GoF = NULL}}, when all the paths are drawn usign solid black lines.
}

\author{Luigi Augugliaro (\email{luigi.augugliaro@unipa.it})}

\references{
Augugliaro, L., Abbruzzo, A., and Vinciotti, V. (2020a) <\doi{10.1093/biostatistics/kxy043}>.
\eqn{\ell_1}{l1}-Penalized censored Gaussian graphical model.
\emph{Biostatistics} \bold{21}, e1--e16.

Augugliaro, L., Sottile, G., and Vinciotti, V. (2020b) <\doi{10.1007/s11222-020-09945-7}>.
The conditional censored graphical lasso estimator.
\emph{Statistics and Computing} \bold{30}, 1273--1289.

Friedman, J.H., Hastie, T., and Tibshirani, R. (2008) <\doi{10.1093/biostatistics/kxm045}>.
Sparse inverse covariance estimation with the graphical lasso.
\emph{Biostatistics} \bold{9}, 432--441.

Stadler, N. and Buhlmann, P. (2012) <\doi{10.1007/s11222-010-9219-7}>.
Missing values: sparse inverse covariance estimation and an extension to sparse regression.
\emph{Statistics and Computing} \bold{22}, 219--235.
}

\seealso{
\code{\link{cglasso}}, \code{\link{AIC.cglasso}} and \code{\link{BIC.cglasso}}.
}

\examples{
set.seed(123)

################################################################################
# Model 1: censored glasso estimator (Augugliaro \emph{and other}, 2020a)
# Y ~ N(0, Sigma) and probability of left/right censored values equal to 0.05
################################################################################
n <- 100L
p <- 10L
rho <- 0.3
Sigma <- outer(1L:p, 1L:p, function(i, j) rho^abs(i - j))
Z <- rcggm(n = n, Sigma = Sigma, probl = 0.05, probr = 0.05)
out <- cglasso(. ~ ., data = Z)

# coefficient profile plot of the off-diagonal entries of the precision matrix
plot(out, what = "Theta")
# the same as
#plot(out, what = "Theta", penalty = "rho")

# the output of a goodness-of-fit function can be used to identify the 
# non-zero estimates
plot(out, what = "Theta", GoF = BIC)

# letting 'GoF = NULL' the previous option is disabled
plot(out, what = "Theta", GoF = NULL)

# coefficient profile plot of the diagonal entries of the precision matrix
plot(out, what = "diag(Theta)")

# coefficient profile plot of the expected values of the response variables
plot(out, what = "b0")

################################################################################
# Model 2: conditional censored glasso estimator (Augugliaro \emph{and other}, 2020b)
# Y ~ N(b0 + XB, Sigma)  and probability of left/right censored values equal to 0.05
################################################################################
n <- 100L
p <- 10L
q <- 5L
b0 <- runif(p)
B <- matrix(runif(q * p), nrow = q, ncol = p)
X <- matrix(rnorm(n * q), nrow = n, ncol = q)
rho <- 0.3
Sigma <- outer(1L:p, 1L:p, function(i, j) rho^abs(i - j))
Z <- rcggm(n = n, b0 = b0, X = X, B = B, Sigma = Sigma, probl = 0.05, probr = 0.05)
out <- cglasso(. ~ ., data = Z)

# conditional profile plot of the estimater partial correlation coefficients versus
# the used values of the tunin parameter rho and given the first lambda-value
plot(out, what = "Theta", penalty = "rho", given = 1L)
out$lambda[1L]

lambda.id <- c(2, 4)
plot(out, what = "Theta", penalty = "rho", given = lambda.id)
out$lambda[lambda.id]

# in this case a sequence of ten conditional profile plots is returned, tha is a 
# plot for each lambda-value.
plot(out, what = "Theta", penalty = "rho")

# optionally, the user can specify the conditional profile plots using the model
# formula
plot(out, what = Theta ~ rho | lambda.id)
lambda.id

plot(out, what = Theta ~ rho)

# the output of a goodness-of-fit function can be used to identify the 
# non-zero estimates
plot(out, what = Theta ~ rho | 10, GoF = BIC)

# letting 'GoF = NULL' the previous option is disabled
plot(out, what = Theta ~ rho | 10, GoF = NULL)
}

\keyword{graphs}
