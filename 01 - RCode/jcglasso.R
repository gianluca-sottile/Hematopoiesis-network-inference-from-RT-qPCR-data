##### jcglasso function, used to fit the Joint Conditional Graphical Lasso Estimator with partially observed data ######

#' fit a GLM with lasso or elasticnet regularization
#' 
#' jcglasso fits the conditional graphical lasso model to datasets with censored and/or missing values across different conditions
#' by using the group lasso penalty.
#' 
#' \code{jcglasso} is the main model-fitting function and can be used to fit a broad range of extensions of the joint glasso estimator (Danaher \emph{and other}, 2011). It is specifically proposed to study datasets with censored and/or missing response values. To help the user, the \code{jcglasso} function has been designed to automatically select the most suitable extension by using the information stored in the \sQuote{\code{\link{datajcggm}}} object passed through \code{data}. 
#' 
#' @param formula an object of class \sQuote{\code{formula}} (or one that can be coerced to that class): a symbolic description of the model to be fitted
#' @param data an \R object of S3 class \sQuote{\code{datacggm}}, that is, the output of the function \sQuote{\code{datacggm}}.
#' @param subset an optional vector specifying a subset of observations to be used in the fitting process.
#' @param contrasts an optional list. See the \code{contrasts.arg} of \code{\link[stats]{model.matrix.default}}.
#' @param diagonal logical. Should diagonal entries of the concentration matrix be penalized? Default is \sQuote{\code{diagonal = FALSE}}.
#' @param weights.B an optional list of length \code{K} with \eqn{q\times p}{q x p} dimensional matrix of non-negative weights used to penalize the regression coefficients (intercepts are unpenalized). This matrix can be also used to specify the unpenalized regression coefficients (\sQuote{\code{weights.B[i, j] = 0}}) or the structural zeros (\sQuote{\code{weights.B[i, j] = +Inf}}). By default, all weights are set equal to 1, meaning that the penalized regression coefficients are unweighted.
#' @param weights.Tht an optional list of length \code{K} with symmetric matrix of non-negative weights used to penalize the partial regression coefficients. This matrix can be used to specify the unpenalized partial correlation coefficients (\sQuote{\code{weights.Tht[i, j] = 0}}) or the structural zeros in the precision matrix (\sQuote{\code{weights.Tht[i, j] = +Inf}}). By default, the off-diagonal entries of the matrix \sQuote{\code{weights.Tht}}  are set equal to 1, meaning that the partial correlation coefficients are unweighted, whereas the diagonal entries are equal to 0, that is the diagonal entries of the precision matrix are unpenalized.
#' @param penalty the penalty used across the \code{K} conditions. Default is \sQuote{\code{penalty = "group"}}. WARNING: fused lasso is not yet implemented.
#' @param nrho integer. The number of \eqn{\rho}{rho}-values used to penalize the partial correlation coefficients. By default \sQuote{\code{nrho = 10}}.
#' @param rho.min.ratio the smallest \eqn{\rho}{rho}-value is defined as a fraction of \sQuote{\code{rho.max}} (i.e., the smallest \eqn{\rho}{rho}-value for which all the estimated partial correlation coefficients are equal to zero). The default depends on the sample size \sQuote{\eqn{n}{n}} relative to the number of response variables \sQuote{\eqn{p}{p}}. If \sQuote{\eqn{p < n}{p < n}}, the default is \sQuote{1.0E-6} otherwise the value \sQuote{1.0E-2} is used as default. A very small value of \sQuote{\code{rho.min.ratio}} will lead to a saturated fitted model in the \sQuote{\eqn{p < n}{p < n}} case.
#' @param rho an optional user supplied decreasing sequence of \eqn{\rho}{rho}-values. By default \code{jcglasso} computes a sequence of \eqn{\rho}{rho}-values using \code{nrho} and \code{rho.min.ratio}. If \code{rho} is supplied then \code{nrho} and \code{rho.min.ratio} are overwritten.
#' @param nlambda integer. The number of \eqn{\lambda}{lambda}-values used to penalize the regression coefficients. By default \sQuote{\code{nlambda = 10}}.
#' @param lambda.min.ratio the smallest \eqn{\lambda}{lambda}-value is defined as a fraction of \sQuote{\code{lambda.max}} (i.e., the smallest \eqn{\lambda}{lambda}-value for which all the estimated regression coefficients are equal to zero). The default depends on the sample size \sQuote{\eqn{n}{n}} relative to the number of predictors \sQuote{\eqn{q}{q}}. If \sQuote{\eqn{q < n}{q < n}}, default is \sQuote{1.0E-6} otherwise the value \sQuote{1.0E-2} is used as default. A very small value of \sQuote{\code{lambda.min.ratio}} will lead to a saturated fitted model in the \sQuote{\eqn{q < n}{q < n}} case.
#' @param lambda an optional user-supplied decreasing sequence of \eqn{\lambda}{lambda}-values. By default \code{jcglasso} computes a sequence of \eqn{\lambda}{lambda}-values using \code{nlambda} and \code{lambda.min.ratio}. If \code{lambda} is supplied then \code{nlambda} and \code{lambda.min.ratio} are overwritten.
#' @param nu the value used to penalize the partial correlation coefficients of \code{X}, i.e., \eqn{\Omega}{Omega}
#' @param alpha1 the mixing parameter used in \eqn{\rho}{rho} to weigh between lasso and group-lasso penalty for \eqn{\Theta}{Theta}
#' @param alpha2 the mixing parameter used in \eqn{\lambda}{lambda} to weigh between lasso and group-lasso penalty for \code{B}
#' @param alpha3 the mixing parameter used in \eqn{\nu}{nu} to weigh between lasso and group-lasso penalty for \eqn{\Omega}{Omega}
#' @param maxit.em maximum number of iterations of the EM algorithm. Default is \code{1.0E+5}.
#' @param thr.em threshold for the convergence of the EM algorithm. Default value is \code{1.0E-4}.
#' @param maxit.bcd maximum number of iterations of the internal algorithms. Default is \code{1.0E+5}.
#' @param thr.bcd threshold for the convergence of the internal algorithms. Default is \code{1.0E-4}.
#' @param trace integer for printing information out as iterations proceed: \code{trace = 0} no information is printed out; \code{trace = 1} minimal information is printed on screen; \code{trace = 2} detailed information is printed on screen.
#' @param offset this can be used to specify an a priori known component to be included in the linear predictor during fitting. 
#' @param covar2corr logical. Experimental flag to convert precision matrices to correlation after fitting data.
#' @param truncate defaults to \code{1.0E-6}. At convergence, all values of theta below this number will be set to zero.
#' @return an object of S3 class \dQuote{\code{jcglasso}}, i.e., a named list containing the following components: \item{call}{the call that produced this object.} 
#' \item{Zipt}{an array of dimension \sQuote{\code{n x p x K x nlambda x nrho}} storing the \sQuote{working response matrices}, that is, \code{Zipt[, , , i, j]} is used as response array to compute the multilasso estimator (Sottile \emph{and other}, 2024).} 
#' \item{B}{an array of dimension \sQuote{\code{(q + 1) x p x K x nlambda x nrho}} storing the penalized estimate of the regression coefficient matrices. The accessor function \sQuote{\code{\link{coef.jcglasso}}} can be used to extract the desired estimates.} 
#' \item{mu}{an array of dimension \sQuote{\code{n x p x K x nlambda x nrho}} storing the fitted values. The accessor function \sQuote{\code{\link{fitted.jcglasso}}} can be used to extract the desired fitted values.} 
#' \item{R}{an array of dimension \sQuote{\code{n x p x K x nlambda x nrho}} storing the \sQuote{working residuals}, that is, \code{R[, , i, j]} is defined as the difference between 
#' \code{Zipt[, , , i, j]} and \code{mu[, , , i, j]}. The accessor function \sQuote{\code{\link{residuals.cglasso}}} can be used to extract the desidered residual matrix.} 
#' \item{S}{an array of dimension \sQuote{\code{(p + q) x (p + q) x K x nlambda x nrho}} storing the \sQuote{working empirical covariance matrices}, that is, \sQuote{\code{S[, , , i, j]}} is used to compute the joint glasso estimator.} 
#' \item{Sgm}{an array of dimension \sQuote{\code{(p + q) x (p + q) x K x nlambda x nrho}} storing the estimated covariance matrices. The accessor function \sQuote{\code{\link{coef}}} can be used to extract the desired estimates.} 
#' \item{Tht}{an array of dimension \sQuote{\code{(p + q) x (p + q) x K x nlambda x nrho}} storing the estimated precision matrices. The accessor funtion \sQuote{\code{\link{coef}}} can be used to extract the desired estimates.} 
#' \item{Omega}{an array of dimension \sQuote{\code{q x q x K x nlambda x nrho}} storing the estimated precision matrices. The accessor funtion \sQuote{\code{\link{coef}}} can be used to extract the desired estimates.} 
#' \item{dfB}{an array of dimension \sQuote{\code{(p + 1) x K x nlambda x nrho}} storing the number of estimated non-zero regression coefficients. Only for internal purpose.} 
#' \item{dfTht}{an array of dimension \sQuote{\code{K x nlambda x nrho}} storing the number of estimated non-zero (off-diagonal) partial correlation coefficients of \code{Y}. Only for internal purpose.}
#' \item{dfOmg}{an array of dimension \sQuote{\code{K x nlambda x nrho}} storing the number of estimated non-zero (off-diagonal) partial correlation coefficients of \code{X}. Only for internal purpose.} 
#' \item{InfoStructure}{a named list whose elements contain information about the estimated networks. Only for internal purpose.} \item{nit}{an array of dimension \sQuote{\code{2 x nlambda x nrho}} storing the number of EM steps.}
#' \item{Z}{the \sQuote{\code{datajcggm}} object used to compute the estimator.} \item{diagonal}{the flag used to specify if the diagonal entries of the precision matrix are penalized.} \item{weights.B}{the array of non-negative weights used for the regression coefficients.} 
#' \item{weights.Tht}{the array of non-negative weights used for the precision matrix.} \item{nlambda}{the number of \eqn{\lambda}{\lambda}-values used.} \item{lambda.min.ratio}{the value used to compute the smallest \eqn{\lambda}{lambda}-value.} 
#' \item{lambda}{the sequence of \eqn{\lambda}{lambda}-values used to fit the model.} \item{nrho}{the number of \eqn{\rho}{\rho}-values used.} \item{rho.min.ratio}{the value used to compute the smallest \eqn{\rho}{rho}-value.} \item{rho}{the sequence of \eqn{\rho}{rho}-values used to fit the model.}
#' \item{nu}{the value of \eqn{\nu}{nu} used to fit the model.} \item{alpha1}{the mixing parameter used in \eqn{\rho}{rho} to fit the model.} \item{alpha2}{the mixing parameter used in \eqn{\lambda}{lambda} to fit the model.} \item{alpha3}{the mixing parameter usedin in \eqn{\nu}{nu} to fit the model.}
#' \item{connected}{an array of dimension \sQuote{\code{(p + q) x nlambda x nrho}} storing the connected components for the joint glasso step.} \item{penalty}{the penalty used.} \item{model}{a description of the fitted model.} \item{maxit.em}{maximum number of iterations of the EM algorithm.}
#' \item{thr.em}{threshold for the convergence of the EM algorithm.} \item{maxit.bcd}{maximum number of iterations of the internal algorithms.} \item{thr.bcd}{threshold for the convergence of the internal algorithms.} \item{conv}{a description of the error that has occurred.} 
#' \item{subrout}{the name of the Fortran subroutine where the error has occurred (for internal debug only).} \item{trace}{the integer used for printing information on screen.} \item{nobs}{the sample size} \item{nresp}{the number of response variables used to fit the model.} 
#' \item{npred}{the number of predictors used to fit the model.}
#' 
#' @author Gianluca Sottile, Luigi Augugliaro\cr Maintainer: Gianluca Sottile \email{gianluca.sottile@@unipa.it}
#' 
#' @keywords models regression
#' 
#' @examples
#' 
#' n <- 200L                # sample size (integer)
#' p <- 25L                 # number of response variables (integer)
#' q <- 10L                 # number of predictors (integer)
#' K <- 2L                  # number of subpopulation (integer)
#' ncompPrecision <- 5L     # number of connected components in the precision matrices Omega e Theta for each k (integer)
#' ncompB <- 2L             # number of connected components in the coefficient regression matrix B for each k (integer)
#' percConnectedCompK <- .5 # percentages of connected components for Omega, Theta and B across k (numeric vector of length k - 1)
#' perc.NA <- 0.25          # percentage of missing-at-random data for each column
#' perc.X.na <- 0.25        # percentage of X columns with missing-at-random data
#' perc.Y.na <- 0.25        # percentage of Y columns with missing-at-random data
#' up <- 40.0               # upper censored value for the response matrix
#' perc.cens <- 0.25        # percentage of censored data for each column
#' perc.Y.cens <- 0.25      # percentage of Y columns with censored data
#' Thtmin <- 0.3            # minimum value of the precision matrix Theta
#' Thtmax <- 0.5            # maximum value of the precision matrix Theta
#' Omgmin <- 0.3            # minimum value of the precision matrix Omega
#' Omgmax <- 0.5            # maximum value of the precision matrix Omega
#' Bmin <- 0.5              # minimum value of coefficient regression matrix B
#' Bmax <- 0.7              # maximum value of coefficient regression matrix B
#' 
#' set.seed(1234)           # random seed for reproducibility
#' sim <- rjcggm(n = n, p = p, q = q, K = K, ncompPrecision = ncompPrecision, ncompB = ncompB, 
#'               percConnectedCompK = percConnectedCompK, perc.NA = perc.NA, perc.X.na = perc.X.na, perc.Y.na = perc.Y.na, 
#'               up = up, perc.cens = perc.cens, perc.Y.cens = perc.Y.cens, Thtmin = Thtmin, Thtmax = Thtmax, 
#'               Omgmin = Omgmin, Omgmax = Omgmax, Bmin = Bmin, Bmax = Bmax)
#' 
#' ##### fits the joint graphical lasso model to X datasets with censored and/or missing values              
#' model1 <- jcglasso(data = sim$ZX, nrho = 25L, rho.min.ratio = .01, alpha1 = .5, trace = 1L)
#' model1
#' 
#' ##### choice of the optimal nu parameter fixed alpha = 0.5
#' par(mfrow = c(2, 2))
#' # by using penalised estimates
#' plot(AIC(model1))
#' plot(BIC(model1))
#' 
#' # by using maximum likelihood estimates
#' qfun <- QFun2(model1, mle = TRUE)
#' plot(AIC(model1, Qfun = qfun))
#' plot(BIC(model1, Qfun = qfun))
#' 
#' ##### best model obtained minimizing the BIC criterion
#' model1best <- select.jcglasso(model1, GoF = BIC(model1, Qfun = qfun))
#' model1best
#' 
#' ##### simple grayscale plot of the precision matrices Omega
#' plot(model1best)
#' 
#' ##### fits the joint conditional graphical lasso model to the datasets with censored and/or missing values
#' ##### fixing nu to its optimal value and alpha3 = 0.5
#' model2 <- jcglasso(data = sim$Z, nrho = 5L, rho.min.ratio = .01, alpha1 = .5,
#'                    nlambda = 5L, lambda.min.ratio = .01, alpha2 = .5,
#'                    nu = model1best$rho, alpha3 = model1best$alpha1, trace = 1L)
#'                    
#' ##### choice of the optimal lambda and rho parameters fixed alpha1 = alpha2 = 0.5
#' ##### and fixing nu to its optimal value and alpha3 = 0.5
#' 
#' par(mfrow = c(2, 2)
#' # by using penalised estimates
#' plot(AIC(model2))
#' plot(BIC(model2))
#' 
#' # by using maximum likelihood estimates
#' qfun2 <- QFun2(model2, mle = TRUE)
#' plot(AIC(model2, Qfun = qfun2))
#' plot(BIC(model2, Qfun = qfun2))
#' 
#' ##### best model obtained minimizing the BIC criterion
#' model2best <- select.jcglasso(model2, GoF = BIC(model2, Qfun = qfun2))
#' model2best
#' 
#' ##### simple grayscale plot of the precision matrices Theta
#' plot(model2best)
#' 
#' #### perform post-hoc maximum likelihood refitting of a selected joint 
#' #### conditional graphical lasso model with censored and/or missing values.
#' model2bestMLE <- jcggm(model2best)
#' model2bestMLE
#' 
#' #### extract the main components of the model fit
#' coef(model2bestMLE, type = "B")
#' coef(model2bestMLE, type = "Theta")
#' coef(model2bestMLE, type = "Omega")
#' 
#' #### returns a named list of graphs using the results of an R object of class ‘jcglasso’ or ‘jcggm’.
#' graphs <- to_graph2(model2bestMLE)
#' graphs
#' 
#' par(mfrow = c(1, 2))
#' # plot of an undirected graph representing the conditional dependence structure among the p 
#' # response variables (i.e., Theta), highlighting the top 3 connected components
#' plot(graphs, type = "Gyy", highlight.connections = 3L)
#' # plot of a directed graph representing the effetcs of the q predictors onto the p 
#' # response variables (i.e., B), highlighting the top 3 connected components
#' plot(graphs, type = "Gxy", highlight.connections = 3L)
#' # plot of an undirected graph representing the conditional dependence structure among the q covariates (i.e., Omega)
#' plot(graphs, type = "Gxx")
#' # plot of both the undirected graph representing the conditional dependence structure among the p response variables (i.e.,Theta) 
#' # and the directed graph representing the effects of the q predictors on the p response variables (i.e., B).
#' plot(graphs, type = "conditional")
#' # overall plot of both the undirected graph representing the conditional dependence structure among the p response variables (i.e., Theta) 
#' # and the undirected graph representing the conditional dependence structure among the q covariates (i.e., Omega), 
#' # plus the directed graph representing the effects of the q predictors on the p response variables (i.e., B).
#' plot(graphs, type = "bipartite")
jcglasso <- function (formula, data, subset, contrasts = NULL, diagonal = FALSE, weights.B = NULL, 
                      weights.Tht = NULL, penalty = c("group", "fused"), 
                      nrho, rho.min.ratio, rho, nlambda, lambda.min.ratio, lambda, nu = NULL, 
                      alpha1 = 0.5, alpha2 = 0.5, alpha3 = 0.5,
                      maxit.em = 1E+5, thr.em = 1E-4, maxit.bcd = 1E+5, thr.bcd = 1E-04, 
                      trace = 0L, offset = NULL, covar2corr = FALSE, truncate = 1E-6) {
  this.call <- match.call()
  penalty <- match.arg(penalty)
  zero <- 1e-06
  
  # if(missing(rho1) | missing(rho2)) stop("Please, insert non-negative values for rho1 and rho2!!")
  # if(is.datajcggm(data)) data <- list(data)
  if (!is.datajcggm(data)) stop(sQuote("data"), " is not an object of class ", sQuote("datajcggm"))
  K <- length(data)
  
  X <- Y <- Z <- vector(mode = "list", length = K)
  for(k in seq_len(K)) {
    # testing 'formula'
    if (missing(formula)) formula <- . ~ . # stop(sQuote("formula"), " is missing")
    # testing 'data'
    # if (!is.datajcggm(data[[k]])) stop(sQuote("data"), " is not an object of class ", sQuote("datajcggm"))
    
    ##### formula2datacggm #####
    # testing LHS 'formula'
    fmlTerms <- function(fml, pos) paste(deparse(fml[[pos]], width.cutoff = 500L), collapse = " ")
    if (fmlTerms(formula, 2L) == ".")
      formula <- formula(paste0(paste0("cbind(", paste(colNames2(data)[[k]]$Y, collapse = ", "), ")"), " ~ ", fmlTerms(formula, 3L)))
    if (as.character(formula[[2L]])[1L] != "cbind")
      stop("Please use ", sQuote("cbind"), " to specify LHS in ", sQuote("formula"), " object")
    Y.LHS <- unlist(lapply(formula[[2L]][-1L], deparse))
    if (any(table(Y.LHS) > 1L))
      stop("repeated response variables are not permitted in ", sQuote("formula"), " object")
    noVars <- !is.element(Y.LHS, colNames2(data)[[k]]$Y)
    if (any(noVars)) stop("Following variables are not stored as rensponse variables: ", Y.LHS[noVars])
    if (is.null(getMatrix2(data, name = "X")[[k]])) {
      if (fmlTerms(formula, 3L) == ".") formula <- update(formula, . ~ 1)
      mt <- terms(formula)
      if (length(attr(mt, which = "term.labels")) != 0L)
        stop("Predictors are not stored in ", sQuote("data"))
      if (!as.logical(attr(mt, "intercept"))) {
        warning("Current version does not fit models without intercept term, thus it is added to the current formula.")
        formula <- update(formula, . ~ 1)
        mt <- terms(formula)
      }
      data.df <- data.frame(getMatrix2(data, name = "Y")[[k]][, Y.LHS, drop = FALSE])
      nResp <- length(Y.LHS)
      nPred <- 0L
    } 
    else {
      FML.RHS <- attr(terms(formula, data = data.frame(getMatrix2(data, name = "X")[[k]])), "term.labels")
      if (length(FML.RHS) > 0L) {
        FML.RHS <- paste(FML.RHS, collapse = " + ")
        formula <- formula(paste0(fmlTerms(formula, 2L), " ~ ", FML.RHS))
        X.RHS <- all.vars(formula[[3L]])
        noVars <- !is.element(X.RHS, colNames2(data)[[k]]$X)
        if (any(noVars)) stop("Following variables are not stored as predictors: ", X.RHS[noVars])
        data.df <- data.frame(getMatrix2(data, name = "Y")[[k]][, Y.LHS, drop = FALSE], 
                              data.frame(getMatrix2(data, name = "X")[[k]][, X.RHS, drop = FALSE]))
        nPred <- length(X.RHS)
      } else {
        data.df <- data.frame(getMatrix2(data, name = "Y")[[k]][, Y.LHS, drop = FALSE])
        nPred <- 0L
      }
      nResp <- length(Y.LHS)
    }
    mt <- terms(formula)
    if (!as.logical(attr(mt, "intercept"))) {
      warning("Current version does not fit models without intercept term, thus it has been added to the current formula.")
      formula <- update(formula, . ~ . + 1)
      mt <- terms(formula)
    }
    
    # creating model.frame #
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "subset"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf$na.action <- na.pass
    mf$drop.unused.levels <- TRUE
    mf$formula <- formula
    mf$data <- data.df
    mf[[1L]] <- quote(stats::model.frame)
    mf <- eval(mf, parent.frame())
    # updating 'datacggm' obejct
    Y[[k]] <- model.response(mf, type = "numeric")
    lo <- lower2(data)[[k]]$Y[colnames(Y[[k]])]
    up <- upper2(data)[[k]]$Y[colnames(Y[[k]])]
    if (length(attr(mt, which = "term.labels")) != 0L) {
      # Z[[k]] <- datajcggm(Y = Y[[k]], lo = lo, up = up)
    # } 
    # else {
      X[[k]] <- model.matrix(mt, mf, contrasts)[, -1L, drop = FALSE]
      lo <- c(lo, lower2(data)[[k]]$X[colnames(X[[k]])])
      up <- c(up, upper2(data)[[k]]$X[colnames(X[[k]])])
      # Z[[k]] <- datajcggm(Y = Y[[k]], X = X[[k]], lo = lo, up = up)
    }
  }
  X.null <- all(sapply(X, is.null))
  Z <- if(X.null) 
    datajcggm(Y = Y, lo = lo, up = up)
  else datajcggm(Y = Y, X = X, lo = lo, up = up)

  n <- nobs(Z)
  p <- nresp(Z)
  # if(all(p == p[1])) { p <- p[1] } else { stop("variabili diverse!??!") }
  q <- npred(Z)
  # if(all(q == q[1])) { q <- q[1] } else { stop("variabili diverse!??!") }
  if (p == 1L) stop("number of response variables is equal to ", sQuote(1))
  if(alpha1 < 0 | alpha1 > 1 | is.null(alpha1)) stop(sQuote("alpha1"), "must be in (0, 1)")
  if((alpha2 < 0 | alpha2 > 1 | is.null(alpha2)) & !X.null) stop(sQuote("alpha2"), "must be in (0, 1)")
  if((alpha3 < 0 | alpha3 > 1 | is.null(alpha3)) & !X.null) stop(sQuote("alpha3"), "must be in (0, 1)")
  
  xnames <- colNames2(Z)[[1]]$X
  ynames <- colNames2(Z)[[1]]$Y
  
  if (!is.vector(diagonal)) stop(sQuote("diagonal"), " is not a vector")
  if (length(diagonal) > 1L) stop(sQuote("diagonal"), " is not an object of length ", sQuote(1))
  if (!is.logical(diagonal)) stop(sQuote("diagonal"), " is not a logical object")
  if (!is.null(weights.B)) {
    for(k in seq_len(K)) {
      if (X.null) stop("Argument ", sQuote("weights.B"), " can not be used because ", sQuote("X"), " is missing")
      if (!is.matrix(weights.B[[k]])) stop(sQuote("weights.B"), " is not a matrix")
      weights.B.dim <- c(q, p)
      if (any(dim(weights.B[[k]]) != weights.B.dim)) 
        stop(sQuote("weights.B"), " is not a matrix with 'dim' attribute equal to ", 
             weights.B.dim[1L], " ", weights.B.dim[2L])
      if (any(weights.B[[k]] < 0)) 
        stop("negative weights in ", sQuote("weights.B "), " are not allowed")
      if (all(weights.B[[k]] == 0)) 
        stop("all entries in ", sQuote("weights.B "), " are equal to zero")
      weights.B[[k]][weights.B[[k]] == +Inf] <- .Machine$double.xmax
      rownames(weights.B[[k]]) <- xnames
      colnames(weights.B[[k]]) <- ynames
    }
  }
  else {
    if (!X.null) {
      weights.B <- vector(mode = "list", length = K)
      for(k in seq_len(K)) 
        weights.B[[k]] <- matrix(1, nrow = q, ncol = p, dimnames = list(xnames, ynames))
    }
  }
  if (!is.null(weights.Tht)) {
    for(k in seq_len(K)) {
      if (!is.matrix(weights.Tht[[k]])) stop(sQuote("weights.Tht"), " is not a matrix")
      weights.Tht.dim <- c(p, p)
      if (any(dim(weights.Tht[[k]]) != weights.Tht.dim))
        stop(sQuote("weights.Tht"), " is not a matrix with dim attribute equal to ",
             weights.Tht.dim[1L], " ", weights.Tht.dim[2L])
      if (any(weights.Tht[[k]] < 0))
        stop("negative weights in ", sQuote("weights.Tht "), " are not allowed")
      if (all(weights.Tht[[k]] == 0))
        stop("all entries in ", sQuote("weights.Tht "), " are equal to zero")
      weights.Tht[[k]][weights.Tht[[k]] == +Inf] <- .Machine$double.xmax
      if (!all(weights.Tht[[k]] == t(weights.Tht[[k]]))) stop(sQuote("weights.Tht"), " is not a symmetric matrix")
      if (!diagonal & any(diag(weights.Tht[[k]]) != 0))
        stop(sQuote("diagonal = FALSE"), " but some diagonal entry in ",
             sQuote("weights.Tht"), " is not zero")
      if (diagonal & all(diag(weights.Tht[[k]]) == 0))
        stop(sQuote("diagonal = TRUE"), " but all diagonal entries in ",
             sQuote("weights.Tht"), " are zero")
      rownames(weights.Tht[[k]]) <- ynames
      colnames(weights.Tht[[k]]) <- ynames
    }
  }
  else {
    weights.Tht <- vector(mode = "list", length = K)
    for(k in seq_len(K)) {
      weights.Tht[[k]] <- matrix(1, nrow = p, ncol = p)
      diag(weights.Tht[[k]]) <- ifelse(diagonal, 1, 0)
      rownames(weights.Tht[[k]]) <- ynames
      colnames(weights.Tht[[k]]) <- ynames
    }
  }
  
  if (!missing(nlambda) & !missing(lambda)) 
    warning("Argumnet ", sQuote("nlambda"), " is overwritten using length(lambda)")
  if (!missing(lambda.min.ratio) & !missing(lambda)) 
    warning("Argumnet ", sQuote("lambda.min.ratio"), " is overwritten using lambda")
  if (!missing(nlambda)) {
    if (X.null)
      stop("Argument ", sQuote("nlambda"), " can not be used because ", sQuote("X"), " is missing")
    if (!is.vector(nlambda)) 
      stop(sQuote("nlambda"), " is not a vector")
    if (length(nlambda) != 1) 
      stop(sQuote("nlambda"), " is not an object of length ", sQuote(1))
    if (abs(as.integer(nlambda) - nlambda) > 0) 
      stop(sQuote("nlambda"), " is not an integer value")
    if (nlambda <= 0) 
      stop(sQuote("nrnlambdaho1"), " is not a strictly positive integer value")
  }
  else nlambda <- ifelse(X.null, 1L, 10L)
  if (!missing(lambda.min.ratio)) {
    if (X.null)
      stop("Argument ", sQuote("lambda.min.ratio"), " can not be used because ", sQuote("X"), " is missing")
    if (!is.vector(lambda.min.ratio)) 
      stop(sQuote("lambda.min.ratio"), " is not a vector")
    if (length(lambda.min.ratio) != 1) 
      stop(sQuote("lambda.min.ratio"), " is not an object of length ", sQuote(1))
    if (lambda.min.ratio < 0 | lambda.min.ratio > 1) 
      stop(sQuote("lambda.min.ratio"), " does not belong to the closed interval [0, 1]")
  }
  else {
    if(X.null) lambda.min.ratio <- 0
    else lambda.min.ratio <- ifelse(all(p < n), zero, 0.01)
  }
  if (!missing(lambda)) {
    if (X.null)
      stop("Argument ", sQuote("lambda"), " can not be used because ", sQuote("X"), " is missing")
    if (!is.vector(lambda)) 
      stop(sQuote("lambda"), " is not a vector")
    if (any(lambda < 0)) 
      stop("some entry in ", sQuote("lambda"), " is negative")
    lambda <- sort(lambda, decreasing = TRUE)
    id <- lambda <= zero
    if (any(id)) 
      lambda <- c(lambda[!id], zero)
    nlambda <- length(lambda)
    lambda.min.ratio <- lambda[nlambda]/lambda[1L]
  }
  else {
    if(X.null) lambda <- 0
    else lambda <- double(nlambda)
  }
  
  if (!missing(nrho) & !missing(rho)) 
    warning("Argumnet ", sQuote("nrho"), " is overwritten using length(rho)")
  if (!missing(rho.min.ratio) & !missing(rho)) 
    warning("Argumnet ", sQuote("rho.min.ratio"), " is overwritten using rho")
  if (!missing(nrho)) {
    if (!is.vector(nrho)) 
      stop(sQuote("nrho"), " is not a vector")
    if (length(nrho) != 1) 
      stop(sQuote("nrho"), " is not an object of length ", sQuote(1))
    if (abs(as.integer(nrho) - nrho) > 0) 
      stop(sQuote("nrho"), " is not an integer value")
    if (nrho <= 0) 
      stop(sQuote("nrho"), " is not a strictly positive integer value")
  }
  else nrho <- 10L
  if (!missing(rho.min.ratio)) {
    if (!is.vector(rho.min.ratio)) 
      stop(sQuote("rho.min.ratio"), " is not a vector")
    if (length(rho.min.ratio) != 1) 
      stop(sQuote("rho.min.ratio"), " is not an object of length ", sQuote(1))
    if (rho.min.ratio < 0 | rho.min.ratio > 1) 
      stop(sQuote("rho.min.ratio"), " does not belong to the closed interval [0, 1]")
  }
  else rho.min.ratio <- ifelse(all(p < n), zero, 0.01)
  if (!missing(rho)) {
    if (!is.vector(rho)) 
      stop(sQuote("rho"), " is not a vector")
    if (any(rho < 0)) 
      stop("some entry in ", sQuote("rho"), " is negative")
    rho <- sort(rho, decreasing = TRUE)
    id <- rho <= zero
    if (any(id)) 
      rho <- c(rho[!id], zero)
    nrho <- length(rho)
    rho.min.ratio <- rho[nrho]/rho[1L]
  }
  else rho <- double(nrho)
  
  if (!is.vector(maxit.em)) stop(sQuote("maxit.em"), " is not a vector")
  if (length(maxit.em) != 1) 
    stop(sQuote("maxit.em"), " is not an object of length ", sQuote(1))
  if (abs(as.integer(maxit.em) - maxit.em) > 0) 
    stop(sQuote("maxit.em"), " is not an object of type ", dQuote("integer"))
  if (maxit.em <= 0) stop(sQuote("maxit.em"), " is not a positive integer")
  if (!is.vector(thr.em)) stop(sQuote("thr.em"), " is not a vector")
  if (length(thr.em) != 1) 
    stop(sQuote("thr.em"), " is not an object of length ", sQuote(1))
  if (thr.em <= 0) stop(sQuote("thr.em"), " is not a positive value")
  if (!is.vector(maxit.bcd)) stop(sQuote("maxit.bcd"), " is not a vector")
  if (length(maxit.bcd) != 1) 
    stop(sQuote("maxit.bcd"), " is not an object of length ", sQuote(1))
  if (abs(as.integer(maxit.bcd) - maxit.bcd) > 0) 
    stop(sQuote("maxit.bcd"), " is not an object of type ", dQuote("integer"))
  if (maxit.bcd <= 0) stop(sQuote("maxit.bcd"), " is not a positive integer")
  if (!is.vector(thr.bcd)) stop(sQuote("thr.bcd"), " is not a vector")
  if (length(thr.bcd) != 1) 
    stop(sQuote("thr.bcd"), " is not an object of length ", sQuote(1))
  if (thr.bcd <= 0) stop(sQuote("thr.bcd"), " is not a positive value")
  if (!is.vector(trace)) stop(sQuote("trace"), " is not a vector")
  if (length(trace) != 1) 
    stop(sQuote("trace"), " is not an object of length ", sQuote(1))
  if (is.logical(trace)) 
    stop(sQuote("trace"), " is not an object of type ", dQuote("integer"))
  if (abs(as.integer(trace) - trace) > 0) 
    stop(sQuote("trace"), " is not an object of type ", dQuote("integer"))
  if (!is.element(trace, c(0L, 1L, 2L))) 
    stop("not allowed value in ", sQuote("trace"), ". Please, choice ", 
         sQuote("0"), ", ", sQuote("1"), " or ", sQuote("2"))
  # if(alpha1 == 1) penalty <- "lasso"
  out <- jcglasso.fit(Z = Z, diagonal = diagonal, weights.B = weights.B, weights.Tht = weights.Tht, 
                      penalty = penalty, nrho = nrho, rho.min.ratio = rho.min.ratio, rho = rho, 
                      nlambda = nlambda, lambda.min.ratio = lambda.min.ratio, lambda = lambda, 
                      nu = nu, alpha1 = alpha1, alpha2 = alpha2, alpha3 = alpha3, 
                      maxit.em = maxit.em, thr.em = thr.em, maxit.bcd = maxit.bcd, 
                      thr.bcd = thr.bcd, trace = trace, offset = offset, covar2corr = covar2corr, truncate = truncate)
  if (out$conv != "Ok") {
    msg <- paste("jcglasso does not converge. Subroutine",
                 sQuote(out$subrout), "returns message:\n", sQuote(out$conv))
    warning(msg)
  }
  for(k in seq_len(K)) {
    if (!X.null) weights.B[[k]][weights.B[[k]] == .Machine$double.xmax] <- +Inf
    out$wTht[, , k][out$wTht[, , k] == .Machine$double.xmax] <- +Inf
  }
  if (!X.null) weights.B <- list2array(weights.B)
  InfoStructure <- list(Adj = out$Adj, 
                        Id = out$Id, np = out$nP, InfoP = out$InfoP,
                        ncomp = out$ncomp, Ck = out$Ck, pk = out$pk,
                        lo = out$lo, up = out$up, zm = out$zm, zv = out$zv,
                        id_X = out$id_X, id_Y = out$id_Y)
  if (all(sapply(seq_len(K), function(k) Z[[k]]$Info$Pattern[1L, "i"] > n[k])))
    model <- "joint glasso"
  else {
    id.mar <- any(sapply(seq_len(K), function(k) Z[[k]]$Info$Pattern[1L, "nmar"] > 0))
    id.cens <- any(sapply(seq_len(K), function(k) Z[[k]]$Info$Pattern[1L, c("nlc", "nrc")] > 0))
    if (id.mar & !id.cens)
      model <- "joint missglasso"
    if (!id.mar & id.cens)
      model <- "joint censored glasso"
    if (id.mar & id.cens)
      model <- "joint hybrid glasso"
  }
  if (q > 0L) model <- paste("conditional", model)
  
  out.jcglasso <- list(call = this.call, Zipt = out$Zipt, B = out$B,
                       mu = out$mu, R = out$R, S = out$S, Sgm = out$Sgm, Tht = out$Tht,
                       Omega = out$Omega, dfB = out$dfB, dfTht = out$dfTht, dfOmg = out$dfOmg,
                       InfoStructure = InfoStructure, nit = out$nit, Z = Z, diagonal = diagonal, 
                       weights.B = weights.B, weights.Tht = out$wTht,
                       nlambda = nlambda, lambda.min.ratio = lambda.min.ratio, lambda = out$lambda, 
                       nrho = nrho, rho.min.ratio = rho.min.ratio, rho = out$rho, 
                       nu = ifelse(is.null(nu), 0, nu), alpha1 = alpha1, alpha2 = alpha2, alpha3 = alpha3, 
                       connected = out$connected, penalty = penalty,
                       model = model, maxit.em = maxit.em, thr.em = thr.em,
                       maxit.bcd = maxit.bcd, thr.bcd = thr.bcd, conv = out$conv,
                       offset = out$offset, covar2corr = covar2corr, truncate = truncate,
                       subrout = out$subrout, trace = trace, nobs = n, nresp = p, npred = q)
  
  class(out.jcglasso) <- "jcglasso"
  out.jcglasso
}

jcglasso.fit <- function (Z, diagonal, weights.B, weights.Tht, penalty, 
                          nrho, rho.min.ratio, rho, nlambda, lambda.min.ratio, lambda, nu, 
                          alpha1, alpha2, alpha3,
                          maxit.em, thr.em, maxit.bcd, thr.bcd, trace, 
                          offset = NULL, covar2corr = FALSE, truncate = 1e-6) {
  
  storage.mode(diagonal) <- "integer"
  storage.mode(nlambda) <- "integer"
  storage.mode(lambda.min.ratio) <- "double"
  storage.mode(lambda) <- "double"
  storage.mode(nrho) <- "integer"
  storage.mode(rho.min.ratio) <- "double"
  storage.mode(rho) <- "double"
  storage.mode(maxit.em) <- "integer"
  storage.mode(thr.em) <- "double"
  storage.mode(maxit.bcd) <- "integer"
  storage.mode(thr.bcd) <- "double"
  storage.mode(trace) <- "integer"
  
  # if(is.datajcggm(Z)) Z <- list(Z)
  if (!is.datajcggm(Z)) stop(sQuote("data"), " is not an object of class ", sQuote("datajcggm"))
  K <- length(Z)
  storage.mode(K) <- "integer"
  if(K == 1) {alpha1 <- alpha2 <- alpha3 <- 1.0}
  storage.mode(alpha1) <- "double"
  storage.mode(alpha2) <- "double"
  storage.mode(alpha3) <- "double"
  
  n <- nobs(Z)
  storage.mode(n) <- "integer"
  weights <- proportions(n)
  storage.mode(weights) <- "double"
  
  p <- nresp(Z)
  # if(all(p == p[1])) { p <- p[1] } else { stop("variabili diverse!??!") }
  storage.mode(p) <- "integer"
  
  q <- npred(Z)
  # if(all(q == q[1])) { q <- q[1] } else { stop("variabili diverse!??!") }
  storage.mode(q) <- "integer"
  
  k_lab <- paste0("class", seq_len(K))
  rho_lab <- paste0("rho_", seq_len(nrho))
  lambda_lab <- paste0("lambda_", seq_len(nlambda))
  
  ynames <- colNames2(Z)[[1]]$Y
  xnames <- colNames2(Z)[[1]]$X
  znames <- c(ynames, xnames)
  
  X.null <- q == 0L
  
  seq_k <- seq_len(K)
  id_Y <- seq_len(p)
  id_X <- if(!X.null) seq_len(q) + p else NULL
  dim_Z <- p + q
  id_Z <- seq_len(dim_Z)
  storage.mode(dim_Z) <- "integer"
  
  weights.Tht <- list2array(weights.Tht)
  storage.mode(weights.Tht) <- "double"
  if(!X.null) {
    # weights.B2 <- lapply(seq_k, \(k) rbind(0.0, weights.B[[k]]))
    # weights.B2 <- list2array(weights.B2)
    # storage.mode(weights.B2) <- "double"
    weights.B <- list2array(weights.B)
    storage.mode(weights.B) <- "double"
  }
  
  ##### computing output and working objects #####
  Zmat <- array(0.0, dim = c(max(n), dim_Z, K),
                dimnames = list(NULL, znames, k_lab))
  zm <- zv <- loz <- upz <- matrix(0.0, nrow = dim_Z, ncol = K,
                                   dimnames = list(znames, k_lab))
  Id <- InfoP <- vector(mode = "list", length = K)
  nP <- integer(K)
  Zipt <- array(0.0, dim = c(max(n), dim_Z, K, nlambda, nrho),
                dimnames = list(NULL, response = znames, 
                                pop = k_lab, lambda = lambda_lab, rho = rho_lab))
  Zipt_n <- array(0.0, dim = c(max(n), dim_Z, K),
                  dimnames = list(NULL, znames, k_lab))
  Zipt_lo <- Zipt_up <- matrix(0.0, nrow = dim_Z, ncol = K,
                               dimnames = list(znames, k_lab))
  B <- array(0.0, dim = c(q + 1L, p, K, nlambda, nrho), 
             dimnames = list(coef = c("Int.", xnames), response = ynames, 
                             pop = k_lab, lambda = lambda_lab, rho = rho_lab))
  B_n <- array(0.0, dim = c(q + 1L, p, K), 
               dimnames = list(coef = c("Int.", xnames), response = ynames, pop = k_lab))
  mu <- R <- array(0.0, dim = c(max(n), dim_Z, K, nlambda, nrho), 
                   dimnames = list(NULL, response = znames, 
                                   pop = k_lab, lambda = lambda_lab, rho = rho_lab))
  mu_n <- mu[, , , 1L, 1L]
  dim(mu_n) <- c(max(n), dim_Z, K)
  R <- array(0.0, dim = c(max(n), p, K, nlambda, nrho), 
             dimnames = list(NULL, response = ynames, 
                             pop = k_lab, lambda = lambda_lab, rho = rho_lab))
  R_n <- R[, , , 1L, 1L]
  dim(R_n) <- c(max(n), p, K)
  Adj <- S <- Sgm <- Tht <- array(0.0, dim = c(dim_Z, dim_Z, K, nlambda, nrho), 
                                  dimnames = list(response = znames, response = znames, 
                                                  pop = k_lab, lambda = lambda_lab, rho = rho_lab))
  S_n <- S[, , , 1L, 1L]
  dim(S_n) <- c(dim_Z, dim_Z, K)
  Sgm_n <- Sgm[, , , 1L, 1L]
  dim(Sgm_n) <- c(dim_Z, dim_Z, K)
  Tht_n <- Tht[, , , 1L, 1L]
  dim(Tht_n) <- c(dim_Z, dim_Z, K)
  
  dfB <- array(0L, dim = c(p + 1L, K, nlambda, nrho), 
               dimnames = list(df = c(ynames, "Tot."), 
                               pop = k_lab, lambda = lambda_lab, rho = rho_lab))
  dfTht <- array(0L, dim = c(K, nlambda, nrho),
                 dimnames = list(pop = k_lab, lambda = lambda_lab, rho = rho_lab))
  
  if(!X.null) {
    Thtxx <- array(0.0, dim = c(q, q, K, nlambda, nrho), 
                   dimnames = list(response = xnames, response = xnames, 
                                   pop = k_lab, lambda = lambda_lab, rho = rho_lab))
    
    wThtxx <- attr(weights.Tht, "Thtxx")
    if(is.null(wThtxx)) {
      wThtxx <- array(1.0, dim = c(q, q, K),
                      dimnames = list(xnames, xnames, k_lab))
      for(k in seq_k) wThtxx[, , k] <- wThtxx[, , k] - diag(1.0, q, q)
    } else wThtxx <- list2array(wThtxx)
    weights.Tht <- list2array(lapply(seq_k, \(k) blockdiag(weights.Tht[, , k], matrix(wThtxx[, , k], q, q))))
    
    dfOmg <- array(0L, dim = c(K, nlambda, nrho),
                   dimnames = list(pop = k_lab, lambda = lambda_lab, rho = rho_lab))
    
  } 
  else {
    Thtxx <- NULL
    dfOmg <- NULL
  }
  
  ncomp <- matrix(0L, nrow = nlambda, ncol = nrho,
                  dimnames = list(lambda = lambda_lab, rho = rho_lab))
  pk <- Ck <- array(0L, dim = c(dim_Z, nlambda, nrho), 
                    dimnames = list(NULL, lambda = lambda_lab, rho = rho_lab))
  
  T1o <- zm
  T2o <- S_n
  T1 <- T1o[, 1]
  T2 <- T2o[, , 1]
  
  for(k in seq_k){
    Zmat[seq_len(n[k]), , k] <- cbind(getMatrix2(Z, name = "Y", ordered = TRUE)[[k]], 
                                      getMatrix2(Z, name = "X", ordered = TRUE)[[k]])
    Zmat[seq_len(n[k]), , k][is.na(Zmat[seq_len(n[k]), , k])] <- 0
    
    zm[, k] <- unlist(ColMeans2(Z)[[k]], use.names = FALSE)
    zv[, k] <- unlist(ColVars2(Z)[[k]], use.names = FALSE)
    loz[, k] <- Z[[k]]$Info$lo
    upz[, k] <- Z[[k]]$Info$up
    Id[[k]] <- event2(Z, ordered = TRUE)[[k]]
    InfoP[[k]] <- Z[[k]]$Info$Pattern
    storage.mode(InfoP[[k]]) <- "integer"
    nP[[k]] <- dim(InfoP[[k]])[1L]
  }
  
  nit.tot <- array(0L, dim = c(2L, nlambda, nrho), 
                   dimnames = list(steps = c("EM", "nit"), 
                                   lambda = lambda_lab, rho = rho_lab))
  
  cnnctd <- array(0L, dim = c(dim_Z, nlambda, nrho), 
                  dimnames = list(response = znames, 
                                  lambda = lambda_lab, rho = rho_lab))
  conv <- integer(1)
  subrout <- integer(1)
  
  is.null.offset <- is.null(offset)
  if(is.null.offset) {
    offset <- matrix(0.0, nrow = max(n), ncol = K)
    m_offset <- rep(0.0, K)
  }
  else {
    m_offset <- sapply(offset, mean)
    offset <- list2array(offset)
  }
  
  tmpSV <- starting_values(Zmat, zm, zv, loz, upz, Id, n, K, dim_Z, covar2corr,
                           offset, m_offset, id_X, id_Y, id_Z, X.null, Thtxx,
                           Zipt_n, Zipt_lo, Zipt_up, T1o, T2o,
                           B_n, mu_n, R_n, S_n, Sgm_n, Tht_n)
  Zipt_n <- tmpSV$Zipt_n
  Zipt_lo <- tmpSV$Zipt_lo
  Zipt_up <- tmpSV$Zipt_up
  T1o <- tmpSV$T1o
  T2o <- tmpSV$T2o
  B_n <- tmpSV$B_n
  mu_n <- tmpSV$mu_n
  R_n <- tmpSV$R_n
  S_n <- tmpSV$S_n
  Sgm_n <- tmpSV$Sgm_n
  Tht_n <- tmpSV$Tht_n
  Thtxx <- tmpSV$Thtxx
  
  if(all(rho == 0)) {
    U <- outer(id_Y, id_Y, "<")
    if(alpha1 < 1 & penalty == "group") 
      rho_max <- max(sqrt(rowSums(sapply(seq_k, function(k) 
        pmax(weights[k] * abs(S_n[id_Y, id_Y, k][U]), 0))^2)))
    else {
      rho_max <- max(sapply(seq_k, function(k) 
        max(abs(weights[k]*S_n[id_Y, id_Y, k][U]))))
    }
    rho_min <- rho_max * rho.min.ratio
    rho <- exp(seq(from = log(rho_max), to = log(rho_min), length = nrho))
    # rho <- seq(from = rho_max, to = rho_min, length = nrho)
  }
  if(is.null(nu)) {
    nu <- rho
  } else nu <- rep(nu, nrho)
  if(!X.null) {
    if(all(lambda == 0)) {
      if(alpha2 < 1 & penalty == "group") 
        # lambda_max <- max(sqrt(rowSums(sapply(seq_k, function(k)
        #   pmax(weights[k] * abs(S_n[id_X, id_Y, k] %*% Tht_n[id_Y, id_Y, k] / matrix(diag(Tht_n[id_Y, id_Y, k]), nrow = q, ncol = p, byrow = TRUE)), 0))^2)))
        lambda_max <- max(sqrt(rowSums(sapply(seq_k, function(k)
          pmax(weights[k] * abs(S_n[id_X, id_Y, k] %*% Tht_n[id_Y, id_Y, k]), 0))^2)))
      else {
        # lambda_max <- max(sapply(seq_k, function(k) 
        #   max(weights[k] * abs(S_n[id_X, id_Y, k] %*% Tht_n[id_Y, id_Y, k] / matrix(diag(Tht_n[id_Y, id_Y, k]), nrow = q, ncol = p, byrow = TRUE)))))
        lambda_max <- max(sapply(seq_k, function(k) 
          max(weights[k] * abs(S_n[id_X, id_Y, k] %*% Tht_n[id_Y, id_Y, k]))))
      }
      lambda_min <- lambda_max * lambda.min.ratio
      lambda <- exp(seq(from = log(lambda_max), to = log(lambda_min), length = nlambda))
      # lambda <- seq(from = lambda_max, to = lambda_min, length = nlambda)
    }
  } 
  
  X <- array(0.0, c(max(n), q, K))
  
  for(h in seq_len(nlambda)) {
    for(o in seq_len(nrho)) {
      nit <- double(2L)
      
      ##### starting EM algorithm #####
      if(trace == 2) {
        if(!X.null){
          cat("\n*************************************************************************\n",
              "Fitting jcglasso model number ", 
              formatC(o + nrho * (h - 1), digits = 0, width = 5, format = "d"),
              "\n\t\t      rho = ",
              formatC(rho[o], digits = 6, width = 9, format = "f"),
              "\n\t\t   alpha1 = ",
              formatC(alpha1, digits = 3, width = 8, format = "f"),
              "\n\t\t   lambda = ",
              formatC(lambda[h], digits = 6, width = 9, format = "f"),
              "\n\t\t   alpha2 = ",
              formatC(alpha2, digits = 3, width = 8, format = "f"),
              "\n\t\t       nu = ", 
              formatC(nu[o], digits = 6, width = 9, format = "f"),
              "\n\t\t   alpha3 = ",
              formatC(alpha3, digits = 3, width = 8, format = "f"),
              "\n")
        } 
        else {
          cat("\n*************************************************************************\n",
              "Fitting jcglasso model number ", 
              formatC(o + nrho * (h - 1), digits = 0, width = 5, format = "d"),
              "\n\t\t      rho = ",
              formatC(rho[o], digits = 6, width = 9, format = "f"),
              "\n\t\t   alpha1 = ",
              formatC(alpha1, digits = 3, width = 6, format = "f"),
              "\n")
        }
      }
      dOmg <- 0.0
      for(ii in 1:maxit.em) {
        B_o <- B_n
        Tht_o <- Tht_n
        if(!X.null) {
          Thtxx_o <- Thtxx[, , , h, o]
          dim(Thtxx_o) <- c(q, q, K)
        }
        
        ##### computing E step #####
        for(k in seq_k) {
          temp <- .Fortran(cglasso:::C_e_step_v2, n = n[k], p = dim_Z, Y = Zmat[seq_len(n[k]), , k], 
                           lo = loz[, k], up = upz[, k], nP = nP[k], 
                           InfoP = InfoP[[k]], T1o = T1o[, k], T2o = T2o[, , k], 
                           mu_n = mu_n[seq_len(n[k]), , k], Sgm = Sgm_n[, , k], Tht = Tht_n[, , k], 
                           Yipt_lo = Zipt_lo[, k], Yipt_up = Zipt_up[, k], Yipt_n = Zipt_n[seq_len(n[k]), , k], 
                           T1 = T1, T2 = T2, conv = conv)
          if(temp$conv != 0) { 
            subrout <- 1
            stop("error in E step!") 
          }
          if(trace == 2) cat("\n\tE-step completed!")
          
          zm[, k] <- temp$T1 / n[k]
          Zipt_n[seq_len(n[k]), , k] <- temp$Yipt_n
          
          S_n[, , k] <- temp$T2
          B_n[1L, , k] <- zm[id_Y, k] - m_offset[k]
          mu_n[seq_len(n[k]), id_Y, k] <- outer(offset[seq_len(n[k]), k], B_n[1L, , k], FUN = "+")
          R_n[seq_len(n[k]), , k] <- Zipt_n[seq_len(n[k]), id_Y, k] - mu_n[seq_len(n[k]), id_Y, k]
          
          if(!X.null) {
            mu_n[seq_len(n[k]), id_X, k] <- rep(zm[id_X, k], each = n[k])
            X[seq_len(n[k]), , k] <- Zipt_n[seq_len(n[k]), id_X, k] - mu_n[seq_len(n[k]), id_X, k]
          }
        }
        if(trace == 2) cat("\n\tE-step completed!")
        
        ##### fitting multilasso model #####
        if(!X.null) {
          if(trace == 2) {
            cat("\n\tM-step:\n\t       fitting multivariate sparse group lasso model with ",
                "\n\t       lambda = ", formatC(lambda[h], digits = 6, width = 9, format = "f"),
                "and alpha2 = ", formatC(alpha2, digits = 3, width = 4, format = "f"),
                "\n\n")
          }
          
          # browser()
          
          # tmpB <- apg(grad.b, prox.sparsegrouplasso,
          #             opts = list(p = p, q = q, n = rep(n, q), f = rep(weights, q),
          #                         A = do.call(blockdiag, array2list(X)),
          #                         b = do.call(blockdiag, array2list(R_n)),
          #                         Tht = do.call(blockdiag, array2list(Tht_n[id_Y, id_Y, , drop = FALSE])),
          #                         weights = do.call(blockdiag, array2list(weights.B)),
          #                         lambda = lambda[h], alpha = alpha2,
          #                         xm = zm[id_X, ], ym = zm[id_Y, ],
          #                         groups = list(idrow = rep(1:K, each = q), idcol = rep(1:K, each = p), k = K),
          #                         X_INIT = do.call(blockdiag, array2list(B_n[-1L, , , drop = FALSE])),
          #                         MAX_ITERS = maxit.bcd, EPS = thr.bcd, QUIET = ifelse(trace > 1L, 0L, 1L)))
          # dfB[, , h, o] <- tmpB$df
          # nit[2L] <- nit[2L] + tmpB$nit[1L]
          # B_n <- tmpB$x
          
          tmpB <- .Fortran(cglasso:::C_apg, p = p*K, q = q*K, n = max(n)*K, k = K,
                           nk = as.double(rep(n, each = q)), fk = as.double(rep(weights, each = q)),
                           A = do.call(blockdiag, array2list(X)), b = do.call(blockdiag, array2list(R_n)),
                           Tht = do.call(blockdiag, array2list(Tht_n[id_Y, id_Y, , drop = FALSE])),
                           weights = do.call(blockdiag, array2list(weights.B)),
                           lambda = lambda[h], alpha = alpha2, xm = zm[id_X, ], ym = zm[id_Y, ],
                           x = do.call(blockdiag, array2list(B_n[-1L, , , drop = FALSE])),
                           beta = B_n,
                           maxit = maxit.bcd, thr = thr.bcd, trace = trace,
                           df = matrix(0L, p + 1L, K), nit = integer(1), conv = integer(1))
          dfB[, , h, o] <- tmpB$df
          nit[2L] <- nit[2L] + tmpB$nit[1L]
          B_n <- tmpB$beta
        }
        for(k in seq_k){
          if(!X.null) {
            mu_n[seq_len(n[k]), id_Y, k] <- sweep(as.matrix(Zipt_n[seq_len(n[k]), id_X, k]) %*% B_n[-1L, , k], 2, B_n[1L, , k], "+")
            R_n[seq_len(n[k]), , k] <- Zipt_n[seq_len(n[k]), id_Y, k] - mu_n[seq_len(n[k]), id_Y, k] - offset[seq_len(n[k]), k]
          }
          YM <- crossprod(Zipt_n[seq_len(n[k]), , k], mu_n[seq_len(n[k]), , k])
          S_n[, , k] <- (S_n[, , k] + crossprod(mu_n[seq_len(n[k]), , k]) - YM - t(YM)) / n[k]
          if(covar2corr) S_n[, , k] <- cov2cor(S_n[, , k])
        }
        
        ##### computing theta and S values for connected components #####
        if(trace == 2) {
          if(X.null) {
            cat("\n")
            cat("\tM-step:\n\t       fitting joint glasso model with",
                "rho = ", formatC(rho[o], digits = 6, width = 9, format = "f"),
                "and alpha1 = ", formatC(alpha1, digits = 3, width = 4, format = "f"),
                "\n\n")
          } 
          else {
            cat("\n\t       fitting joint glasso model with",
                "\n\t       rho = ", formatC(rho[o], digits = 6, width = 9, format = "f"),
                "and alpha1 = ", formatC(alpha1, digits = 3, width = 4, format = "f"),
                "\n\t        nu = ", formatC(nu[o], digits = 6, width = 9, format = "f"), 
                "and alpha3 = ", formatC(alpha3, digits = 3, width = 4, format = "f"),
                "\n\n")
          }
        }
        
        tmpTht <- .Fortran(cglasso:::C_admm_tht_sub, p = p, N = K, fk = weights, S = S_n[id_Y, id_Y, , drop = FALSE],
                           wTht = weights.Tht[id_Y, id_Y, , drop = FALSE], pendiag = diagonal, rho = rho[o], alpha = alpha1,
                           maxit = maxit.bcd, thr = thr.bcd, Tht = Tht_n[id_Y, id_Y, , drop = FALSE], k = integer(1),
                           Ck = integer(p), pk = integer(p), nit = integer(1), df = integer(K),
                           conv = integer(1), trace = trace)
        nit[2] <- nit[2] + tmpTht$nit
        Tht_n[id_Y, id_Y, ] <- tmpTht$Tht
        
        if(!X.null) { 
          tmpOmg <- .Fortran(cglasso:::C_admm_tht_sub, p = q, N = K, fk = weights, S = S_n[id_X, id_X, , drop = FALSE], 
                             wTht = weights.Tht[id_X, id_X, , drop = FALSE], pendiag = diagonal, rho = nu[o], alpha = alpha3, 
                             maxit = maxit.bcd, thr = thr.bcd, Tht = Tht_n[id_X, id_X, , drop = FALSE], k = integer(1),
                             Ck = integer(q), pk = integer(q), nit = integer(1), df = integer(K), 
                             conv = integer(1), trace = trace)
          nit[2] <- nit[2] + tmpOmg$nit
          Thtxx[, , , h, o] <- tmpOmg$Tht
        }
        
        nit[1] <- ii
        
        for(k in seq_k) {
          Adj[id_Y, id_Y, k, h, o] <- 1*(Tht_n[id_Y, id_Y, k] != 0)
          if(!X.null) { 
            Adj[id_X, id_X, k, h, o] <- 1*(Thtxx[, , k, h, o] != 0)
            Tht_n[id_X, id_X, k] <- Thtxx[, , k, h, o] + B_n[-1L, , k] %*% Tht_n[id_Y, id_Y, k] %*% t(matrix(B_n[-1L, , k], q, p))
            Tht_n[id_X, id_Y, k] <- -(B_n[-1L, , k] %*% Tht_n[id_Y, id_Y, k])
            Tht_n[id_Y, id_X, k] <- t(Tht_n[id_X, id_Y, k])
          }
          Sgm_n[, , k] <- solve(Tht_n[, , k])
        } 
        
        dB <- if(!X.null)
          sum(weights * sapply(seq_len(K), function(k) norm(B_n[, , k] - B_o[, , k], type = "F") / (p + dfB[p + 1, k, h, o])))
        else sum(weights * sapply(seq_len(K), function(k) norm(B_n[, , k] - B_o[, , k], type = "2") / (p + dfB[p + 1, k, h, o])))
        dTht <- sum(weights * sapply(seq_len(K), function(k) norm(Tht_n[id_Y, id_Y, k] - Tht_o[id_Y, id_Y, k], type = "F") / (p + tmpTht$df[k])))
        if(!X.null) dOmg <- sum(weights * sapply(seq_len(K), function(k) norm(matrix(Thtxx[, , k, h, o], q, q) - matrix(Thtxx_o[, , k], q, q), "F") / (q + tmpOmg$df[k])))
        
        if(trace == 2) {
          cat("\tM-step completed!\n")
          if(!X.null) {
            cat("\n\tChecking convergence criterion (threshold = ", 
                formatC(thr.em, digits = 6, width = 8, format = "f"),
                ")\n\t||B_old - B_new|_F / dfB  = ",
                formatC(dB, digits = 7, width = 9, format = "f"),
                "\n\t||Tht_old - Tht_new||_F / dfTht = ",
                formatC(dTht, digits = 7, width = 9, format = "f"), 
                "\n\t||Omg_old - Omg_new||_F / dfOmg = ", 
                formatC(dOmg, digits = 7, width = 9, format = "f"), 
                "\n")
            if(max(dB, dTht, dOmg) <= thr.em) cat("\tConvergence criterion is met!\n\n")
          }
          else {
            cat("\n\tChecking convergence criterion (threshold = ", 
                formatC(thr.em, digits = 6, width = 8, format = "f"),
                ")\n\t||B_old - B_new|_F / dfB  = ",
                formatC(dB, digits = 7, width = 9, format = "f"),
                "\n\t||Tht_old - Tht_new||_F / dfTht = ",
                formatC(dTht, digits = 7, width = 9, format = "f"), 
                "\n")
            if(max(dB, dTht) <= thr.em) cat("\tConvergence criterion is met!\n\n")
          }
        } 
        if(max(dB, dTht, dOmg) <= thr.em) break
      }
      if(ii >= maxit.em) {
        conv <- 1
        subrout <- 3
        stop("error in M step (maximum number of iteration attained)!")
      }
      if(trace == 1) {
        if(!X.null){
          cat("\njcglasso model number ", 
              formatC(o + nrho * (h - 1), digits = 0, width = 5, format = "d"),
              "\n\t        rho = ",
              formatC(rho[o], digits = 6, width = 9, format = "f"),
              "\n\t     alpha1 = ",
              formatC(alpha1, digits = 3, width = 6, format = "f"),
              "\n\t     lambda = ",
              formatC(lambda[h], digits = 6, width = 9, format = "f"),
              "\n\t     alpha2 = ",
              formatC(alpha2, digits = 3, width = 6, format = "f"),
              "\n\t         nu = ",
              formatC(nu[o], digits = 6, width = 9, format = "f"),
              "\n\t     alpha3 = ",
              formatC(alpha3, digits = 3, width = 6, format = "f"),
              "\n\t   EM-steps = ",
              formatC(nit[1], digits = 0, width = 5, format = "d"),
              "\n\t      steps = ",
              formatC(nit[2], digits = 0, width = 5, format = "d"),
              "\n"
          )
        } 
        else {
          cat("\njcglasso model number ", 
              formatC(o + nrho * (h - 1), digits = 0, width = 5, format = "d"),
              "\n\t        rho = ",
              formatC(rho[o], digits = 6, width = 9, format = "f"),
              "\n\t     alpha1 = ",
              formatC(alpha1, digits = 3, width = 6, format = "f"),
              "\n\t   EM-steps = ",
              formatC(nit[1], digits = 0, width = 5, format = "d"),
              "\n\t      steps = ",
              formatC(nit[2], digits = 0, width = 5, format = "d"),
              "\n"
          )
        }
      }
      
      ##### buiding output components #####
      ncomp[h, o] <- tmpTht$k #+ tmpOmg$k
      Ck[id_Y, h, o] <- tmpTht$Ck
      pk[id_Y, h, o] <- tmpTht$pk
      cnnctd[id_Y, h, o] <- tmpTht$Ck
      if(!X.null) {
        Ck[id_X, h, o] <- tmpOmg$Ck
        pk[id_X, h, o] <- tmpOmg$pk
        cnnctd[id_X, h, o] <- tmpOmg$Ck
      }
      for(k in seq_k) {
        Zipt[, , k, h, o] <- Zipt_n[, , k]
        B[, , k, h, o] <- B_n[, , k]
        mu[, , k, h, o] <- mu_n[, , k]
        R[, , k, h, o] <- R_n[, , k]
        S[, , k, h, o] <- S_n[, , k]
        Sgm[, , k, h, o] <- Sgm_n[, , k]
        Tht[, , k, h, o] <- Tht_n[, , k]
        diag(Adj[, , k, h, o]) <- 0
        dfTht[k, h, o] <- tmpTht$df[k]
        if(!X.null) {
          Adj[id_X, id_Y, k, h, o] <- 1*(B[-1, , k, h, o] != 0)
          Adj[id_Y, id_X, k, h, o] <- t(Adj[id_X, id_Y, k, h, o])
          dfOmg[k, h, o] <- dfOmg[k, h, o] + tmpOmg$df[k]
        }
        
        row.order <- order(Z[[k]]$Info$order)
        Zipt[seq_len(n[k]), , k, h, o] <- Zipt[seq_len(n[k]), , k, h, o][row.order, ,drop = FALSE]
        mu[seq_len(n[k]), , k, h, o] <- mu[seq_len(n[k]), , k, h, o][row.order, , drop = FALSE]
        R[seq_len(n[k]), , k, h, o] <- R[seq_len(n[k]), , k, h, o][row.order, , drop = FALSE]
      }
      nit.tot[, h, o] <- nit
    }
  }
  
  names(Z) <- k_lab
  
  out <- list(n = n, p = p, q = q, id_X = id_X, id_Y = id_Y, Z = Zmat, 
              Id = Id, nP = nP, InfoP = InfoP, ncomp = ncomp, Ck = Ck, pk = pk, 
              lo = loz, up = upz, zm = zm, zv = zv, wTht = weights.Tht, 
              nrho = nrho, rhoratio = rho.min.ratio, rho = rho,
              nlambda = nlambda, lambdaratio = lambda.min.ratio, lambda = lambda, 
              alpha1 = alpha1, alpha2 = alpha2, alpha3 = alpha3, pendiag = diagonal, connected = cnnctd, 
              maxit_em = maxit.em, thr_em = thr.em, maxit_bcd = maxit.bcd, thr_bcd = thr.bcd, 
              Zipt = Zipt, B = B, mu = mu, R = R, S = S, Sgm = Sgm, Tht = Tht, 
              Adj = Adj, Omega = Thtxx, dfB = dfB, dfTht = dfTht, dfOmg = dfOmg,
              nit = nit.tot, conv = conv, offset = offset,
              subrout = subrout, trace = trace)
  
  out$conv <- switch(as.character(out$conv),
                     "-1" = "memory allocation error",
                     "0" = "Ok",
                     "1" = "maximum number of iterations has been exceeded",
                     "2" = "error in E-step",
                     "3" = "matrix inversion failed")
  
  return(out)
}

##### auxiliary function used by jcglasso algorithm #####
starting_values <- function(Zmat, zm, zv, loz, upz, Id, n, K, dim_Z, covar2corr,
                            offset, m_offset, id_X, id_Y, id_Z, X.null, Thtxx,
                            Zipt_n, Zipt_lo, Zipt_up, T1o, T2o,
                            B_n, mu_n, R_n, S_n, Sgm_n, Tht_n) {
  seq_k <- seq_len(K)
  ##### computing Zipt_n, Zipt_lo and Zipt_up #####
  Zipt_lo <- zm - 3 * sqrt(zv)
  Zipt_up <- zm + 3 * sqrt(zv)
  Zipt_n <- Zmat
  
  for(k in seq_k) {
    ##### computing starting values #####
    for(i in id_Z) {
      id <- Id[[k]][, i]
      T1o[i, k] <- sum(Zipt_n[seq_len(n[k]), i, k][id == 0L])
      
      if(any(id == 9L)) Zipt_n[seq_len(n[k]), i, k][id == 9L] <- zm[i, k]
      if(any(id == -1L)) {
        z <- (loz[i, k] - zm[i, k]) / sqrt(zv[i, k])
        tmean <- zm[i, k] - sqrt(zv[i, k]) * dnorm(z) / pnorm(z)
        tmean <- ifelse(tmean <= Zipt_lo[i, k], Zipt_lo[i, k], tmean)
        Zipt_n[seq_len(n[k]), i, k][id == -1L] <- tmean
      }
      if(any(id == 1L)) {
        z <- (upz[i, k] - zm[i, k]) / sqrt(zv[i, k])
        tmean <- zm[i, k] + sqrt(zv[i, k]) * dnorm(z) / pnorm(z)
        tmean <- ifelse(tmean >= Zipt_up[i, k], Zipt_up[i, k], tmean)
        Zipt_n[seq_len(n[k]), i, k][id == 1L] <- tmean
      }
      
      for(j in i:dim_Z) {
        if(any(id <- Id[[k]][, i] == 0L & Id[[k]][, j] == 0L)) {
          T2o[i, j, k] <- sum(Zipt_n[id, i, k] * Zipt_n[id, j, k])
          T2o[j, i, k] <- T2o[i, j, k]
        }
      }
    }
    
    zm[, k] <- colMeans(Zipt_n[seq_len(n[k]), , k])
    B_n[1L, , k] <- zm[id_Y, k] - m_offset[k]
    
    ##### computing mu_n, R_n and S_n, Tht_n #####
    mu_n[seq_len(n[k]), id_Y, k] <- outer(offset[seq_len(n[k]), k], B_n[1L, , k], FUN = "+")
    R_n[seq_len(n[k]), , k] <- Zipt_n[seq_len(n[k]), id_Y, k] - mu_n[seq_len(n[k]), id_Y, k]
    S_n[id_Y, id_Y, k] <- crossprod(R_n[seq_len(n[k]), , k]) / n[k]
    if(!X.null) {
      mu_n[seq_len(n[k]), id_X, k] <- rep(zm[id_X, k], each = n[k])
      S_n[id_X, id_X, k] <- crossprod(Zipt_n[seq_len(n[k]), id_X, k] - mu_n[seq_len(n[k]), id_X, k]) / n[k]
      S_n[id_X, id_Y, k] <- crossprod(Zipt_n[seq_len(n[k]), id_X, k], R_n[seq_len(n[k]), , k]) / n[k]
      S_n[id_Y, id_X, k] <- t(S_n[id_X, id_Y, k])
      Thtxx[, , k, 1, 1] <- diag(1 / diag(as.matrix(S_n[id_X, id_X, k])))
    }
    if(covar2corr) S_n[, , k] <- cov2cor(S_n[, , k]) 
    Sgm_n[, , k] <- diag(diag(S_n[, , k]))
    Tht_n[, , k] <- diag(1 / diag(S_n[, , k]))
  }
  
  list(Zipt_n = Zipt_n, Zipt_lo = Zipt_lo, Zipt_up = Zipt_up, T1o = T1o, 
       T2o = T2o, B_n = B_n, mu_n = mu_n, R_n = R_n, S_n = S_n, 
       Sgm_n = Sgm_n, Tht_n = Tht_n, Thtxx = Thtxx)
}

list2array <- function(x) {
  if(!is.list(x)) stop("x is not a list!")
  K <- length(x)
  if(is.matrix(x[[1]])) {
    n <- sapply(x, nrow)
    p <- ncol(x[[1]])
    y <- array(0, dim = c(max(n), p, K), dimnames = list(NULL, colnames(x[[1]]), names(x)))
    for(k in seq_len(K)){ y[seq_len(n[k]), , k] <- x[[k]] }
  } else {
    n <- length(x[[1]])
    y <- matrix(unlist(x), nrow = n, ncol = K, dimnames = list(names(x[[1]]), names(x)))
  }
  y
}
array2list <- function(x) {
  if(!is.array(x)) stop("x is not a list!")
  dimX <- dim(x)
  K <- tail(dimX, 1)
  y <- vector(mode = "list", length = K)
  for(k in seq_len(K)) y[[k]] <- if(length(dimX) == 2) x[, k] else na.omit(matrix(x[,,k], dimX[1], dimX[2]))
  names(y) <- tail(dimnames(x), 1)[[1]]
  y
}
blockdiag <- function(...) {
  args <- list(...)
  if(length(args) == 1L) {
    ret <- args[[1L]]
  }
  else {
    nc <- sapply(args, ncol)
    cumnc <- cumsum(nc)
    NC <- sum(nc)
    rowfun <- function(m, zbefore, zafter) {
      cbind(matrix(0, ncol = zbefore, nrow = nrow(m)), 
            m, matrix(0, ncol = zafter, nrow = nrow(m)))
    }
    ret <- rowfun(args[[1]], 0, NC - ncol(args[[1]]))
    for (i in 2:length(args)) {
      ret <- rbind(ret, rowfun(args[[i]], cumnc[i - 1], NC - cumnc[i]))
    }
  }
  ret
}

apg <- function(grad_f, prox_h, opts) {
  
  # Set default parameters
  X_INIT <- if(is.null(opts[["X_INIT"]])) 
    matrix(0.0, ncol(opts[["A"]]), ncol(opts[["b"]]))
  else opts[["X_INIT"]]
  MAX_ITERS <- if(is.null(opts[["MAX_ITERS"]])) # maximum iterations before termination
    1E5
  else opts[["MAX_ITERS"]] 
  EPS <- if(is.null(opts[["EPS"]])) # tolerance for termination
    1E-6
  else opts[["EPS"]] 
  QUIET <- if(is.null(opts[["QUIET"]])) # if false writes out information every 100 iters
    TRUE 
  else opts[["QUIET"]]
  USE_RESTART <- TRUE # use adaptive restart scheme
  ALPHA <- 1.01 # step-size growth factor
  BETA <- 0.5 # step-size shrinkage factor
  GEN_PLOTS <- TRUE # if true generates plots of norm of proximal gradient
  USE_GRA <- FALSE # if true uses UN-accelerated proximal gradient descent (typically slower)
  STEP_SIZE <- NULL # starting step-size estimate, if not set then apg makes initial guess
  FIXED_STEP_SIZE <- FALSE # don't change step-size (forward or back tracking), uses initial step-size throughout, only useful if good STEP_SIZE set
  
  # # Replace the default parameters by the ones provided in opts if any
  # for (u in c("X_INIT", "USE_RESTART", "MAX_ITERS", "EPS", "ALPHA", "BETA", "QUIET", "GEN_PLOTS", 
  #             "USE_GRA", "STEP_SIZE", "FIXED_STEP_SIZE")) {
  #   eval(parse(text = paste('if (exists("', u, '", where = opts)) ', u, ' <- opts[["', u, '"]]', sep = '')))
  # }
  
  # Initialization
  x <- X_INIT
  y <- x
  g <- grad_f(y, opts)
  if ((nrm_g <- norm_vec(g)) < EPS) return(list(x = x, t = 0))
  
  theta <- 1
  
  # Initial step size
  if (is.null(STEP_SIZE)) {
    # Barzilai-Borwein step-size initialization:
    t <- 1 / nrm_g
    x_hat <- x - t*g
    g_hat <- grad_f(x_hat, opts)
    t <- abs(sum((x - x_hat)*(g - g_hat)) / sum((g - g_hat)^2))
  } 
  else {
    t <- STEP_SIZE
  }
  
  # Main loop
  for (k in seq_len(MAX_ITERS)) {
    # if (!QUIET && (k %% 10 == 0)) message(paste0('iter num ', k, ', norm(tGk): ', err1, ', step-size: ', t))
    
    x_old <- x
    y_old <- y
    
    # The proximal gradient step (update x)
    x <- prox_h(y - t * g, t, opts)
    
    # The error for the stopping criterion
    err1 <- norm_vec(y - x) / max(1, norm_vec(x))
    if (err1 < EPS) break
    
    # Update theta for acceleration
    theta <- if(!USE_GRA) 2 / (1 + sqrt(1 + 4 / (theta^2))) else 1
    
    # Update y
    if (USE_RESTART && sum((y - x) * (x - x_old)) > 0) {
      x <- x_old
      y <- x
      theta <- 1
    } 
    else {
      y <- x + (1 - theta) * (x - x_old)
    }
    
    # New gradient
    g_old <- g
    g <- grad_f(y, opts)
    
    # Update stepsize by TFOCS-style backtracking
    if (!FIXED_STEP_SIZE) {
      t_hat <- 0.5 * sum((y - y_old)^2) / abs(sum((y - y_old) * (g_old - g)))
      t <- min(ALPHA * t, max(BETA * t, t_hat))
    }
    if (!QUIET) message(paste0('iter num ', k, ', norm(tGk): ', err1, ', step-size: ', t))
  }
  if (!QUIET) {
    message(paste('iter num ', k,', norm(tGk): ', err1, ', step-size: ', t, sep = ""))
    if (k == MAX_ITERS) message(paste('Warning: maximum number of iterations reached'))
    message('Terminated')
  }
  
  df <- matrix(0L, opts[["p"]] + 1L, opts[["groups"]]$k)
  out <- array(0.0, c(opts[["q"]] + 1L, opts[["p"]], opts[["groups"]]$k))
  for(i in 1:opts[["groups"]]$k) {
    out[1L, , i] <- opts[["ym"]][, i] - opts[["xm"]][, i] %*% x[opts[["groups"]]$idrow == i, opts[["groups"]]$idcol == i]
    out[-1L, , i] <- x[opts[["groups"]]$idrow == i, opts[["groups"]]$idcol == i]
    df[1:opts[["p"]], i] <- colSums(out[-1L, , i] != 0)
    df[opts[["p"]] + 1L, i] <- sum(df[, i])
  }
  
  # Return solution, degrees of freedom, step size and number of iterations
  return(list(x = out, df = df, t = t, nit = k))
}
norm_vec <- function(x) sqrt(sum(x^2))
prox.lasso <- function(x, t, opts) {
  sign(x) * pmax(abs(x) - t * opts[["lambda"]] * opts[["weights"]], 0)
}
prox.grouplasso <- function(x, t, opts) {
  nrm <- 0
  for (i in 1:opts[["groups"]]$k) nrm <- nrm + x[opts[["groups"]]$idrow == i, opts[["groups"]]$idcol == i]^2
  nrm <- kronecker(diag(opts[["groups"]]$k), sqrt(nrm))
  notshrunk <- (nrm > (t * opts[["lambda"]] * opts[["weights"]])) * 1
  nrm <- nrm + (1 - notshrunk)
  nrm <- (1 - (t * opts[["lambda"]] * opts[["weights"]]) / nrm) * notshrunk
  x * nrm
}
prox.sparsegrouplasso <- function(x, t, opts) {
  x <- sign(x) * pmax(abs(x) - t * opts[["alpha"]] * opts[["lambda"]] * opts[["weights"]], 0)
  nrm <- 0
  for (i in 1:opts[["groups"]]$k) nrm <- nrm + x[opts[["groups"]]$idrow == i, opts[["groups"]]$idcol == i]^2
  nrm <- kronecker(diag(opts[["groups"]]$k), sqrt(nrm))
  notshrunk <- (nrm > (t * (1 - opts[["alpha"]]) * opts[["lambda"]] * opts[["weights"]])) * 1
  nrm <- nrm + (1 - notshrunk)
  nrm <- (1 - (t * (1 - opts[["alpha"]]) * opts[["lambda"]] * opts[["weights"]]) / nrm) * notshrunk
  x * nrm
}
grad.b <- function(x, opts) {
  MU <- opts[["A"]] %*% x
  R <- opts[["b"]] - MU
  - opts[["f"]] * crossprod(opts[["A"]], R) %*% opts[["Tht"]] / opts[["n"]]
}

gradB <- function(obj, lambda.id, rho.id) {
  K <- length(obj$Z)
  
  n <- obj$nobs
  lambda <- obj$lambda[lambda.id]
  id_X <- obj$InfoStructure$id_X
  id_Y <- obj$InfoStructure$id_Y
  
  if(K == 1) {
    X <- obj$Zipt[, id_X, K, lambda.id, rho.id]
    R <- obj$R[, , K, lambda.id, rho.id]
    Tht <- coef(obj, type = "Theta", lambda.id = lambda.id, rho.id = rho.id)
    grdB <- t(X) %*% R %*% Tht / n
    
    B <- coef(obj, type = "B", lambda.id = lambda.id, rho.id = rho.id)
    eps <- (grdB - lambda * sign(B[-1, , drop = FALSE])) * (B[-1, , drop = FALSE] != 0)
  }
  else{
    seqK <- seq_len(K)
    fk <- proportions(n)
    alpha <- obj$alpha2
    grdB <- list()
    nrmB <- 0
    for(k in seqK) {
      X <- obj$Zipt[seq_len(n[k]), id_X, k, lambda.id, rho.id]
      R <- obj$R[seq_len(n[k]), , k, lambda.id, rho.id]
      Tht <- coef(obj, type = "Theta", class.id = k, lambda.id = lambda.id, rho.id = rho.id)
      grdB[[k]] <- fk[k] * (t(X) %*% R %*% Tht) / n[k]
      
      B <- coef(obj, type = "B", class.id = k, lambda.id = lambda.id, rho.id = rho.id)
      nrmB <- nrmB + B[-1, , drop = FALSE]^2
    }
    nrmB <- sqrt(nrmB)
    
    eps <- list()
    for(k in 1:K) {
      B <- coef(obj, type = "B", class.id = k, lambda.id = lambda.id, rho.id = rho.id)
      eps[[k]] <- (grdB[[k]] - 
                     (1 - alpha) * lambda * B[-1, , drop = FALSE] / nrmB - 
                     alpha * lambda * sign(B[-1, , drop = FALSE])) * (B[-1, , drop = FALSE] != 0)
    }
  }
  maxEps <- if(is.list(eps)) 
    sapply(eps, \(.eps) max(abs(.eps), na.rm = TRUE)) 
  else max(abs(eps), na.rm = TRUE)
  list("eps" = eps, "maxEps" = maxEps)
}
gradTht <- function(obj, lambda.id, rho.id) {
  K <- length(obj$Z)
  
  rho <- obj$rho[rho.id]
  
  n <- obj$nobs
  p <- obj$nres
  id_X <- obj$InfoStructure$id_X
  id_Y <- obj$InfoStructure$id_Y
  
  if(K == 1) {
    Tht <- coef(obj, type = "Theta", lambda.id = lambda.id, rho.id = rho.id)
    S <- obj$S[id_Y, id_Y, K, lambda.id, rho.id]
    
    grdTht <- solve(Tht) - S 
    eps <- (grdTht - rho * sign(Tht)) * (Tht != 0)
  }
  else{
    seqK <- seq_len(K)
    fk <- proportions(n)
    alpha <- obj$alpha1
    
    grdTht <- list()
    nrmTht <- 0
    for(k in seqK) {
      Tht <- coef(obj, type = "Theta", class.id = k, lambda.id = lambda.id, rho.id = rho.id)
      S <- obj$S[id_Y, id_Y, k, lambda.id, rho.id]
      
      grdTht[[k]] <- fk[k] * (solve(Tht) - S)
      nrmTht <- nrmTht + Tht^2
    }
    nrmTht <- sqrt(nrmTht)
    
    eps <- list()
    for(k in 1:K) {
      Tht <- coef(obj, type = "Theta", class.id = k, lambda.id = lambda.id, rho.id = rho.id)
      eps[[k]] <- (grdTht[[k]] - 
                     (1 - alpha) * rho * Tht / nrmTht - 
                     alpha * rho * sign(Tht)) * (Tht != 0)
    }
  }
  maxEps <- if(is.list(eps)) 
    sapply(eps, \(.eps) max(abs(.eps[outer(1:p, 1:p, "<")]), na.rm = TRUE)) 
  else max(abs(eps[outer(1:p, 1:p, "<")]), na.rm = TRUE)
  list("eps" = eps, "maxEps" = maxEps)
}

##### jcglasso S3 methods #####
print.jcglasso <- function (x, digits = 3L, ...){
  K <- length(x$Z)
  p <- x$nresp
  q <- x$npred
  nrho <- x$nrho
  rho <- x$rho
  nlambda <- x$nlambda
  lambda <- x$lambda
  alpha1 <- x$alpha1
  alpha2 <- x$alpha2
  alpha3 <- x$alpha3
  nu <- x$nu
  
  dots <- list(...)
  if (is.null(dots$print.gap)) dots$print.gap <- 2L
  if (is.null(dots$quote)) dots$quote <- FALSE
  if (is.null(dots$row.names)) dots$row.names <- FALSE
  
  if (nrho == 1L | nlambda == 1L) {
    df.B <- apply(x$dfB, 2, function(x) as.vector(t(x[p + 1L, , ])) + p)
    if(is.vector(df.B)) df.B <- t(df.B)
    df.B <- rowSums(df.B)
    df.Tht <- apply(x$dfTht, 1, function(x) as.vector(t(x)) + p)
    if(is.vector(df.Tht)) df.Tht <- t(df.Tht)
    df.Tht <- rowSums(df.Tht)
    df.Omg <- if(!is.null(x$dfOmg)) {
      df.Omg <- apply(x$dfOmg, 1, function(x) as.vector(t(x)) + q)
      if(is.vector(df.Omg)) df.Omg <- t(df.Omg)
      df.Omg <- rowSums(df.Omg)
    } else 0L
    df <- df.B + df.Tht + df.Omg
    df.max <- p * (p + 1L) / 2 + (q + 1L) * p
    if(!is.null(x$dfOmg)) df.max <- df.max + q * (q + 1L) / 2
    df.max <- df.max * K
    df.per <- formatC(round(df / df.max * 100, digits = digits), format = 'f', digits = digits)
    df.per <- paste("(", df.per, "%)", sep = "")
    tbl <- data.frame(rho = rho, lambda = lambda, nu = if(is.null(nu)) rho else nu, 
                      alpha1 = alpha1, alpha2 = alpha2, alpha3 = alpha3,
                      df.Tht, df.B, df.Omg, 
                      df, df.per)
    names(tbl)[11L] <- "(df%)"
    if (q == 0L) tbl <- tbl[, -c(2L, 3L, 5L, 6L, 8L, 9L), drop = FALSE]
    cat("\nCall:  ", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
    do.call(function(...) print.data.frame(tbl, digits = digits, ...), dots)
    cat("\n---\n")
  } 
  else {
    lambda <- rep(lambda, each = nrho)
    rho <- rep(rho, times = nlambda)
    
    df.B <- apply(x$dfB, 2, function(x) as.vector(t(x[p + 1L, , ])) + p)
    df.B <- rowSums(df.B)
    df.Tht <- apply(x$dfTht, 1, function(x) as.vector(t(x)) + p)
    df.Tht <- rowSums(df.Tht)
    df.Omg <- if(!is.null(x$dfOmg)) {
      df.Omg <- apply(x$dfOmg, 1, function(x) as.vector(t(x)) + q)
      df.Omg <- rowSums(df.Omg)
    } else 0L
    df <- df.B + df.Tht + df.Omg
    df.max <- p * (p + 1) / 2 + (q + 1L) * p
    if(!is.null(x$dfOmg)) df.max <- df.max + q * (q + 1) / 2
    df.max <- df.max * K
    df.per <- formatC(round(df / df.max * 100, digits = digits), format = 'f', digits = digits)
    df.per <- paste("(", df.per, "%)", sep = "")
    tbl <- data.frame(rho = rho, lambda = lambda, nu = if(is.null(nu)) rho else nu, 
                      alpha1 = alpha1, alpha2 = alpha2, alpha3 = alpha3,
                      df.Tht, df.B, df.Omg, 
                      df, df.per)
    names(tbl)[11L] <- "(df%)"
    # names(tbl) <- c("rho1", "rho2", "df", "", "N. Comp.")
    if (q == 0L) tbl <- tbl[, -c(2L, 3L, 5L, 6L, 8L, 9L), drop = FALSE]
    cat("\nCall:  ", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
    #tbl.list <- split(tbl, f = rep(seq_len(nlambda), each = nrho), drop = FALSE)
    if (nlambda <= nrho) f <- rep(seq_len(nlambda), each = nrho)
    else f <- rep(seq_len(npen1), times = npen2)
    tbl.list <- split(tbl, f = f, drop = FALSE)
    do.call(function(...) print.listof(tbl.list, digits = digits, ...), dots)
    cat("---\n")
  }
  
  cat("\nmodel:", sQuote(x$model))
  cat("\n    nObs:", x$nobs)
  cat("\n   nResp:", x$nresp)
  cat("\n   nPred:", x$npred, "\n\n")
  invisible(tbl)
}

coef.jcglasso <- function(object, type = c("all", "B", "Theta", "Omega"), 
                          class.id, rho.id, lambda.id, drop = TRUE, ...) {
  type <- match.arg(type)
  K <- length(object$nobs)
  nrho <- object$nrho
  nlambda <- object$nlambda
  id_X <- object$InfoStructure$id_X
  id_Y <- object$InfoStructure$id_Y
  if (!is.logical(drop)) stop(sQuote("drop"), " is not an object of type ", sQuote("integer"))
  if (missing(class.id)) class.id <- seq_len(K)
  else {
    if(!is.vector(class.id)) stop(sQuote("class.id"), " is not a vector")
    if(any(abs(as.integer(class.id) - class.id) > 0)) stop(sQuote("class.id"), " is not an object of type ", dQuote("integer"))
    if(min(class.id) <= 0) stop("some entry in ", sQuote("class.id"), " is not a positive integer")
    if(max(class.id) > K) stop("some entry in ", sQuote("class.id"), " is larger than ", sQuote(class.id))
  }
  
  if (missing(rho.id)) rho.id <- seq_len(nrho)
  else {
    if(!is.vector(rho.id)) stop(sQuote("rho.id"), " is not a vector")
    if(any(abs(as.integer(rho.id) - rho.id) > 0)) stop(sQuote("rho.id"), " is not an object of type ", dQuote("integer"))
    if(min(rho.id) <= 0) stop("some entry in ", sQuote("rho.id"), " is not a positive integer")
    if(max(rho.id) > nrho) stop("some entry in ", sQuote("rho.id"), " is larger than ", sQuote(nrho))
  }
  if (missing(lambda.id)) lambda.id <- seq_len(nlambda)
  else {
    if(!is.vector(lambda.id)) stop(sQuote("lambda.id"), " is not a vector")
    if(any(abs(as.integer(lambda.id) - lambda.id) > 0)) stop(sQuote("lambda.id"), " is not an object of type ", dQuote("integer"))
    if(min(lambda.id) <= 0) stop("some entry in ", sQuote("lambda.id"), " is not a positive integer")
    if(max(lambda.id) > nlambda) stop("some entry in ", sQuote("lambda.id"), " is larger than ", sQuote(nlambda))
  }
  
  if (type == "all") {
    out.coef <- list(B = object$B[, , class.id, lambda.id, rho.id, drop = FALSE],
                     Theta = object$Tht[id_Y, id_Y, class.id, lambda.id, rho.id, drop = FALSE],
                     Omega = object$Omega[, , class.id, lambda.id, rho.id, drop = FALSE])
    if(drop) out.coef <- lapply(out.coef, drop)
  } else {
    out.coef <- switch(type,
                       B = object$B[, , class.id, lambda.id, rho.id, drop = FALSE],
                       Theta = object$Tht[id_Y, id_Y, class.id, lambda.id, rho.id, drop = FALSE],
                       Omega = object$Omega[, , class.id, lambda.id, rho.id, drop = FALSE])
    if(drop) out.coef <- drop(out.coef)
  }
  out.coef
}

fitted.jcglasso <- function(object, class.id, rho.id, lambda.id, drop = TRUE, ...) {
  K <- length(object$nobs)
  nrho <- object$nrho
  nlambda <- object$nlambda
  id_X <- object$InfoStructure$id_X
  id_Y <- object$InfoStructure$id_Y
  if (!is.logical(drop)) stop(sQuote("drop"), " is not an object of type ", sQuote("integer"))
  if (missing(class.id)) class.id <- seq_len(K)
  else {
    if(!is.vector(class.id)) stop(sQuote("class.id"), " is not a vector")
    if(any(abs(as.integer(class.id) - class.id) > 0)) stop(sQuote("class.id"), " is not an object of type ", dQuote("integer"))
    if(min(class.id) <= 0) stop("some entry in ", sQuote("class.id"), " is not a positive integer")
    if(max(class.id) > K) stop("some entry in ", sQuote("class.id"), " is larger than ", sQuote(class.id))
  }
  
  if (missing(rho.id)) rho.id <- seq_len(nrho)
  else {
    if(!is.vector(rho.id)) stop(sQuote("rho.id"), " is not a vector")
    if(any(abs(as.integer(rho.id) - rho.id) > 0)) stop(sQuote("rho.id"), " is not an object of type ", dQuote("integer"))
    if(min(rho.id) <= 0) stop("some entry in ", sQuote("rho.id"), " is not a positive integer")
    if(max(rho.id) > nrho) stop("some entry in ", sQuote("rho.id"), " is larger than ", sQuote(nrho))
  }
  if (missing(lambda.id)) lambda.id <- seq_len(nlambda)
  else {
    if(!is.vector(lambda.id)) stop(sQuote("lambda.id"), " is not a vector")
    if(any(abs(as.integer(lambda.id) - lambda.id) > 0)) stop(sQuote("lambda.id"), " is not an object of type ", dQuote("integer"))
    if(min(lambda.id) <= 0) stop("some entry in ", sQuote("lambda.id"), " is not a positive integer")
    if(max(lambda.id) > nlambda) stop("some entry in ", sQuote("lambda.id"), " is larger than ", sQuote(nlambda))
  }
  out.fitted <- object$mu[, id_Y, class.id, lambda.id, rho.id]
  if(drop) out.fitted <- drop(out.fitted)
  out.fitted
}

residuals.jcglasso <- function(object, type = c("observed", "working"), class.id, rho.id, lambda.id, drop = TRUE, ...) {
  type <- match.arg(type)
  K <- length(object$nobs)
  nrho <- object$nrho
  nlambda <- object$nlambda
  id_Y <- object$InfoStructure$id_Y
  nobs <- object$nobs
  
  if (!is.logical(drop)) stop(sQuote("drop"), " is not an object of type ", sQuote("integer"))
  if (missing(class.id)) class.id <- seq_len(K)
  else {
    if(!is.vector(class.id)) stop(sQuote("class.id"), " is not a vector")
    if(any(abs(as.integer(class.id) - class.id) > 0)) stop(sQuote("class.id"), " is not an object of type ", dQuote("integer"))
    if(min(class.id) <= 0) stop("some entry in ", sQuote("class.id"), " is not a positive integer")
    if(max(class.id) > K) stop("some entry in ", sQuote("class.id"), " is larger than ", sQuote(class.id))
  }
  
  if (missing(rho.id)) rho.id <- seq_len(nrho)
  else {
    if(!is.vector(rho.id)) stop(sQuote("rho.id"), " is not a vector")
    if(any(abs(as.integer(rho.id) - rho.id) > 0)) stop(sQuote("rho.id"), " is not an object of type ", dQuote("integer"))
    if(min(rho.id) <= 0) stop("some entry in ", sQuote("rho.id"), " is not a positive integer")
    if(max(rho.id) > nrho) stop("some entry in ", sQuote("rho.id"), " is larger than ", sQuote(nrho))
  }
  if (missing(lambda.id)) lambda.id <- seq_len(nlambda)
  else {
    if(!is.vector(lambda.id)) stop(sQuote("lambda.id"), " is not a vector")
    if(any(abs(as.integer(lambda.id) - lambda.id) > 0)) stop(sQuote("lambda.id"), " is not an object of type ", dQuote("integer"))
    if(min(lambda.id) <= 0) stop("some entry in ", sQuote("lambda.id"), " is not a positive integer")
    if(max(lambda.id) > nlambda) stop("some entry in ", sQuote("lambda.id"), " is larger than ", sQuote(nlambda))
  }
  
  out.residuals <- object$R
  if (type == "observed") {
    for(k in class.id){
      R <- event2(object$Z)[[k]]
      if(all(R == 0)) break
      residuals.dim <- dim(out.residuals)
      for (i in 1:residuals.dim[4L]) {
        for (j in 1:residuals.dim[5L]) 
          out.residuals[seq_len(nobs[k]), , k, i, j][R[, id_Y] != 0] <- NA
      }
    }
  }
  
  out.residuals <- out.residuals[, , class.id, lambda.id, rho.id, drop = FALSE]
  if(drop) out.residuals <- drop(out.residuals)
  out.residuals
}

summary.jcglasso <- function(object, GoF = AIC, print.all = TRUE, digits = 3L,  ...){
  if (!is.element(class(GoF), c("function", "GoF2")))
    stop (sQuote("GoF"), " is not a goodness-of-fit function (AIC or BIC) neither an object of class ", sQuote("GoF"))
  dots <- list(...)
  n <- object$nobs
  p <- object$nresp
  q <- object$npred
  K <- length(n)
  if (is.function(GoF)) {
    if (is.null(dots$type)) dots$type <- ifelse(q == 0L, "FD", "CC")
    GoF.name <- deparse(substitute(GoF))
    if (!is.element(GoF.name, c("AIC", "BIC")))
      stop(sQuote(GoF.name), " is not a valid function. Please, use ", sQuote("AIC"), " or ", sQuote("BIC"))
    GoF <- switch(GoF.name,
                  AIC = do.call(function(...) AIC(object, ...), dots),
                  BIC = do.call(function(...) BIC(object, ...), dots))
  }
  if (!is.vector(print.all)) stop(sQuote("print.all"), " is not a vector")
  if (length(print.all) != 1L) stop(sQuote("print.all"), " is not a vector of length ", sQuote("1"))
  if (!is.logical(print.all)) stop(sQuote("print.all"), " is not an object of type ", sQuote("logical"))
  rho <- object$rho
  nrho <- object$nrho
  lambda <- object$lambda
  nlambda <- object$nlambda
  nu <- object$nu
  alpha1 <- object$alpha1
  alpha2 <- object$alpha2
  alpha3 <- object$alpha3
  if (is.null(dots$print.gap)) dots$print.gap <- 2L
  if (is.null(dots$quote)) dots$quote <- FALSE
  if (is.null(dots$row.names)) dots$row.names <- FALSE
  
  if (nrho == 1L | nlambda == 1L) {
    df.B <- apply(object$dfB, 2, function(x) as.vector(t(x[p + 1L, , ])) + p)
    if(is.vector(df.B)) df.B <- t(df.B)
    df.B <- rowSums(df.B)
    df.Tht <- apply(object$dfTht, 1, function(x) as.vector(t(x)) + p)
    if(is.vector(df.Tht)) df.Tht <- t(df.Tht)
    df.Tht <- rowSums(df.Tht)
    df.Omg <- if(!is.null(object$dfOmg)) {
      df.Omg <- apply(object$dfOmg, 1, function(x) as.vector(t(x)) + q)
      if(is.vector(df.Omg)) df.Omg <- t(df.Omg)
      df.Omg <- rowSums(df.Omg)
    } else 0L
    df <- df.B + df.Tht + df.Omg
    df.max <- p * (p + 1L) / 2 + (q + 1L) * p
    if(!is.null(object$dfOmg)) df.max <- df.max + q * (q + 1L) / 2
    df.max <- df.max * K
    df.per <- formatC(round(df / df.max * 100, digits = digits), format = 'f', digits = digits)
    df.per <- paste("(", df.per, "%)", sep = "")
    val <- drop(GoF$value_gof)
    rnk <- rank(val)
    rnk_min <- which.min(rnk)
    rnk <- as.character(rnk)
    rnk[-rnk_min] <- paste(rnk[-rnk_min], "  ")
    rnk[rnk_min] <- paste(rnk[rnk_min], "<-")
    lambda.id <- ifelse (nlambda == 1L, 1L, rnk_min)
    rho.id <- ifelse (nrho == 1L, 1L, rnk_min)
    tbl <- data.frame(rho = rho, lambda = lambda, nu = if(is.null(nu)) rho else nu, 
                      alpha1 = alpha1, alpha2 = alpha2, alpha3 = alpha3, 
                      df.Tht, df.B, df.Omg, 
                      df, df.per, val, Rank = rnk)
    names(tbl)[11L] <- "(df%)"
    names(tbl)[12L] <- GoF$type
    if (q == 0L) tbl <- tbl[, -c(2L, 3L, 5L, 6L, 8L, 9L), drop = FALSE]
    cat("\nCall:  ", paste(deparse(object$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
    if (print.all) {
      do.call(function(...) print.data.frame(tbl, digits = digits, ...), dots)
      cat("\n---\n")
    }
  } 
  else {
    lambda <- rep(lambda, each = nrho)
    lambda.id <- rep(seq_len(nlambda), each = nrho)
    rho <- rep(rho, times = nlambda)
    rho.id <- rep(seq_len(nrho), times = nlambda)
    
    df.B <- apply(object$dfB, 2, function(x) as.vector(t(x[p + 1L, , ])) + p)
    df.B <- rowSums(df.B)
    df.Tht <- apply(object$dfTht, 1, function(x) as.vector(t(x)) + p)
    df.Tht <- rowSums(df.Tht)
    df.Omg <- if(!is.null(object$dfOmg)) {
      df.Omg <- apply(object$dfOmg, 1, function(x) as.vector(t(x)) + q)
      df.Omg <- rowSums(df.Omg)
    } else 0L
    df <- df.B + df.Tht + df.Omg
    df.max <- p * (p + 1) / 2 + (q + 1L) * p
    if(!is.null(object$dfOmg)) df.max <- df.max + q * (q + 1) / 2
    df.max <- df.max * K
    df.per <- formatC(round(df / df.max * 100, digits = digits), format = 'f', digits = digits)
    df.per <- paste("(", df.per, "%)", sep = "")
    
    val <- as.vector(t(GoF$value_gof))
    rnk <- rank(val)
    rnk_min <- which.min(rnk)
    rnk <- as.character(rnk)
    rnk[-rnk_min] <- paste(rnk[-rnk_min], "  ")
    rnk[rnk_min] <- paste(rnk[rnk_min], "<-")
    lambda.id <- lambda.id[rnk_min]
    rho.id <- rho.id[rnk_min]
    tbl <- data.frame(rho = rho, lambda = lambda, nu = if(is.null(nu)) rho else nu, 
                      alpha1 = alpha1, alpha2 = alpha2, alpha3 = alpha3, 
                      df.Tht, df.B, df.Omg, 
                      df, df.per, val, Rank = rnk)
    names(tbl)[11L] <- "(df%)"
    names(tbl)[12L] <- GoF$type
    if (q == 0L) tbl <- tbl[, -c(2L, 3L, 5L, 6L, 8L, 9L), drop = FALSE]
    cat("\nCall:  ", paste(deparse(object$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
    tbl.list <- split(tbl, f = rep(seq_len(nlambda), each = nrho))
    if (print.all) do.call(function(...) print.listof(tbl.list, digits = digits, ...), dots)
  }
  
  lbl <- c("model", "nClasses", 
           "nObs", "nResp", "nPred", 
           "penalty", 
           "rho", "rho.id", "alpha1",
           "lambda", "lambda.id", "alpha2",
           "nu", "alpha3", 
           GoF$type, 
           "df.Tht", "df.B", "df.Omg", "df")  
  
  lbl <- paste("\n", format(lbl, justify = "right"), ":", sep = "")
  cat("\n===============================================================")
  cat("\n\nSummary of the Selected Model\n")
  if (!is.null(object$model)) cat(lbl[1L], sQuote(object$model))
  cat(lbl[2L], K)
  cat(lbl[3L], n)
  cat(lbl[4L], p)
  cat(lbl[5L], q)
  cat(lbl[6L], object$penalty)
  cat(lbl[7], tbl[rnk_min, "rho"])
  cat(lbl[8L], rho.id)
  cat(lbl[9L], alpha1)
  if (q > 0L) {
    cat(lbl[10L], tbl[rnk_min, "lambda"]) 
    cat(lbl[11L], lambda.id)
    cat(lbl[12L], alpha2)
    cat(lbl[13L], tbl[rnk_min, "nu"])
    cat(lbl[14L], alpha3)
  }
  cat(lbl[15L], tbl[rnk_min, GoF$type])
  cat(lbl[16L], tbl[rnk_min, "df.Tht"])
  if (q > 0L) {
    cat(lbl[17L], tbl[rnk_min, "df.B"])
    cat(lbl[18L], tbl[rnk_min, "df.Omg"])
  }
  cat(lbl[19L], tbl[rnk_min, "df"])
  cat("\n\n===============================================================\n\n")
  invisible(list(table = tbl, rho.id = rho.id, lambda.id = lambda.id))
}

plot.jcglasso <- function(x, what = c("Theta", "diag(Theta)", "b0", "B"), penalty = ifelse(nrho >= nlambda, "rho", "lambda"),
                          given = NULL, GoF = AIC, add.labels, matplot.arg1, matplot.arg2, labels.arg, abline.arg, mtext.arg, save.plot,
                          grdev = pdf, grdev.arg, digits = 4L, ...) {
  p <- x$nresp
  q <- x$npred
  K <- length(x$nobs)
  nrho <- x$nrho
  rho <- x$rho
  nlambda <- x$nlambda
  lambda <- x$lambda
  
  # if (nrho1 == 1L & nrho2 == 1L) stop("plot method is not available because ", sQuote("nrho2 = 1"), " and ", sQuote("nrho1 = 1"))
  if (nrho == 1L & nlambda == 1L) { plot.jcglasso.single(x, ...); return(invisible(NULL)) }
  if (inherits(what, "formula")) {
    xy <- function2xy.v2(what)
    penalty <- xy[["penalty"]]
    what <- xy[["what"]]
    given <- xy[["given"]]
  }
  # testing arguments
  penalty <- match.arg(penalty, c("rho", "lambda"))
  if (penalty == "lambda" & nlambda == 1L) stop("plot method is not available because ", sQuote("nlambda = 1"))
  if (penalty == "rho" & nrho == 1L) stop("plot method is not available because ", sQuote("nrho = 1"))
  what <- match.arg(what)
  if (what == "B" & q == 0L) stop("plot method is not available because the number of predictors is zero")
  if (is.null(given)) given <- seq_len(ifelse(penalty == "rho", nlambda, nrho))
  else {
    given.max <- ifelse(penalty == "rho", nlambda, nrho)
    if (!is.vector(given)) stop (sQuote("given"), "is not a vector")
    if (any(abs(as.integer(given) - given) > 0)) stop(sQuote("given"), " is not an object of type ", dQuote("integer"))
    if (min(given) <= 0) stop("some entry in ", sQuote("given"), " is not a positive integer")
    if (any(given > given.max)) stop("some entry in ", sQuote("given"), " is larger than ", sQuote(given.max))
  }
  ngiven <- length(given)
  # testing 'matplot.arg1'
  if (missing(matplot.arg1)) matplot.arg1 <- vector(mode = "list")
  if (is.null(matplot.arg1$type)) matplot.arg1$type <- "l"
  if (is.null(matplot.arg1$col)) matplot.arg1$col <- "black"
  if (is.null(matplot.arg1$lty)) matplot.arg1$lty <- 1L
  if (!is.list(matplot.arg1)) stop(sQuote("matplot.arg1"), " is not an object of type ", dQuote("list"))
  if (is.null(names(matplot.arg1))) stop(sQuote("matplot.arg1"), " is not a named list")
  # testing 'matplot.arg2'
  if (missing(matplot.arg2)) matplot.arg2 <- vector(mode = "list")
  if (is.null(matplot.arg2$type)) matplot.arg2$type <- "l"
  if (is.null(matplot.arg2$col)) matplot.arg2$col <- "gray70"
  if (is.null(matplot.arg2$lty)) matplot.arg2$lty <- 2L
  if (!is.list(matplot.arg2)) stop(sQuote("matplot.arg2"), " is not an object of type ", dQuote("list"))
  if (is.null(names(matplot.arg2))) stop(sQuote("matplot.arg2"), " is not a named list")
  # testing 'labels.arg'
  if (missing(labels.arg)) labels.arg <- vector(mode = "list")
  if (is.null(labels.arg$pos)) labels.arg$pos <- 4L
  if (!is.list(labels.arg)) stop(sQuote("labels.arg"), " is not an object of type ", dQuote("list"))
  if (is.null(names(labels.arg))) stop(sQuote("labels.arg"), " is not a named list")
  # testing 'abline.arg'
  if (missing(abline.arg)) abline.arg <- vector(mode = "list")
  if (is.null(abline.arg$lwd)) abline.arg$lwd <- 2L
  if (is.null(abline.arg$lty)) abline.arg$lty <- 2L
  if (is.null(abline.arg$col)) abline.arg$col <- 2L
  if (!is.list(abline.arg)) stop(sQuote("abline.arg"), " is not an object of type ", dQuote("list"))
  if (is.null(names(abline.arg))) stop(sQuote("abline.arg"), " is not a named list")
  # testing 'mtext.arg'
  if (missing(mtext.arg)) mtext.arg <- vector(mode = "list")
  else {
    if (!is.list(mtext.arg)) stop(sQuote("mtext.arg"), " is not an object of type ", dQuote("list"))
    if (is.null(names(mtext.arg))) stop(sQuote("mtext.arg"), " is not a naned list")
  }
  # testing 'save.plot'
  if (missing(save.plot)) save.plot <- FALSE
  if (length(save.plot) != 1L) stop(sQuote(save.plot), " is not an object of length ", sQuote("1"))
  if (!inherits(save.plot, c("logical", "character"))) stop(sQuote(save.plot), " is not an object of type ", sQuote("logical"), "or ", , sQuote("character"))
  if(inherits(save.plot, "character")) {
    oldPath <- getwd()
    newPath <- save.plot
    setwd(newPath)
    on.exit(setwd(oldPath), add = TRUE, after = TRUE)
    save.plot <- TRUE
  }
  # testing 'grdev'
  if (!is.function(grdev)) stop(sQuote("grdev"), " is not a function")
  else {
    grdev.name <- deparse(substitute(grdev))
    if (!is.element(grdev.name, c("pdf", "postscript", "svg", "bmp", "jpeg", "png", "tiff")))
      stop(sQuote(grdev.name), " is not a valid graphics device for exporting plots. Please, use one of the following functions:\n",
           sQuote("pdf"), ", ", sQuote("postscript"), ", ", sQuote("svg"), ", ", sQuote("bmp"), ", ",
           sQuote("jpeg"), ", ", sQuote("png"), " or ", sQuote("tiff"), ", ")
    if (grdev.name == "postscript") grdev.name <- "ps"
  }
  if (missing(grdev.arg)) grdev.arg <- NULL
  else {
    if (!is.list(grdev.arg)) stop(sQuote("grdev.arg"), " is not an object of type ", sQuote("list"))
  }
  # testing 'digits'
  if (!is.vector(digits)) stop(sQuote("digits"), "is not a vector")
  if (length(digits) > 1L) stop(sQuote("digits"), "is not a vector of length ", sQuote(1))
  if (abs(as.integer(digits) - digits) > 0) stop(sQuote("digits"), " is not an object of type ", dQuote("integer"))
  if (digits < 0L) stop(sQuote("digits"), " is not a positive integer")
  # section: goodness-of-fit
  dots <- list(...)
  if (is.function(GoF)) {
    GoF.name <- deparse(substitute(GoF))
    if (!is.element(GoF.name, c("AIC", "BIC")))
      stop(sQuote(GoF.name), " is not a valid function. Please, use ", sQuote("AIC"), " or ", sQuote("BIC"))
    GoF <- switch(GoF.name,
                  AIC = do.call(function(...) AIC(x, ...), dots),
                  BIC = do.call(function(...) BIC(x, ...), dots))
  }
  # testing 'add.labels'
  if (missing(add.labels)) add.labels <- !is.null(GoF)
  if (!is.logical(add.labels)) stop(sQuote(add.labels), " is not an object of type ", sQuote("logical"))
  if (length(add.labels) != 1L) stop(sQuote(add.labels), " is not an object of length ", sQuote("1"))
  if (add.labels & is.null(GoF)) warning("labels are not added to the plot because GoF is NULL" )
  # starting main code
  tp <- x[[penalty]]
  tp.given <- x[[ifelse(penalty == "rho", "lambda", "rho")]][given]
  tp.given <- formatC(round(tp.given, digits = digits), format = 'f', digits = digits)
  coef.name <- what
  if (coef.name == "diag(Theta)") coef.name <- "Theta"
  if (coef.name == "b0") coef.name <- "B"
  if(save.plot) {
    nplot <- ifelse(what != "B", ngiven, ngiven * p)
    file.names <- paste0(what, "_path_", formatC(seq_len(nplot), width = nchar(nplot), format = "d", flag = "0"), ".", grdev.name)
    cat("Exporting plots\n")
    pb <- txtProgressBar(min = 0L, max = nplot, style = 3L)
    on.exit(close(pb), add = TRUE, after = TRUE)
  } 
  else {
    if (ngiven > 1L | what == "B") {
      op <- par(no.readonly = TRUE)
      on.exit(par(op), add = TRUE)
      devAskNewPage(TRUE)
    }
  }
  if (what == "Theta") U <- .row(c(p, p)) < .col(c(p, p))
  if (add.labels) {
    if (what == "b0") lbls <- colNames2(x$Z[[1]])$Y
    if (what == "B") lbls <- colNames2(x$Z[[1]])$X
    if (what == "diag(Theta)") {
      lbls <- seq_len(p)
      lbls <- paste(lbls, lbls, sep = ",")
      lbls <- mapply(function(x) bquote(hat(theta)[.(x)]), lbls)
      names(lbls) <- NULL
      lbls <- sapply(lbls, as.expression)
    }
    if (what == "Theta") {
      lbls <- paste(unlist(mapply(seq, from = 1, to = 1:(p-1))), unlist(mapply(rep, x = 2:p, each = 1:(p-1))), sep = ",")
      lbls <- mapply(function(x) bquote(hat(theta)[.(x)]), lbls)
      names(lbls) <- NULL
      lbls <- sapply(lbls, as.expression)
    }
  }
  hh <- 0
  for (m in seq_len(ngiven)) {
    if (penalty == "lambda") {
      tp.lab <- bquote(lambda ~ " | {" ~ rho %~~% .(tp.given[m]) ~ "}")
      Y <- coef.jcglasso(x, type = coef.name, rho.id = given[m], drop = TRUE)
      if (!is.null(GoF))  minGoF.id <- which.min(GoF$value_gof[, given[m]])
    }
    else {
      if (q != 0) 
        tp.lab <- bquote(rho ~ " | {" ~ lambda %~~% .(tp.given[m]) ~ "}")
      else tp.lab <- bquote(rho)
      Y <- coef.jcglasso(x, type = coef.name, lambda.id = given[m], drop = TRUE)
      if (!is.null(GoF)) minGoF.id <- which.min(GoF$value_gof[given[m], ])
    }
    ##################
    # plotting 'b0'
    if (what == "b0") {
      hh <- hh + 1L
      if (save.plot) {
        if (is.null(grdev.arg)) grdev(file = file.names[hh])
        else do.call(function(...) grdev(file = file.names[hh], ...), grdev.arg)
      } else dev.hold()
      
      xlab <- if (is.null(matplot.arg1$xlab)) tp.lab else matplot.arg1$xlab
      ylab <- if (is.null(matplot.arg1$ylab))  "" else matplot.arg1$ylab
      main <- if (is.null(matplot.arg1$main)) ifelse(q == 0L, "Expected Values", "Intercepts") else matplot.arg1$main
      arg.name <- setdiff(names(matplot.arg1), c("x", "y", "xlab", "ylab", "main"))
      for(k in seq_len(K)){
        if (q == 0L) yy <- t(Y[ , k, ])
        else yy <- t(Y[1L, , k, ])
        do.call(what = function(...) matplot(x = tp, y = yy, xlab = xlab, ylab = ylab, main = main, ...), args = matplot.arg1[arg.name])
        if (!is.null(GoF)) {
          if(add.labels) {
            labels <- if(is.null(labels.arg$labels)) lbls else labels.arg$labels
            arg.name2 <- setdiff(names(labels.arg), "labels")
            do.call(what = function(...) text(x = tp[minGoF.id], y = yy[minGoF.id, ], labels = labels, ...), args = labels.arg[arg.name2])
          }
          abline.arg$v <- tp[minGoF.id]
          do.call(what = abline, args = abline.arg)
          mtext.arg$text <- GoF$type
          mtext.arg$at <- tp[minGoF.id]
          do.call(what = mtext, args = mtext.arg)
        }
      }
      
      if (save.plot) setTxtProgressBar(pb, hh) else dev.flush()
      if(save.plot) dev.off()
    }
    ##################
    # plotting 'B'
    if (what == "B") {
      varname <- colNames2(x$Z[[1L]])$Y
      for (k in seq_len(p)) {
        hh <- hh + 1L
        if (save.plot) {
          if (is.null(grdev.arg)) grdev(file = file.names[hh])
          else do.call(function(...) grdev(file = file.names[hh], ...), grdev.arg)
        } else dev.hold()
        
        xlab <- if (is.null(matplot.arg1$xlab)) tp.lab else matplot.arg1$xlab
        ylab <- if (is.null(matplot.arg1$ylab)) "Regression Coefficients" else matplot.arg1$ylab
        main <- if (is.null(matplot.arg1$main)) paste0("Response Variable ", varname[k]) else matplot.arg1$main
        arg.name <- setdiff(names(matplot.arg1), c("x", "y", "xlab", "ylab", "main"))
        for(o in seq_len(K)) {
          yy <- Y[-1L, k, o, ]
          if (!is.vector(yy)) yy <- t(yy)
          if (is.null(GoF)) do.call(what = function(...) matplot(x = tp, y = yy, xlab = xlab, ylab = ylab, main = main, ...), args = matplot.arg1[arg.name])
          else {
            matplot(x = tp, y = yy, xlab = xlab, ylab = ylab, main = main, type = "n")
            A <- abs(yy[minGoF.id, ]) > 0
            if(any(!A)) do.call(what = function(...) matpoints(x = tp, y = yy[, !A], ...), args = matplot.arg2[arg.name])
            if(any(A)) {
              do.call(what = function(...) matpoints(x = tp, y = yy[, A], ...), args = matplot.arg1[arg.name])
              if(add.labels) {
                labels <- if (is.null(labels.arg$labels)) lbls[A] else labels.arg$labels
                arg.name2 <- setdiff(names(labels.arg), "labels")
                do.call(what = function(...) text(x = tp[minGoF.id], y = yy[minGoF.id, A], labels = labels, ...), args = labels.arg[arg.name2])
              }
            }
            abline.arg$v <- tp[minGoF.id]
            do.call(what = abline, args = abline.arg)
            mtext.arg$text <- GoF$type
            mtext.arg$at <- tp[minGoF.id]
            do.call(what = mtext, args = mtext.arg)
          }
        }
        
        if (save.plot) setTxtProgressBar(pb, hh) else dev.flush()
        if(save.plot) dev.off()
      }
    }
    ##################
    # plotting 'diag(Theta)'
    if (what == "diag(Theta)") {
      hh <- hh + 1L
      if (save.plot) {
        if (is.null(grdev.arg)) grdev(file = file.names[hh])
        else do.call(function(...) grdev(file = file.names[hh], ...), grdev.arg)
      } else dev.hold()
      
      xlab <- if (is.null(matplot.arg1$xlab)) tp.lab else matplot.arg1$xlab
      ylab <- if (is.null(matplot.arg1$ylab)) "Diagonal Elements" else matplot.arg1$ylab
      main <- if (is.null(matplot.arg1$main)) "Precision Matrix" else matplot.arg1$main
      arg.name <- setdiff(names(matplot.arg1), c("x", "y", "xlab", "ylab", "main"))
      for(k in seq_len(K)){
        yy <- t(apply(Y[, , k, ], 3L, diag))
        do.call(what = function(...) matplot(x = tp, y = yy, xlab = xlab, ylab = ylab, main = main, ...), args = matplot.arg1[arg.name])
        if (!is.null(GoF)) {
          if(add.labels) {
            labels <- if(is.null(labels.arg$labels)) lbls else labels.arg$labels
            arg.name2 <- setdiff(names(labels.arg), "labels")
            do.call(what = function(...) text(x = tp[minGoF.id], y = yy[minGoF.id, ], labels = labels, ...), args = labels.arg[arg.name2])
          }
          abline.arg$v <- tp[minGoF.id]
          do.call(what = abline, args = abline.arg)
          mtext.arg$text <- GoF$type
          mtext.arg$at <- tp[minGoF.id]
          do.call(what = mtext, args = mtext.arg)
        }
      }
      
      if (save.plot) setTxtProgressBar(pb, hh) else dev.flush()
      if(save.plot) dev.off()
    }
    ##################
    # plotting 'Theta'
    if (what == "Theta") {
      hh <- hh + 1L
      if (save.plot) {
        if (is.null(grdev.arg)) grdev(file = file.names[hh])
        else do.call(function(...) grdev(file = file.names[hh], ...), grdev.arg)
      } else dev.hold()
      
      xlab <- if (is.null(matplot.arg1$xlab)) tp.lab else matplot.arg1$xlab
      ylab <- if (is.null(matplot.arg1$ylab)) "Off-Diagonal Elements" else matplot.arg1$ylab
      main <- if (is.null(matplot.arg1$main)) "Precision Matrix" else matplot.arg1$main
      arg.name <- setdiff(names(matplot.arg1), c("x", "y", "xlab", "ylab", "main"))
      for(k in seq_len(K)) {
        yy <- t(apply(Y[, , k, ], 3L, function(M) M[U]))
        if (is.null(GoF)) do.call(what = function(...) matplot(x = tp, y = yy, xlab = xlab, ylab = ylab, main = main, ...), args = matplot.arg1[arg.name])
        else {
          matplot(x = tp, y = yy, xlab = xlab, ylab = ylab, main = main, type = "n")
          A <- abs(yy[minGoF.id, ]) > 0
          if(any(!A)) do.call(what = function(...) matpoints(x = tp, y = yy[, !A], ...), args = matplot.arg2[arg.name])
          if(any(A)) {
            do.call(what = function(...) matpoints(x = tp, y = yy[, A], ...), args = matplot.arg1[arg.name])
            if(add.labels) {
              labels <- if(is.null(labels.arg$labels)) lbls[A] else labels.arg$labels
              arg.name2 <- setdiff(names(labels.arg), "labels")
              do.call(what = function(...) text(x = tp[minGoF.id], y = yy[minGoF.id, A], labels = labels, ...), args = labels.arg[arg.name2])
            }
          }
          abline.arg$v <- tp[minGoF.id]
          do.call(what = abline, args = abline.arg)
          mtext.arg$text <- GoF$type
          mtext.arg$at <- tp[minGoF.id]
          do.call(what = mtext, args = mtext.arg)
        }
      }
      
      if (save.plot) setTxtProgressBar(pb, hh) else dev.flush()
      if(save.plot) dev.off()
    }
  }
  invisible(NULL)
}

##### auxiliary function for plot.jcglasso S3 method #####
function2xy.v2 <- function (frml) {
  x <- match.arg(all.vars(frml)[2L], c("rho", "lambda"))
  y <- as.character(terms(frml)[[2L]])
  if (length(y) == 2L) {
    if (y[1L] != "diag" | y[2L] != "Theta") 
      stop(sQuote(paste0(y[1L], "(", y[2L], ")")), " is not available as response. Please, use ", 
           sQuote("diag(Theta)"))
    y <- "diag(Theta)"
  }
  given <- as.character(terms(frml)[[3L]])
  if (length(given) == 1L) 
    given <- NULL
  else {
    if (given[3L] == "rho") {
      if (x == "rho") 
        stop("You can not condition on ", sQuote(given[3L]), 
             ". Please, use ", sQuote("lambda"))
      return(list(x = x, y = y, given = NULL))
    }
    if (given[3L] == "lambda") {
      if (x == "lambda") 
        stop("You can not condition on ", sQuote(given[3L]), 
             ". Please, use ", sQuote("rho"))
      return(list(x = x, y = y, given = NULL))
    }
    given.terms <- gsub("[[:blank:]]", "", given[3L])
    given.terms <- unlist(strsplit(given.terms, split = "\n", 
                                   fixed = TRUE))
    if (length(given.terms) == 1L) 
      given <- eval(parse(text = given.terms))
    else {
      given.terms <- given.terms[2L]
      given.name <- unlist(strsplit(given.terms, "=", fixed = TRUE))[1L]
      if (!is.element(given.name, c("lambda.id", "rho.id"))) 
        stop("Please, use ", sQuote(ifelse(x == "rho", "lambda.id", "rho.id")), " instead of ", sQuote(given.name))
      if ((x == "rho" & given.name == "rho.id") | (x == "lambda" & given.name == "lambda.id")) 
        stop("You can not condition on ", sQuote(given.name), 
             ". Please, use ", sQuote(ifelse(given.name == "rho.id", "lambda.id", "rho.id")))
      given <- eval(parse(text = given.terms))
    }
  }
  out <- list(penalty = x, what = y, given = given)
  out
}

plot.jcglasso.single <- function(x, type = c("Gyy", "Gxy", "Gxx"), ...) {
  K <- length(x$Z)
  type <- match.arg(type)
  id_X <- x$InfoStructure$id_X
  id_Y <- x$InfoStructure$id_Y
  if (is.null(x$npred > 0L) & is.element(type, c("Gxy", "Gxx")))
    stop(sQuote(type), " is not available. Please, use type = ", dQuote("Gyy"))
  mfrow <- if(K == 1) c(1, 2) else c(2, K)
  par(mfrow = mfrow, pty = "s", omi = c(0.3, 0.3, 0.3, 0.3), mai = c(0.3, 0.3, 0.3, 0.3))
  if(type == "Gyy") {
    apply(x$InfoStructure$Adj, 3L, function(x) {
      diag(x[id_Y, id_Y, 1L, 1L]) <- 1
      image(x[id_Y, id_Y, 1L, 1L], col = gray.colors(256), main = "Adjacency Matrix")
    })
    g <- apply(x$InfoStructure$Adj, 3L, function(x) graph.adjacency(x[id_Y, id_Y, 1L, 1L], 
                                                                    mode = "undirected", diag = FALSE))
    layout.grid <- lapply(g, layout.fruchterman.reingold)
    sapply(seq_len(length(g)), function(i) plot(g[[i]], layout = layout.grid[[i]], edge.color = "gray50", 
                                                vertex.color = "red", vertex.size = 3, vertex.label = NA, 
                                                main = "Graph Pattern"))
  } 
  if(type == "Gxx") {
    apply(x$InfoStructure$Adj, 3L, function(x) {
      diag(x[id_X, id_X, 1L, 1L]) <- 1
      image(x[id_X, id_X, 1L, 1L], col = gray.colors(256), main = "Adjacency Matrix")
    })
    g <- apply(x$InfoStructure$Adj, 3L, function(x) graph.adjacency(x[id_X, id_X, 1L, 1L], 
                                                                    mode = "undirected", diag = FALSE))
    layout.grid <- lapply(g, layout.fruchterman.reingold)
    sapply(seq_len(length(g)), function(i) plot(g[[i]], layout = layout.grid[[i]], edge.color = "gray50", 
                                                vertex.color = "red", vertex.size = 3, vertex.label = NA, 
                                                main = "Graph Pattern"))
  }
  if(type == "Gxy") {
    apply(x$InfoStructure$Adj, 3L, 
          function(x) image(x[id_X, id_Y, 1L, 1L], col = gray.colors(256), main = "Adjacency Matrix"))
    g <- apply(x$InfoStructure$Adj, 3, function(x) {
      out <- x[id_X, id_Y, 1L, 1L]
      nmsX <- rownames(out)
      nmsY <- colnames(out)
      eB <- data.frame(from = rep(nmsX, ncol(out)), to = rep(nmsY, each = nrow(out)))[c(out) != 0, ]
      vB <- data.frame(id = c(nmsY, nmsX))
      graph_from_data_frame(d = eB, directed = TRUE, vertices = vB)
    })
    layout.grid <- lapply(g, layout.fruchterman.reingold)
    sapply(seq_len(length(g)), function(i) plot(g[[i]], layout = layout.grid[[i]], edge.color = "gray50", 
                                                vertex.color = "red", vertex.size = 3, vertex.label = NA, 
                                                main = "Graph Pattern"))
  }
  
  invisible(NULL)
}

