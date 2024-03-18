# Sparse inference of the human hematopoietic system from heterogeneous and partially observed genomic

## Required R packages:

CRAN packages:
- cglasso (you need to install by source the version 2.0.5 in the subfolder **02 - Simulations/SimulationStudy2/CGLASSO**)
- VennDiagram
- MASS
- huge
- caret
- ggplot2
- glmnet
- grpreg
- glasso
- BDgraph
- microbenchmark

## 01 - RCode
This folder contains six R files with the main functions needed to implement the jcglasso estimator:
1. datajcggm.R, which implements the new class object
2. jcglasso.R containing the main functions for the proposed estimator
3. jcggm.R containing functions to perform a post-hoc maximum likelihood fit of a selected jcglasso model
4. gof.R containing the functions to compute the goodness-of-fit measure for model selection
5. to_graph.R containing routines to create an igraph object and plot the estimated network
6. rjcggm.R containing routines to create simulated data
7. test_jcglasso_function.R contains a simple example of how the proposed model works.

## 02 - Simulation
This folder contains two subfolders, each corresponding to a block of the simulation study. Each subdirectory contains script files to reproduce the simulation study and supplementary files with relative plots and, where compatible, a folder containing the figure.

## 03 - Analysis
This folder contains the main script to reproduce the whole analysis divided into:
1. to preprocess the raw data (contained in the subfolder data) as proposed by Psaila et al. (2016);
2. to describe data and to check censoring assumptions;
3. to estimate the optimal jcglasso model after selecting, firstly, the three main tuning parameters nu, lambda and rho and then the mixing parameters alpha;
4. to analyse networks for only the experimentally validated edges.

The auxiliary subfolder contains all the files relevant for separating transcripts into nuclear activities (responses) and membrane receptor activities (covariates), converting transcript IDs into IPA names and highlighting experimentally validated relationships between the selected transcripts. Finally, the figs subfolder contains all the figures calculated throughout the analysis.

## Example
```{r message = FALSE, warning = FALSE}
library(cglasso)
library(MASS)
library(huge)

source("01 - RCode/datajcggm.R")
source("01 - RCode/jcglasso.R")
source("01 - RCode/jcggm.R")
source("01 - RCode/gof.R")
source("01 - RCode/to_graph.R")
source("01 - RCode/rjcggm.R")

n <- 200L                # sample size (integer)
p <- 25L                 # number of response variables (integer)
q <- 10L                 # number of predictors (integer)
K <- 2L                  # number of subpopulation (integer)
ncompPrecision <- 5L     # number of connected components in the precision matrices Omega e Theta for each k (integer)
ncompB <- 2L             # number of connected components in the coefficient regression matrix B for each k (integer)
percConnectedCompK <- .5 # percentages of connected components for Omega, Theta and B across k (numeric vector of length k - 1)
perc.NA <- 0.25          # percentage of missing-at-random data for each column
perc.X.na <- 0.25        # percentage of X columns with missing-at-random data
perc.Y.na <- 0.25        # percentage of Y columns with missing-at-random data
up <- 40.0               # upper censored value for the response matrix
perc.cens <- 0.25        # percentage of censored data for each column
perc.Y.cens <- 0.25      # percentage of Y columns with censored data
Thtmin <- 0.3            # minimum value of the precision matrix Theta
Thtmax <- 0.5            # maximum value of the precision matrix Theta
Omgmin <- 0.3            # minimum value of the precision matrix Omega
Omgmax <- 0.5            # maximum value of the precision matrix Omega
Bmin <- 0.5              # minimum value of coefficient regression matrix B
Bmax <- 0.7              # maximum value of coefficient regression matrix B

set.seed(1234)           # random seed for reproducibility

sim <- rjcggm(n = n, p = p, q = q, K = K, ncompPrecision = ncompPrecision, ncompB = ncompB, 
              percConnectedCompK = percConnectedCompK, perc.NA = perc.NA, perc.X.na = perc.X.na, perc.Y.na = perc.Y.na, 
              up = up, perc.cens = perc.cens, perc.Y.cens = perc.Y.cens, Thtmin = Thtmin, Thtmax = Thtmax, 
              Omgmin = Omgmin, Omgmax = Omgmax, Bmin = Bmin, Bmax = Bmax)

##### datajcggm object with only the X matrix
sim$ZX

##### datajcggm object with both the X and Y matrices
sim$Z

##### produces summaries of an object of class ‘datajcggm’ for the first 
summary(sim$Z)

##### fits the joint graphical lasso model to X datasets with censored and/or missing values
model1 <- jcglasso(data = sim$ZX, nrho = 25L, rho.min.ratio = .01, alpha1 = .5, trace = 1L)
model1

##### choice of the optimal nu parameter fixed alpha = 0.5
par(mfrow = c(2, 2))
# by using penalised estimates
plot(AIC(model1))
plot(BIC(model1))

# by using maximum likelihood estimates
qfun <- QFun2(model1, mle = TRUE)
plot(AIC(model1, Qfun = qfun))
plot(BIC(model1, Qfun = qfun))

##### best model obtained minimizing the BIC criterion
summary(model1, GoF = BIC(model1, Qfun = qfun), print.all = FALSE)
model1best <- select.jcglasso(model1, GoF = BIC(model1, Qfun = qfun))
model1best

##### simple grayscale plot of the precision matrices Omega
plot(model1best)

##### fits the joint conditional graphical lasso model to the datasets with censored and/or missing values
##### fixing nu to its optimal value and alpha3 = 0.5
model2 <- jcglasso(data = sim$Z, nrho = 5L, rho.min.ratio = .01, alpha1 = .5,
                   nlambda = 5L, lambda.min.ratio = .01, alpha2 = .5,
                   nu = model1best$rho, alpha3 = model1best$alpha1, trace = 1L)

##### check for KKT conditions
gradB(model2, lambda.id = 5, rho.id = 5)$maxEps
gradTht(model2, lambda.id = 5, rho.id = 5)$maxEps

##### choice of the optimal lambda and rho parameters fixed alpha1 = alpha2 = 0.5
##### and fixing nu to its optimal value and alpha3 = 0.5
par(mfrow = c(2, 2))
# by using penalised estimates
plot(AIC(model2))
plot(BIC(model2))

# by using maximum likelihood estimates
qfun2 <- QFun2(model2, mle = TRUE)
plot(AIC(model2, Qfun = qfun2))
plot(BIC(model2, Qfun = qfun2))

##### best model obtained minimizing the BIC criterion
summary(model2, GoF = BIC(model2, Qfun = qfun2), print.all = FALSE)
model2best <- select.jcglasso(model2, GoF = BIC(model2, Qfun = qfun2))
model2best

##### check for KKT conditions
gradB(model2best, lambda.id = 1, rho.id = 1)$maxEps
gradTht(model2best, lambda.id = 1, rho.id = 1)$maxEps

##### simple grayscale plot of the precision matrices Theta
plot(model2best)

#### perform post-hoc maximum likelihood refitting of a selected joint 
#### conditional graphical lasso model with censored and/or missing values.
model2bestMLE <- jcggm(model2best)
model2bestMLE

#### extract the main components of the model fit
coef(model2bestMLE, type = "B")
coef(model2bestMLE, type = "Theta")
coef(model2bestMLE, type = "Omega")

#### returns a named list of graphs using the results of an R object of class ‘jcglasso’ or ‘jcggm’.
graphs <- to_graph2(model2bestMLE)
graphs

par(mfrow = c(1, 2))
# plot of an undirected graph representing the conditional dependence structure among the p 
# response variables (i.e., Theta), highlighting the top 3 connected components
plot(graphs, type = "Gyy", highlight.connections = 3L)
# plot of a directed graph representing the effetcs of the q predictors onto the p 
# response variables (i.e., B), highlighting the top 3 connected components
plot(graphs, type = "Gxy", highlight.connections = 3L)
# plot of an undirected graph representing the conditional dependence structure among the q covariates (i.e., Omega)
plot(graphs, type = "Gxx")
# plot of both the undirected graph representing the conditional dependence structure among the p response variables (i.e.,Theta) 
# and the directed graph representing the effects of the q predictors on the p response variables (i.e., B).
plot(graphs, type = "conditional")
# overall plot of both the undirected graph representing the conditional dependence structure among the p response variables (i.e., Theta) 
# and the undirected graph representing the conditional dependence structure among the q covariates (i.e., Omega), 
# plus the directed graph representing the effects of the q predictors on the p response variables (i.e., B).
plot(graphs, type = "bipartite")
```

