# Sparse inference of the human hematopoietic system from heterogeneous and partially observed genomic
The source R code for jcglasso. 

## Required R packages:

CRAN packages:
- cglasso
- JGL
- VennDiagram
- MASS
- huge
- caret
- ggplot2
- glmnet
- grpreg
- glasso

## 01 - RCode
This folder contains six R files with the main functions needed to implement the jcglasso estimator:
1. datajcggm.R, which implements the new class object
2. jcglasso.R containing the main functions for the proposed estimator
3. jcggm.R containing functions to perform a post-hoc maximum likelihood fit of a selected jcglasso model
4. gof.R containing the functions to compute the goodness-of-fit measure for model selection
5. to_graph.R containing routines to create an igraph object and plot the estimated network
6. test_jcglasso_function.R contains a simple example of how the proposed model works.

## 02 - Simulation
This folder contains two subfolders, each of which refers to a block of the simulation study. In each subfolder, there are script files to reproduce the simulation study with relative plots and, where compatible, a folder with the figure.

## 03 - Analysis
This folder contains the main script to reproduce the whole analysis divided into:
1. to preprocess the raw data (contained in the subfolder data) as proposed by Psaila et al. (2016);
2. to descript data and to check censoring assumptions;
3. to estimate the optimal jcglasso model after selecting, firstly, the three main tuning parameters nu, lambda and rho and then the mixing parameters alpha;
4. to analyse networks for only the experimentally validated edges.

In the subfolder auxiliary, there are all the files relevant to dividing transcripts into nuclear activities (responses) and membrane receptor activities (covariates), converting transcript IDs to IPA names and highlighting experimentally validated relationships between the selected transcripts. Finally, in the subfolder figs are stored all the figures computed during the whole analysis.
