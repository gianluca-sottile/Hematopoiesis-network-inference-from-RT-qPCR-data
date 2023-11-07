# Sparse inference of the human hematopoietic system from heterogeneous and partially observed genomic
The source R code for jcglasso. 

## Required R packages:

CRAN packages:
- cglasso
- JGL
- VennDiagram
- MASS
- huge
- lattice
- ggplot2

## 01 - RCode
This folder contains six R files with the main functions needed to implement the jcglasso estimator:
1. datajcggm.R, which implements the new class object
2. jcglasso.R containing the main functions for the proposed estimator
3. jcggm.R containing functions to perform a post-hoc maximum likelihood fit of a selected jcglasso model
4. gof.R containing the functions to compute the goodness-of-fit measure for model selection
5. to_graph.R containing routines to create an igraph object and plot the estimated network
6. test_jcglasso_function.R contains a simple example of how the proposed model works.

## 02 - Simulation
This folder contains 4 subfolders, each of which refers to a block of the simulation study, i.e., 2.1 (a, b, c) and 2.2. In each subfolder, there are script files to reproduce the simulation study with relative plots, one folder in which are located the .RData files storing the results, and where compatible a folder with the figure.

## 03 - Analysis
This folder contains the main script to reproduce the whole analysis divided in:
1. to preprocess the raw data (contained in the subfolder data) as proposed by Psaila et al. (2016);
2. to descript data and to check censoring assumption;
3. to estimate the optimal jcglasso model after seclecting the three main tuning parameters nu, lambda and rho (keeping fixed the alpha' parameters to 0.75);
4. to analyse networks for only the experimentally validated edges.

In the subfolder auxiliary there are all the files relevant to divide transcripts in nuclear activities (responses) and membrane recpetor activities (covariates), to convert transcript IDs to IPA names and to highlight experimentally validated relationships between the selected transctipts. Finally, in the subfolder figs are stored all the figures computed during the whole anlysis.
