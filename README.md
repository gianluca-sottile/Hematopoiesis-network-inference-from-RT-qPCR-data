# Sparse-inference-of-the-human-hematopoietic-system-from-heterogeneous-and-partially-observed-genomic
The source R code for jcglasso. 

## Required R packages:

CRAN packages:
- cglasso
- JGL
- VennDiagram
- MASS
- huge
- lattice

## 01 - RCode
This folder contains a single R file with the main functions necessary for implementing the jcglasso estimator.

## 02 - Simulation
This folder contains three subfolders, each of which refers to a block of the simulation study. In each subfolder, there are script files to reproduce the simulation study with relative plots, one folder in which are located the .RData files storing the results, and where compatible a folder with the figure.

## 03 - Analysis
This folder contains the main script to reproduce the whole analysis divided in:
- to preprocess the raw data (contained in the subfolder data) as proposed by Psaila et al. (2016);
- to descript data and to check censoring assumption;
- to estimate the optimal jcglasso model after seclecting the three main tuning parameters nu, lambda and rho (keeping fixed the alpha' parameters to 0.75);
- to analyse networks for only the experimentally validated edges.

In the subfolder auxiliary there are all the files relevant to divide transcripts in nuclear activities (responses) and membrane recpetor activities (covariates), to convert transcript IDs to IPA names and to highlight experimentally validated relationships between the selected transctipts. Finally, in the subfolder figs are stored all the figures computed during the whole anlysis.
