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

### 0.hy.test.R
This file contains the main functions that are necessary for the implementation of hy.test.

### 1.loading_and_prepocessing.R
This file contains the first script that is necessary to load all the required packages and to preprocess raw data.

### 2.deg_analysis.R
This file contains the second script to run all the differential expression analysis, i.e., Moderated t-test, Significance Analysis of Microarray (sam), Empirical Bayes Analysis of Microarrays (ebam) and hy.test.

### 3.go_enrichment.R
This file contains the third script to run gene onthology enrihment by using all the genes selected in the previous analysis. Genes are selected after a Benjamini & Hochberg correction.

### 4.pubmed_research.R
This file contains the last script for implmenting the hypergeometric test used to the pubmed research to highlight significant terms of GO. Terms are selected after a Benjamini & Hochberg correction.

### TCGA_BRCA75
This folder contains the raw gene expression data for breast cancer. The first .txt contains normal samples and the scond one tumor samples.

### TCGA_KIRC75
This folder contains the raw gene expression data for kidney cancer. The first .txt contains normal samples and the scond one tumor samples.

### Results and outputs
This folder contains two more folders 'brca' and 'kirc'. In each of them, final .RData objects, .csv files and figures of the overall analysis of the paper is included.


