# SpaXSeq - a user-friendly R package for the analysis of spatial transcriptomics from GeoMx and Visium Data

## Installation

SpaXSeq can be installed from Github using the following R commands:
```
# install the package
library("devtools")
install_github("SpaXSeq/SpaXSeq")

# launch the application
library(SpaXSeq)
SpaXSeq:::startSpaXSeq()
```
### Prerequisites

Due to the R version, some packages may encounter an error during installation. If so, please manually install individual packages.

## About

SpaXSeq is a step-by-step interactive web application capable of tuning, exploring, and visualizing spatial transcriptomics data generated from two commercial platforms, Nanostring GeoMx DSP-NGS and 10X Genomics Visium.

![SpaXSeq workflow](/figures/workflow.png)


