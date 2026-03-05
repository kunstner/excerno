FIXES:
replaced defunc tidyverse commands

# Introduction

This package is used as a classifier for determining the origin of a mutation, specifically
for samples that have been preserved using formalin-fixation paraffin-embedding (FFPE).
This preservation process introduces artificial FFPE-caused mutations. The new variants
(mostly C to T type mutations) transform the origin mutation profile, creating difficulty
in detecting the original mutation signature. Excerno quantifies the abundance of the 
FFPE mutational signature and uses Bayes’ formula to filter FFPE artifacts.

Excerno also includes the ability to simulate cancer samples by randomly generating
mutations that corresponds to the FFPE mutational signature and a mutational signature
from the Catalogue Of Somatic Mutations In Cancer.

__Additional data:__

* Simulated cancer samples
* The FFPE signature

# Installation Guide

## Dependencies

Make sure you have all of these packages installed already:

* BiocManager
* BSgenome
* ggplot2
* MutationalPatterns
* stringr
* bedr
* NMF
* BSgenome.Hsapiens.UCSC.hg38
* vcfR
* R.utils
* tidyverse

Notice that for installing some of the packages you need the BiocManager installed, for example

```
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
```

## Installing the package

```
install.packages("devtools")
devtools::install_github("jdavilal/excerno", dependencies = TRUE, build_vignettes = TRUE)
```

# Resources

* Instructions on how to use the library can be found in the [vignette](vignettes/excerno-intro.md). Also trough ```browseVignettes("excerno")```
* The analysis, code, and figures for our [paper](https://github.com/jdavilal/curi_2021/blob/main/paper_analysis.md)
* Our interactive shiny app is [here](https://mitche7.shinyapps.io/excerno/)

# Overview of package workflow

![Workflow image](https://github.com/jdavilal/excerno/blob/master/inst/img/method.png)
