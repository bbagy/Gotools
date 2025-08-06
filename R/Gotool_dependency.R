
#' Install and Load Necessary Packages for Gotools
#'
#' This function automates the process of installing and loading a predefined set of CRAN and Bioconductor packages necessary for the Gotools package. It checks if each package is already installed; if not, it installs the package and then loads it into the R session.
#'
#' @details
#' The function handles two types of packages: CRAN packages and Bioconductor packages. It first installs and loads CRAN packages, followed by Bioconductor packages. For Bioconductor package installation, `BiocManager` is used.
#'
#' If `crayon` is installed, success messages are printed in color. The function is designed to be used at the beginning of a session to ensure all dependencies for Gotools are available.
#'
#' @return
#' This function does not return anything. It invisibly loads the required packages into the R session.
#'
#' @examples
#' Gotool_dependency()
#'
#' @export

Gotool_dependency <- function() {
  # Packages to install from CRAN
  cran_packages <- c("magick", "ape", "Boruta", "car", "cluster", "CLME", "compositions",
                     "cowplot", "crayon", "caret", "colorspace", "digest",
                     "data.table", "devtools", "doParallel", "ellipse", "emmeans",
                     "e1071", "gplots", "ggplot2", "grid", "gridExtra", "ggrepel",
                     "doRNG", "ggalluvial", "ggforce", "Hmisc", "irlba", "huge",  #"igraph",
                     "irr", "lme4", "lmerTest", "nnet", "MLmetrics",
                     "Matrix", "magrittr", "MASS", "missForest", "nlme",
                     "phangorn", "pheatmap", "pkgconfig", "dplyr", "parallel", "pscl",
                     "plotly", "pdftools",  "rlang", "randomForest", #"rfUtilities",
                     "readxl", "RColorBrewer", "ROCR", "reshape", "reshape2", "yarrr",
                     "stringi", "tidyverse", "vegan", "VGAM", "picante", "zoo",
                     "RcppZiggurat", "Rfast", "survival", "withr", "knitr", "kableExtra", "DT")

  # Packages to install from Bioconductor
  bioconductor_packages <- c("phyloseq", "microbiome", "Rhtslib", "dada2", "dplyr",
                             "ggpubr", "ggfortify", "genefilter", "ggpmisc", "S4Vectors",
                             "ShortRead", "illuminaio", "rstatix", "useful", "DECIPHER",
                             "ComplexHeatmap","DESeq2", "ALDEx2","scater","ANCOMBC") #"ComplexHeatmap",

  # Function to install and load CRAN packages
  install_load_cran <- function(package) {
    if (!require(package, character.only = TRUE, quietly = TRUE)) {
      install.packages(package)
    }
    library(package, character.only = TRUE, quietly = TRUE)
  }

  # Function to install and load Bioconductor packages
  install_load_bioc <- function(package) {
    if (!require(package, character.only = TRUE, quietly = TRUE)) {
      if (!require("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager")
      }
      BiocManager::install(package, force = TRUE, ask = FALSE, type = "source")
    }
    library(package, character.only = TRUE, quietly = TRUE)
  }

  # Installing and loading CRAN packages
  for (package in cran_packages) {
    install_load_cran(package)
  }

  # Installing and loading Bioconductor packages
  for (package in bioconductor_packages) {
    install_load_bioc(package)
  }

  # Print success message
  if (require("crayon", quietly = TRUE)) {
    cat(crayon::blue("#--------------------------------------------------------------# \n"))
    cat(crayon::blue("#------       General analysis Of microbiome (Go)        ------# \n"))
    cat(crayon::blue("#------    Quick statistics and visualization tools      ------# \n"))
    cat(crayon::blue("#------                (with R markdown)                 ------# \n"))
    cat(crayon::blue("#--------------------------------------------------------------# \n"))
    cat(crayon::yellow("All the required packages were installed and loaded.\n"))
    cat(crayon::blue("#--------------------------------------------------------------# \n"))
  } else {
    cat("#--------------------------------------------------------------# \n")
    cat("All the required packages were installed and loaded.\n")
    cat("#--------------------------------------------------------------# \n")
  }
}
