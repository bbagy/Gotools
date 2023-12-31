% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Go_perm_V16.R
\name{Go_perm}
\alias{Go_perm}
\title{Global Permutation Analysis for Microbiome Data}
\usage{
Go_perm(
  psIN,
  cate.vars,
  project,
  distance_metrics,
  mul.vars = FALSE,
  name = NULL
)
}
\arguments{
\item{psIN}{A phyloseq object containing microbiome data.}

\item{cate.vars}{A vector of categorical variables for which PERMANOVA will be performed.}

\item{project}{A string representing the project name, used for file naming and directory creation.}

\item{distance_metrics}{A vector of distance metrics to be used in the analysis.}

\item{mul.vars}{A boolean value indicating whether to use multiple variables in the PERMANOVA model.}

\item{name}{An optional string for naming the output files.}
}
\value{
Returns a data frame containing the results of the PERMANOVA analysis, including degrees of freedom, sums of squares, R-squared values, F-model statistics, p-values, and adjusted p-values.
}
\description{
This function performs global permutation (PERMANOVA) analysis on microbiome data.
}
\details{
The function conducts PERMANOVA analysis to understand the impact of various categorical variables on microbial community composition. It supports multiple distance metrics and can handle both single and multiple variables in the model. The function outputs the results in a CSV file and returns them as a data frame.
}
\examples{
# Example usage:
permanova_results <- Go_perm(psIN = my_phyloseq_object,
                             cate.vars = c("TreatmentGroup", "AgeGroup"),
                             project = "MyMicrobiomeStudy",
                             distance_metrics = c("bray", "unifrac"),
                             mul.vars = FALSE,
                             name = "MyAnalysis")

}
