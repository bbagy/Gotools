% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Go_adiv_V2.R
\name{Go_adiv}
\alias{Go_adiv}
\title{Alpha Diversity Analysis for Microbiome Data}
\usage{
Go_adiv(psIN, project, alpha_metrics, name = NULL)
}
\arguments{
\item{psIN}{A phyloseq object containing microbiome data.}

\item{project}{A string representing the project name, used for file naming and directory creation.}

\item{alpha_metrics}{A vector of alpha diversity metrics to be calculated.}
}
\value{
Returns a data frame containing the calculated alpha diversity indices for each sample in the dataset. The result includes both the diversity indices and the corresponding sample metadata.
}
\description{
This function computes various alpha diversity metrics for microbiome data.
}
\details{
The function calculates several alpha diversity indices, including Observed, Chao1, ACE, Shannon, Simpson, InvSimpson, Fisher, and PD (Phylogenetic Diversity: Faith’s PD). It allows for a comprehensive analysis of species richness and evenness within individual samples in a microbiome dataset.
}
\examples{
# Example usage:
adiv_results <- Go_adiv(psIN = my_phyloseq_object,
                        project = "MyMicrobiomeStudy",
                        alpha_metrics = c("Observed", "Shannon", "Simpson"))

}
