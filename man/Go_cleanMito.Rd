% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Go_cleanMito_V1.R
\name{Go_cleanMito}
\alias{Go_cleanMito}
\title{Remove Mitochondrial Reads from Phyloseq Object}
\usage{
Go_cleanMito(psIN, project)
}
\arguments{
\item{psIN}{Phyloseq object containing microbiome data.}

\item{project}{Name of the project, used for output file naming.}
}
\value{
The function generates a phyloseq object with mitochondrial reads removed, a corresponding FASTA file of sequences, and a CSV file containing the processed ASV table. These files are saved in specified directories.
}
\description{
Remove Mitochondrial Reads from Phyloseq Object
}
\details{
This function processes a phyloseq object to remove mitochondrial reads. It specifically targets reads assigned to the order "Rickettsiales", often considered as mitochondrial in microbiome datasets. The function also removes taxa with NA in the Phylum level and filters out taxa with very low counts.
}
\examples{
Go_cleanMito(psIN = my_phyloseq_object,
             project = "MyMicrobiomeProject")

}
