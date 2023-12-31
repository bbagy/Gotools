% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Go_getSigASVs.R
\name{Go_getSigASVs}
\alias{Go_getSigASVs}
\title{Extract Significant ASVs (Amplicon Sequence Variants) and Generate FASTA File}
\usage{
Go_getSigASVs(
  project,
  file_path,
  files = ".csv",
  SigASVs,
  sampledata,
  target.column,
  target.group,
  name = NULL
)
}
\arguments{
\item{project}{String representing the project name.}

\item{file_path}{Path to the directory containing the files.}

\item{files}{File extension filter (default is ".csv").}

\item{SigASVs}{File path of the significant ASVs data.}

\item{sampledata}{Data frame containing sample information.}

\item{target.column}{The target column to filter (e.g., "deseq2", "ancom").}

\item{target.group}{Groups to consider (e.g., "up", "down").}

\item{name}{Optional name for the output file.}
}
\value{
A data frame of significant ASVs with corresponding FASTA files in the specified output directory.
}
\description{
Extract Significant ASVs (Amplicon Sequence Variants) and Generate FASTA File
}
\details{
This function filters ASVs based on significance from DESeq2 or ANCOM analysis. It creates FASTA files for each
significant ASV group (e.g., upregulated or downregulated). This is useful for further sequence-based analysis or
phylogenetic studies.
}
\examples{
Go_getSigASVs(project = "ProjectName",
             file_path = "/path/to/files",
             SigASVs = "/path/to/SigASVs.csv",
             sampledata = sample_dataframe,
             target.column = "deseq2",
             target.group = c("up", "down"))

}
