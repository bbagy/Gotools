% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Go_correlation_V3.R
\name{Go_correlation}
\alias{Go_correlation}
\title{Perform Correlation Analysis on Phyloseq Data}
\usage{
Go_correlation(
  psIN,
  project,
  rank,
  pvalue = 0.05,
  maingroup = NULL,
  top = NULL,
  Table,
  facet = "TAB",
  X = "GROUP",
  Y = "PS",
  method = "spearman",
  padj = TRUE,
  xanlgle = 90,
  ncols = 23,
  name = NULL,
  orders = NULL,
  height,
  width
)
}
\arguments{
\item{psIN}{Phyloseq object containing the data for analysis.}

\item{project}{Name of the project or analysis.}

\item{rank}{Taxonomic rank to consider for correlation analysis.}

\item{pvalue}{P-value threshold for significance in the analysis.}

\item{maingroup}{Main group for correlation analysis.}

\item{top}{Number of top taxa to consider for analysis (if NULL, all taxa are considered).}

\item{Table}{Data frame containing additional data for correlation with Phyloseq data.}

\item{facet}{Faceting variable for the plot.}

\item{X}{Name of the variable to be used on the X-axis of the plot.}

\item{Y}{Name of the variable to be used on the Y-axis of the plot.}

\item{method}{Correlation method to use, such as "spearman", "kendall", or "pearson".}

\item{padj}{Logical value indicating whether to adjust p-values for multiple testing.}

\item{xanlgle}{Angle for x-axis text in the plot.}

\item{ncols}{Number of columns for facet wrapping in the plot.}

\item{name}{Optional name for the analysis.}

\item{orders}{Order of levels in the categorical variables.}

\item{height}{Height of the output plot.}

\item{width}{Width of the output plot.}
}
\value{
Generates a plot of the correlation analysis and saves it as a PDF file.
}
\description{
This function conducts a correlation analysis on data from a Phyloseq object.
It supports various correlation methods and allows adjustment for multiple testing.
The function can filter the top taxa, apply log transformation, and generate correlation plots.
}
\examples{
# Assuming 'ps' is a Phyloseq object and 'additional_data' is a data frame
Go_correlation(ps, "MyProject", "Genus", pvalue = 0.05, "SampleGroup", top = 10,
               Table = additional_data, facet = "SampleType", X = "EnvFactor",
               Y = "Taxa", method = "spearman", padj = TRUE, xanlgle = 90,
               ncols = 5, height = 10, width = 10)

}
