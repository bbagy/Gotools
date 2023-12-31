% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Go_DA_plot_v4.R
\name{Go_DA_plot}
\alias{Go_DA_plot}
\title{Generate Differential Abundance (DA) Plots}
\usage{
Go_DA_plot(
  project,
  file_path,
  files,
  plot = "volcano",
  fdr,
  fc,
  mycols = NULL,
  name,
  overlaps = 10,
  font,
  height,
  width
)
}
\arguments{
\item{project}{Name of the project or analysis.}

\item{file_path}{Directory path containing input files for analysis.}

\item{files}{Pattern to match filenames for analysis.}

\item{plot}{Type of plot to generate: "volcano", "maplot", or "forest".}

\item{fdr}{False Discovery Rate threshold for significance.}

\item{fc}{Fold change threshold for including in the plot.}

\item{mycols}{Optional color palette for the plot.}

\item{name}{Optional name for the analysis.}

\item{overlaps}{Maximum number of overlaps allowed in text repelling.}

\item{font}{Font size for plot text.}

\item{height}{Height of the output plot.}

\item{width}{Width of the output plot.}
}
\value{
Generates DA plots based on specified parameters and saves them as PDF files.
}
\description{
This function generates various plots (volcano, MA, forest) for visualizing differential abundance (DA) data.
}
\details{
The function processes input files containing DA analysis results and generates specified types of plots. It supports custom color schemes and handles multiple comparisons.
}
\examples{
# Example usage:
Go_DA_plot(project = "MyProject",
           file_path = "path/to/data",
           files = "result_pattern.*\\\\.csv",
           plot = "volcano",
           fdr = 0.05,
           fc = 1)

}
