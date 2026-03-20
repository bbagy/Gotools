% Generated manually to match R/Go_lollipopPlot.R
\name{Go_lollipopPlot}
\alias{Go_lollipopPlot}
\title{Generate Ranked Lollipop Plots for DA and ConDA Results}
\usage{
Go_lollipopPlot(
  project,
  file_path = NULL,
  files = NULL,
  result = NULL,
  p_cutoff = 0.05,
  mycols = NULL,
  name = NULL,
  font = 10,
  height = NULL,
  width = 8
)
}
\arguments{
\item{project}{Project name or identifier.}

\item{file_path}{Optional directory path to scan for CSV files.}

\item{files}{Optional filename pattern used together with \code{file_path}.}

\item{result}{Optional result directory returned by a Go DA-family tool or a
ConDA bridge directory.}

\item{p_cutoff}{Numeric cutoff used to keep features. Adjusted p-values
(\code{qval}) are used when available; otherwise raw \code{pval} is used.}

\item{mycols}{Optional two-color vector for negative / positive direction.}

\item{name}{Optional plot label used in output filenames.}

\item{font}{Base font size.}

\item{height}{Deprecated. Height is now determined automatically from the
number of displayed features.}

\item{width}{PDF width.}
}
\value{
The function saves PDF files and returns \code{NULL} invisibly.
}
\description{
This function follows the same directory-first contract as
\code{\link{Go_volcanoPlot}}. It reads CSV result files or bridge CSV files,
recognizes the originating DA tool, and generates ranked lollipop plots into
\verb{pdf/DA_plot} or \verb{pdf/ConDa_plot}.
}
\details{
When \code{result} is supplied, the function scans that directory for CSV
outputs. When \code{file_path} and \code{files} are supplied, it scans the
matching files. Otherwise it searches the dated project table directories for
\code{deseq2}, \code{aldex2}, \code{ancom2}, \code{corncob}, and
\code{maaslin2}.

The function auto-detects the DA tool from the filename, harmonizes effect /
p-value columns, filters features by the supplied p-value cutoff, builds signed
ranking scores, and saves PDF lollipop plots
under \verb{<project_YYMMDD>/pdf/DA_plot/} or
\verb{<project_YYMMDD>/pdf/ConDa_plot/}.
}
\examples{
\dontrun{
Go_lollipopPlot(project = "MyProject")

Go_lollipopPlot(
  project = "MyProject",
  result = "MyProject_250101/table/deseq2",
  p_cutoff = 0.05,
  name = "subsetA"
)
}
}
