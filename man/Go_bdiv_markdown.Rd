% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Go_bdiv_markdown.R
\name{Go_bdiv_markdown}
\alias{Go_bdiv_markdown}
\title{Beta Diversity Analysis}
\usage{
Go_bdiv_markdown(
  psIN,
  cate.vars,
  project,
  orders,
  distance_metrics,
  cate.conf = NULL,
  plot = "PCoA",
  ellipse = "yes",
  statistics = "yes",
  mycols = NULL,
  paired = NULL,
  combination = NULL,
  shapes = NULL,
  ID = NULL,
  facet = NULL,
  name = NULL,
  addnumber = TRUE
)
}
\arguments{
\item{psIN}{Phyloseq object containing the OTU/ASV counts and associated sample data.}

\item{cate.vars}{Vector of column names in the sample data of 'psIN' representing categorical variables.}

\item{project}{Name or identifier for the project or dataset being analyzed.}

\item{orders}{Optional ordering for the levels of the categorical variables.}

\item{distance_metrics}{Vector of distance metrics to be used for beta diversity analysis.}

\item{cate.conf}{Optional vector of column names in the sample data of 'psIN' representing categorical confounding variables.}

\item{plot}{Type of plot to generate, default is "PCoA".}

\item{ellipse}{Option to add ellipses to the plot ("yes" or "no").}

\item{statistics}{Option to perform statistical analysis ("yes" or "no").}

\item{mycols}{Color palette for the plot.}

\item{paired}{Optional parameter for paired analysis.}

\item{combination}{Numeric value for group combination analysis.}

\item{shapes}{Optional vector for shape aesthetics in the plot.}

\item{ID}{Optional parameter for identifying points on the plot.}

\item{facet}{Optional parameter for facetting the plot.}

\item{name}{Optional name for the analysis output.}

\item{addnumber}{Option to add sample size numbers to the plot categories.}
}
\value{
A list of ggplot objects, each representing a beta diversity plot for the specified distance metric and categorical variable.
}
\description{
This function performs beta diversity analysis on a given phyloseq object.
It supports various diversity metrics, group combinations, and visualization options such as PCoA plots with ellipses.
}
\examples{
# Assuming 'ps' is a phyloseq object with appropriate data
results <- Go_bdiv_markdown(ps, cate.vars = c("Group"),
                            project = "my_project",
                            orders = NULL,
                            distance_metrics = c("bray", "unifrac"),
                            plot = "PCoA",
                            ellipse = "yes")

}
