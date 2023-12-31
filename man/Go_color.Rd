% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Go_color_V3.R
\name{Go_color}
\alias{Go_color}
\title{Generate Color Palette for Taxonomic Data}
\usage{
Go_color(cdf, taxaName)
}
\arguments{
\item{cdf}{Data frame containing taxonomic information (phylum level) and corresponding colors.}

\item{taxaName}{Vector of taxonomic names for which colors are assigned.}
}
\value{
A list containing two elements:
\itemize{
\item \verb{$color_table}: A data frame mapping each taxa to its corresponding HSV color code.
\item \verb{$coloring}: A named vector where names are taxa and values are color codes.
}
}
\description{
Generate Color Palette for Taxonomic Data
}
\details{
This function creates a color palette for taxonomic data, primarily at the phylum level. It assigns distinct colors to
each phylum and ensures visual differentiation between various taxonomic groups.
}
\examples{
color_data <- Go_color(cdf = taxonomic_dataframe, taxaName = vector_of_taxa_names)
# Access the color table
color_table <- color_data$color_table
# Access the coloring vector
coloring_vector <- color_data$coloring

}
