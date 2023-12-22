#' Create Phyloseq Object from CSV File
#'
#' This function reads a CSV file, separates the data into numeric and non-numeric components,
#' and constructs a `phyloseq` object. The numeric data is assumed to be OTU (Operational Taxonomic Unit)
#' counts, while the non-numeric data is treated as taxonomic information.
#'
#' @param csv The path to the CSV file containing the microbiome data.
#'   The CSV file should have a specific format: columns with numeric data representing OTUs and
#'   columns with non-numeric data representing taxonomic information.
#'
#' @return A `phyloseq` object constructed from the OTU counts and taxonomic data in the CSV file.
#'
#' @examples
#' # Assuming you have a CSV file 'my_data.csv' in the correct format:
#' ps_object <- Go_tabTops_markdown("path/to/my_data.csv")
#'
#' @export
#' @importFrom phyloseq otu_table tax_table phyloseq
#' @importFrom utils read.csv

Go_tabTops_markdown <- function(csv){

  tab <- read.csv(csv, row.names = NULL, check.names = FALSE)
  
  # Check if a column is numeric
  is_numeric <- sapply(tab, is.numeric)
  
  # Separate the data into two matrices based on the column type
  otu <- as.matrix(tab[, is_numeric])
  tax <- as.matrix(tab[, !is_numeric])
  
  ps <- phyloseq(otu_table(otu, taxa_are_rows=T),  tax_table(tax));ps
  
  #saveRDS(ps, sprintf("%s/ps.tabTops.%s.%s.rds", rds, project,format(Sys.Date(), "%y%m%d")))
  
  return(ps)
  
}
