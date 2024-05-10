
#' Convert CSV Data to Phyloseq Object
#'
#' This function reads a CSV file containing OTU (Operational Taxonomic Unit) data and taxonomy information, then converts it into a phyloseq object.
#'
#' @param csv The file path of the CSV file containing OTU data and taxonomy information.
#' @param project A string representing the name of the project.
#'
#' @return
#' A phyloseq object created from the provided CSV data.
#'
#' @details
#' The function reads the CSV file, splits the data into numeric (OTU) and non-numeric (taxonomy) matrices, and creates a phyloseq object. The phyloseq object is then saved as an RDS file in a directory named '2_rds'.
#'
#' @examples
#' # Example usage:
#' ps_object <- Go_tabTops(csv = "path/to/otu_data.csv", project = "MyMicrobiomeProject")
#'
#' @export

Go_tabTops <- function(csv, project){
  rds <- file.path("2_rds")
  if(!file_test("-d", rds)) dir.create(rds)


  tab <- read.csv(csv, row.names = NULL, check.names = FALSE)

  tab.cleaned <- tab
  for(cleaned in c("__no_feature", "__ambiguous","__too_low_aQual","__not_aligned", "__alignment_not_unique")){
    tab.cleaned <- subset(tab.cleaned, locus_tag != cleaned)
  }

  # Check if a column is numeric
  is_numeric <- sapply(tab.cleaned, is.numeric)

  # Separate the data into two matrices based on the column type
  otu <- as.matrix(tab.cleaned[, is_numeric])
  tax <- as.matrix(tab.cleaned[, !is_numeric])

  ps <- phyloseq(otu_table(otu, taxa_are_rows=T),  tax_table(tax));ps

  saveRDS(ps, sprintf("%s/ps.tabTops.%s.%s.rds", rds,
                      project,format(Sys.Date(), "%y%m%d")))

  return(ps)

}
