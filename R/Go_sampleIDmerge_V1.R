#' Go_sampleIDmerge
#'
#' Merges a Phyloseq object with additional sample metadata and optionally aggregates samples based on a specified column.
#'
#' @param psIN A Phyloseq object containing species count data, sample metadata, and optionally taxonomic data.
#' @param map Either a dataframe or a character string representing the path to a CSV file containing additional sample metadata to merge with the Phyloseq object.
#' @param mergeBy A character string specifying the column in the sample metadata by which to aggregate samples.
#' @param project A character string specifying the project name, used for file naming.
#' @param name An optional character string for additional naming in the output file.
#'
#' @return A Phyloseq object that has been merged and optionally aggregated by sample identifiers specified in the `mergeBy` parameter. The function also saves an RDS file of the merged Phyloseq object.
#'
#' @details
#' The function first checks if the `map` parameter is a character string indicating a file path; if so, it reads the file into a dataframe. It then merges this metadata with the original Phyloseq object's sample data. After merging, the function can aggregate samples based on the specified `mergeBy` column. The resulting Phyloseq object, stripped of its original sample data, is saved as an RDS file in a specified directory.
#'
#' @importFrom phyloseq merge_phyloseq sample_data otu_table tax_table
#' @importFrom utils read.csv
#'
#' @examples
#' \dontrun{
#' # Assuming psIN is a Phyloseq object and map is either a dataframe or a path to a metadata CSV file
#' psIN <- globalPatterns
#' map <- "path/to/sample_metadata.csv"
#' result <- Go_sampleIDmerge(psIN, map, "StudyID", "ExampleProject")
#' print(result)
#' }
#' @export

Go_sampleIDmerge <- function(psIN, map, mergeBy, project, name=NULL){

  # Prepare output directory
  out <- file.path("2_rds")
  if (!dir.exists(out)) dir.create(out)

  # Load map data
  if (is.character(map)){
    mapinput <- read.csv(map, row.names=1, check.names=FALSE)
  } else {
    mapinput <- map
  }

  # Merge Phyloseq object with sample metadata
  psIN1 <- merge_phyloseq(psIN, sample_data(mapinput))

  # Merge samples based on the specified column
  ps1.mergedbysample <- merge_samples(psIN1, group = sample_data(psIN1)[[mergeBy]])

  # Remove sample data from the Phyloseq object
  ps1.mergedbysample_no_sample <- phyloseq(otu_table(ps1.mergedbysample),
                                           tax_table(ps1.mergedbysample))  # sample_data excluded

  # Save the merged Phyloseq object as an RDS file
  saveRDS(ps1.mergedbysample_no_sample, sprintf("%s/ps.mergedbystudyID.%s.%s%s.rds",
                                                out, project,
                                                ifelse(is.null(name), "", paste(name, ".", sep = "")),
                                                format(Sys.Date(), "%y%m%d")))

  return(ps1.mergedbysample_no_sample)
}
