
#' Create Empty Mapping and Metadata Tables
#'
#' This function generates empty mapping and metadata tables for a given phyloseq object and project.
#'
#' @param psIN The phyloseq object for which the mapping and metadata tables are to be created.
#' @param project A string representing the name of the project.
#'
#' @return
#' Generates CSV files with empty mapping and metadata tables, saved in a specified directory.
#'
#' @details
#' The function creates two CSV files:
#'   - An empty mapping table with columns like `SampleID`, `StudyID`, `TreatmentGroup`, `Timepoint`, etc.
#'   - An empty metadata table with columns like `StudyID`, `Variation1`, `Variation2`, etc., and rows for different analyses (e.g., `Go_overview`, `Go_ancombc`).
#' The purpose is to provide templates for users to fill in their specific study details and analysis parameters.
#'
#' @examples
#' # Example usage:
#' Go_emptyMap(psIN = my_phyloseq_object, project = "MyMicrobiomeProject")
#'
#' @export

Go_emptyMap <- function(psIN, project) {
  # Validate input parameters
  if (is.null(psIN) || is.null(project)) {
    stop("Invalid input: 'psIN' and 'project' cannot be NULL.")
  }

  # Function to create and save a data frame
  createAndSaveDataFrame <- function(columnNames, fileName) {
    rowCount <- length(SampleID)
    dataFrame <- data.frame(matrix(ncol = length(columnNames), nrow = rowCount))
    colnames(dataFrame) <- columnNames
    
    dataFrame$SampleID <- SampleID
    dataFrame[is.na(dataFrame)] <- ""
    
    filePath <- file.path(map, sprintf("empty.%s.%s.%s.csv", format(Sys.Date(), "%y%m%d"), project, fileName))
    tryCatch({
      write.csv(dataFrame, quote = FALSE, col.names = NA, row.names = FALSE, file = filePath)
      cat(sprintf("%s is saved in %s.\n", fileName, map))
    }, error = function(e) {
      cat("Error in writing file:", e$message, "\n")
    })
  }

  # Create output directory
  map <- file.path("3_map")
  if (!dir.exists(map)) dir.create(map)

  # Initialize common variables
  SampleID <- sample_names(psIN)

  # Create and save the emptyMap
  emptyMapColumns <- c("SampleID", "StudyID", "TreatmentGroup", "Timepoint", "Description", "etc")
  createAndSaveDataFrame(emptyMapColumns, "mapping")

  # Create and save the SCRubMap
  SCRubMapColumns <- c("SampleID", "is_control", "sample_type", "sample_well")
  createAndSaveDataFrame(SCRubMapColumns, "mapping.SCRub")
}

