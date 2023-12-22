
#' Merge Multiple Data Tables Based on a Common Pattern
#'
#' @param pattern The pattern used to identify and match files for merging.
#' @param file_path The directory path where the files to be merged are located.
#'
#' @details
#' This function aggregates multiple data tables into a single data frame. It searches for files in a specified directory that match a given pattern and merges them row-wise. This function is particularly useful for consolidating similar data spread across multiple files, such as results from repetitive experiments or outputs from batch processes.
#'
#' @return
#' A single data frame consisting of the row-wise merged contents of all files matching the specified pattern in the given directory.
#'
#' @examples
#' merged_data <- Go_mergeTab(pattern = ".csv",
#'                            file_path = "path/to/data/files")
#'
#' @export

Go_mergeTab <- function(pattern, file_path){
  
   # add input files
  path <- file_path
 
  
  filenames <- list.files(path, pattern=pattern);filenames
  sample.names <- sapply(strsplit(filenames, pattern), `[`, 1) ;sample.names
  filenames <- list.files(path, pattern=pattern);filenames
  
  
  cat(sprintf("Files location: %s\n",path))
  cat("=======================================================================\n")
  cat("Merged files:\n")
  cat(sprintf("%s\n",filenames))

  
  # add input files
  df<-{}
  for (sn in sample.names) {
    file <- file.path(path, paste0(sn, pattern))
    df1 <- read.csv(file, row.names=NULL ,check.names=FALSE)
    df <- rbind(df, df1)
  }
  return(df)
}
