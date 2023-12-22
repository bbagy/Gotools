
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

Go_emptyMap <- function(psIN, project){
  
  # out dir
  map <- file.path("3_map") 
  if(!file_test("-d", map)) dir.create(map)
  
  
  # empty map table
  SampleID <- sample_names(psIN)
  # Contamination <- sample_names(psIN)
  StudyID <- sample_names(psIN)
  
  # emptyMap <- data.frame(SampleID, StudyID)
  
  column.names <- c("SampleID",	"StudyID", "TreatmentGroup", "Timepoint","Description","etc")
  col.count <- length(column.names)
  row.count <- length(StudyID)
  
  emptyMap <- data.frame(matrix(ncol = col.count, nrow = row.count))
  colnames(emptyMap) <- c("SampleID",	"StudyID", "TreatmentGroup", "Timepoint","Description","etc")
  
  emptyMap$SampleID <- SampleID
  emptyMap$StudyID <- StudyID
  
  
  emptyMap[is.na(emptyMap)] <- ""
  
  
  cat(sprintf("empty map is saved in %s.\n",map))
  cat("                                                       \n")
  write.csv(emptyMap, quote = FALSE, col.names = NA, row.names = F,
            file=sprintf("%s/empty.%s.%s.mapping.csv",map,format(Sys.Date(), "%y%m%d"), project,sep="/"))
  
  
  # empty metadata table
  column.names <- c("StudyID", "Variation1", "Variation2","etc")
  col.count <- length(column.names)
  
  # 	"Go_overview","Go_ancombc","Go_deseq2","Go_box","Go_bdiv",	"Go_barchart","Go_linear","Go_clme","Go_perm",
  analysis <- c("type",	"baseline",	"Go_reg", "Go_mirkat", "Go_lmem","Confounder")

  row.count <- length(analysis)
  
  emptyMetadata <- data.frame(matrix(ncol = col.count, nrow = row.count))
  colnames(emptyMetadata) <- column.names
  rownames(emptyMetadata) <- analysis


  for(an in analysis){
    if (an == "type"){
      emptyMetadata[c(an), ] <- c("", "factor", "numeric", "factor")
    }else if(an == "baseline"){
      emptyMetadata[c(an), ] <- c("", "control", "before", "male")
    }else{
      emptyMetadata[c(an), ] <- c("no", "no", "yes", "yes")
    }
  }
  
  #cat(sprintf("empty metadata is saved in %s.\n",map))
  #cat("                                                       \n")
  #write.csv(emptyMetadata, quote = FALSE, col.names = NA,  row.names = T,
  #          file=sprintf("%s/emptyControlpanel.%s.%s.csv",map, project,format(Sys.Date(), "%y%m%d"),sep="/"))
} 



