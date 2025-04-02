#' Go_SCRubMap
#'
#' Maps SCRub data to Illumina SampleSheet data, augments it with additional metadata, and outputs plate-level files.
#'
#' @param SCRubfile The filename of the SCRub file located within the specified map directory.
#' @param IlluSamplesheet The filename of the Illumina SampleSheet, assumed to be in the same directory as the SCRub file.
#' @param IlluLocation The full path to another Illumina SampleSheet containing additional metadata not present in the primary SampleSheet.
#' @param Sample_Type Default sample type to be assigned in the augmented data.
#' @param map_dir The directory containing the SCRub file and Illumina SampleSheet. Default is "3_map".
#'
#' @return Writes augmented SCRub data to CSV files organized by plate (as determined by index sets) within the specified map directory. Does not return any data structures into the R environment.
#'
#' @details
#' The function first checks for the existence of the specified files within the given directory. It then reads and processes the SCRub data, mapping it against the primary Illumina SampleSheet to align sample identifiers and merge relevant metadata. Additional metadata is pulled from a second, more detailed SampleSheet if provided. The function also infers control status and sample type from sample identifiers, and finally, it segments the data by index set and saves it into separate CSV files.
#'
#' This function is particularly useful in sequencing data processing pipelines where sample identification and metadata integration are required before further analysis.
#'
#' @importFrom readr read_csv
#' @importFrom dplyr filter select mutate left_join
#' @importFrom utils read.csv write.csv
#'
#' @examples
#' \dontrun{
#' Go_SCRubMap(
#'   SCRubfile = "sample_scrub.csv",
#'   IlluSamplesheet = "sample_sheet.csv",
#'   IlluLocation = "full_sample_sheet.csv",
#'   Sample_Type = "experimental",
#'   map_dir = "path/to/map_dir"
#' )
#' }
#' @export

Go_SCRubMap <- function(SCRubfile, IlluSamplesheet, IlluLocation, Sample_Type, map_dir = "3_map") {
  
  # 파일 경로 수동 지정 가능하게
  scrub_file <- file.path(map_dir, SCRubfile)
  samplesheet_file <- file.path(map_dir, IlluSamplesheet)
  
  if (file.exists(scrub_file) && file.exists(samplesheet_file)) {
    
    data1 <- read.csv(scrub_file, row.names = 1)
    
    if (!"ID" %in% colnames(data1)) {
      data1$ID <- sub("_S\\d+$", "", rownames(data1))
      data1$ID <- gsub("-", "_", data1$ID)
    }
    
    # SampleSheet 불러오기
    lines <- readLines(samplesheet_file)
    data_start <- which(grepl("^\\[Data\\]", lines)) + 1
    
    if (length(lines) >= data_start + 1) {
      samplesheet <- read.csv(
        text = paste(lines[data_start:length(lines)], collapse = "\n"),
        stringsAsFactors = FALSE
      )
      samplesheet$Sample_ID <- gsub("-", "_", samplesheet$Sample_ID)
      
      available_cols <- intersect(colnames(samplesheet), c("Sample_ID", "I7_Index_ID", "I5_Index_ID"))
      
      if (!"Sample_ID" %in% available_cols) {
        stop("Sample_ID not found in Illumina SampleSheet")
      }
      
      sample_match <- match(data1$ID, samplesheet$Sample_ID)
      new_data <- samplesheet[sample_match, available_cols[available_cols != "Sample_ID"], drop = FALSE]
      merged <- cbind(data1, new_data)
      
      # IlluLocation (full samplesheet) 처리
      lines2 <- readLines(IlluLocation)
      data_start2 <- which(grepl("^\\[Data\\]", lines2)) + 1
      samplesheet_full <- read.csv(
        text = paste(lines2[data_start2:length(lines2)], collapse = "\n"),
        stringsAsFactors = FALSE
      )
      
      fill_meta <- function(i7, i5, sheet, field) {
        idx <- which(sheet$I7_Index_ID == i7 & sheet$I5_Index_ID == i5)
        if (length(idx) == 1) {
          return(sheet[[field]][idx])
        } else {
          return(NA)
        }
      }
      
      merged$sample_well <- mapply(fill_meta,
                                   merged$I7_Index_ID,
                                   merged$I5_Index_ID,
                                   MoreArgs = list(sheet = samplesheet_full, field = "Sample_Well"))
      
      merged$Index_sets <- mapply(fill_meta,
                                  merged$I7_Index_ID,
                                  merged$I5_Index_ID,
                                  MoreArgs = list(sheet = samplesheet_full, field = "Index_sets"))
      
      # is_control / sample_type 지정
      id_lower <- tolower(merged$ID)
      merged$is_control <- grepl("neg|pos|negative|positive", id_lower)
      merged$sample_type <- Sample_Type
      merged$sample_type[merged$is_control & grepl("neg|negative", id_lower)] <- "control negative"
      merged$sample_type[merged$is_control & grepl("pos|positive", id_lower)] <- "control positive"
      
      # plate별 저장
      index_set_levels <- unique(na.omit(merged$Index_sets))
      for (iset in index_set_levels) {
        sub_df <- merged[merged$Index_sets == iset, , drop = FALSE]
        plate_file <- file.path(map_dir, sub("\\.csv$", paste0(".filled_", iset, ".csv"), basename(scrub_file)))
        write.csv(sub_df, file = plate_file, row.names = TRUE)
        message(sprintf("Saved plate-level SCRub: %s", basename(plate_file)))
      }
      
    } else {
      stop("Malformed SampleSheet: [Data] section not found")
    }
    
  } else {
    stop("SCRub or SampleSheet file not found.")
  }
}




