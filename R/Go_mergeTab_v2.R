
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

Go_mergeTab <- function(pattern, file_path, addTimePoint = TRUE) {
  # pattern 예: ".Yael_Rejection.csv"
  escape_regex <- function(x) gsub("([][{}()+*^$|\\.^-])", "\\\\\\1", x)

  path <- file_path
  filenames <- list.files(path, pattern = pattern)

  cat(sprintf("Files location: %s\n", path))
  cat("=======================================================================\n")
  cat("Merged files:\n")
  if (length(filenames) == 0) {
    cat("(no files matched)\n")
    return(data.frame())
  }
  cat(sprintf("%s\n", filenames))

  # timepoint은 ".[timepoint][pattern]"에서 [timepoint] 토큰을 의미
  # 즉, 파일명 끝부분 "... .WK1_Post .Yael_Rejection.csv" 의 가운데 토큰
  pat_lit <- escape_regex(pattern)
  time_re <- paste0("^.*\\.([^./]+)", pat_lit, "$")  # 마지막 pattern 바로 앞 토큰 캡처

  out_list <- vector("list", length(filenames))
  warn_idx <- integer(0)

  for (i in seq_along(filenames)) {
    fn <- filenames[i]
    fp <- file.path(path, fn)
    df1 <- read.csv(fp, row.names = NULL, check.names = FALSE, stringsAsFactors = FALSE)

    if (addTimePoint) {
      if (grepl(time_re, fn)) {
        time_val <- sub(time_re, "\\1", fn)
      } else {
        time_val <- NA_character_
        warn_idx <- c(warn_idx, i)
      }
      df1$Timepoint  <- time_val
      df1$SourceFile <- fn
    }
    out_list[[i]] <- df1
  }

  if (length(warn_idx)) {
    warning(sprintf(
      "Timepoint 추출 실패: %s",
      paste(filenames[warn_idx], collapse = ", ")
    ))
  }

  do.call(rbind, out_list)
}
