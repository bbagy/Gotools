#' Go_tabByTimepoint
#'
#' This function reshapes a data frame to summarize the presence or values at different time points for each subject. It can either mark presence with binary indicators (0 or 1) or use a specific variable's values to fill the time points.
#'
#' @param df A data frame containing the data to be analyzed.
#' @param SubjectID A character string specifying the column name in `df` that contains the subject IDs.
#' @param Timepoint A character string specifying the column name in `df` that contains the time points.
#' @param orders A character vector specifying the order of the time points. This determines the factor levels for the time points.
#' @param filled_by A character string specifying the column name in `df` that should be used to fill the time point values. If NULL, presence will be indicated by 0 and 1. Default is NULL.
#'
#' @return A data frame in wide format, summarizing the data by time points for each subject.
#'
#' @examples
#' \dontrun{
#' df <- your_data_frame
#' orders <- c("Baseline", "Month1", "Month3")
#' df_wide <- Go_tabByTimepoint(df, SubjectID = "SubjectID", Timepoint = "Timepoint", orders = orders)
#' df_wide_with_values <- Go_tabByTimepoint(df, SubjectID = "SubjectID", Timepoint = "Timepoint", orders = orders, filled_by = "Measurement")
#' }
#'
#' @export


Go_tabByTimepoint <- function(df, project, SubjectID, Timepoint, orders, filled_by=NULL){
  # out dir 설정
  out <- file.path("3_map")
  if(!file_test("-d", out)) dir.create(out)
  
  # Timepoint 순서 지정
  df[[Timepoint]] <- factor(df[[Timepoint]], levels = orders)
  
  if (is.null(filled_by)) {
    # filled_by가 NULL일 때 0과 1로 표시하고 Sum 열 추가
    df_wide <- df %>%
      group_by(across(all_of(c(SubjectID, Timepoint)))) %>%
      summarize(value = 1, .groups = "drop") %>%
      spread(key = Timepoint, value = value, fill = 0) %>%
      mutate(Sum = rowSums(across(all_of(orders))))  # 모든 timepoint가 모였는지 확인하기 위한 Sum 열 추가
  } else {
    # filled_by가 지정된 경우 해당 열의 값으로 채움
    df_wide <- df %>%
      group_by(across(all_of(c(SubjectID, Timepoint)))) %>%
      summarize(value = dplyr::first(.data[[filled_by]]), .groups = "drop") %>%
      spread(key = Timepoint, value = value, fill = 0)
  }
  
  # 결과 CSV로 저장
  write.csv(df_wide, quote = FALSE, row.names = FALSE,
            file = sprintf("%s/%s.%s%s.table_by_timepoint.csv", out,
                           format(Sys.Date(), "%y%m%d"),
                           ifelse(is.null(filled_by), "", paste(filled_by, ".", sep = "")),
                           project))
  return(df_wide)
}

