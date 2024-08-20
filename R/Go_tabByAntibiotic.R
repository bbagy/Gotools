#' Go_tabByAntibiotic
#'
#' This function adds columns to a data frame indicating whether each subject used specific antibiotics. The function checks the usage of antibiotics by subject ID and adds a new column for each antibiotic with a "_byID" suffix. The updated data frame can also be saved as a CSV file.
#'
#' @param df A data frame containing the data to be analyzed.
#' @param antibiotics A character vector of antibiotic names corresponding to columns in `df`.
#' @param SubjectID A character string specifying the column name in `df` that contains the subject IDs.
#' @param save_CSV A logical value indicating whether to save the updated data frame as a CSV file. Default is FALSE.
#'
#' @return A data frame with additional columns indicating antibiotic usage by subject.
#'
#' @examples
#' \dontrun{
#' df <- your_data_frame
#' antibiotics <- c("Antibiotic1", "Antibiotic2")
#' updated_df <- Go_tabByAntibiotic(df, antibiotics, SubjectID = "SubjectID")
#' }
#' @export

Go_tabByAntibiotic <- function(df, antibiotics, SubjectID, save_CSV=FALSE){
  # out dir
  out <- file.path("3_map")
  if(!file_test("-d", out)) dir.create(out)
  
  # 새로 생성된 열 이름을 저장할 리스트
  new_columns <- c()
  
  # 각 항생제에 대해 _byID 열 추가
  for (antibiotic in antibiotics) {
    # 항생제 사용 여부를 SubjectID별로 확인하고 새로운 열 생성
    new_column_name <- paste0(antibiotic, "_byID")
    df[[new_column_name]] <- ave(df[[antibiotic]], df[[SubjectID]], FUN = function(x) ifelse(any(x == "Yes"), "Yes", "No"))
    
    # 새로 생성된 열 이름을 리스트에 추가
    new_columns <- c(new_columns, new_column_name)
  }
  
  # 결과를 CSV로 저장할지 여부 확인
  if (save_CSV) {
    write.csv(df, quote = FALSE, col.names = NA,
              file=sprintf("%s/%s.%s.updated_data_with_antibiotics.csv",out,
                           format(Sys.Date(), "%y%m%d"),
                           project,sep="/"))
  }
  
  # 새로 생성된 열 이름 출력
  print("Newly created columns:")
  print(new_columns)
  
  return(df)
}




