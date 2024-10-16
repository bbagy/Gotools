#' Go_indexmatch: Dynamically Match and Update Data Frames
#'
#' This function performs a flexible index matching operation between two data frames (`to` and `from`) based on specified matching criteria that can be either row names or column values. It updates the `to` data frame with values from specified columns in the `from` data frame based on the match. Additionally, it generates a CSV file of the updated `to` data frame, incorporating specified naming conventions.
#'
#' @param project A character string specifying the project name, used in the naming convention of the output CSV file.
#' @param to The target data frame to be updated.
#' @param from The source data frame from which values are matched and retrieved.
#' @param matchByTo Specifies the matching criterion for the `to` data frame. This can be "rownames" to use row names for matching or the name of a column in `to`.
#' @param matchByFrom Specifies the matching criterion for the `from` data frame. This can be "rownames" to use row names for matching or the name of a column in `from`.
#' @param columns A character vector specifying which columns in the `from` data frame should be used to update the corresponding columns in the `to` data frame. This assumes the specified columns exist in both data frames.
#' @param name An optional character string to add a specific identifier to the output CSV file name.
#'
#' @return Returns the updated `to` data frame after performing the index match and column updates.
#'
#' @examples
#' # Assume df1 and df2 are existing data frames where df1 needs to be updated based on df2
#' df1 <- data.frame(id = 1:3, value = c("A", "B", "C"))
#' df2 <- data.frame(id = 1:3, value = c("X", "Y", "Z"))
#' # Update df1 based on matching column 'id'
#' updated_df <- Go_indexmatch(project = "MyProject",
#'                             to = df1,
#'                             from = df2,
#'                             matchByTo = "id",
#'                             matchByFrom = "id",
#'                             columns = c("value"),
#'                             name = "UpdatedDF")
#'
#' # Example with row names
#' rownames(df1) <- df1$id
#' rownames(df2) <- df2$id
#' updated_df <- Go_indexmatch(project = "MyProject",
#'                             to = df1,
#'                             from = df2,
#'                             matchByTo = "rownames",
#'                             matchByFrom = "rownames",
#'                             columns = c("value"),
#'                             name = "UpdatedDFUsingRowNames")
#'
#' @export
#'
#' @importFrom utils write.csv

Go_indexmatch <- function(project,
                          to,
                          from,
                          matchByTo,    # Can be "rownames" or a column name in 'to'
                          matchByFrom,  # Can be "rownames" or a column name in 'from'
                          columns,
                          name=NULL){


  to <- as.data.frame(to)
  from <- as.data.frame(from)

  # Determine the matching vector for 'to'
  matchVectorTo <- if(matchByTo == "rownames") rownames(to) else to[,matchByTo]

  # Determine the matching vector for 'from'
  matchVectorFrom <- if(matchByFrom == "rownames") rownames(from) else from[,matchByFrom]

  # Perform the index match based on the determined vectors
  matched_indices <- match(matchVectorTo, matchVectorFrom)

  # Ensure 'columns' parameter is treated correctly, assuming it's always valid column names in 'to'
  if(!is.null(columns) && length(columns) > 0) {
    for(column in columns) {
      # Make sure the column exists in 'from' before attempting to update 'to'
      if(column %in% names(from)) {
        # Update 'to' dataframe with values from 'from' based on matched indices
        # This assumes 'columns' refer to column names in both 'to' and 'from'
        to[, column] <- from[matched_indices, column, drop = FALSE]
      }
    }
  }

  # Write the updated 'to' dataframe to CSV
  write.csv(to, quote = FALSE, row.names = T,
            file=sprintf("%s.indexmatch.%s%s.csv",
                         project,
                         ifelse(is.null(name), "", paste(name, ".", sep = "")),
                         format(Sys.Date(), "%y%m%d")))
  return(to)
}
