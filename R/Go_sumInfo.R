#' Summarize and Organize Text Information
#'
#' This function organizes and returns the paths or vectors of paths for various summary texts related to nucleic acid extraction and sequencing.
#' If no paths are provided, it prints out a guide to using the function.
#'
#' @param Extraction_Summary Path or vector of paths to the extraction summary text.
#' @param Taxa_overview Path or vector of paths to the taxa overview text.
#' @param Alpha_div Path or vector of paths to the alpha diversity summary text.
#' @param Beta_div Path or vector of paths to the beta diversity summary text.
#' @param DA_test Path or vector of paths to the differential abundance test summary text.
#' @param Conclusion Path or vector of paths to the conclusion text.
#'
#' @return A list containing the following elements:
#' \itemize{
#'   \item{extraction_summary}{Path or vector of paths provided for the extraction summary text.}
#'   \item{taxa_overview}{Path or vector of paths provided for the taxa overview text.}
#'   \item{alpha_div}{Path or vector of paths provided for the alpha diversity summary text.}
#'   \item{beta_div}{Path or vector of paths provided for the beta diversity summary text.}
#'   \item{dA_test}{Path or vector of paths provided for the differential abundance test summary text.}
#'   \item{conclusion}{Path or vector of paths provided for the conclusion text.}
#' }
#'
#' @details If no paths are provided, the function will print a guide explaining what each parameter expects and how to use the function.
#' This can be particularly helpful for users setting up the function for the first time.
#'
#' @examples
#' \dontrun{
#' sumInfo <- Go_sumInfo(
#'   Extraction_Summary = "text sentence",
#'   Taxa_overview = "text sentence",
#'   Alpha_div = "text sentence",
#'   Beta_div = "text sentence",
#'   DA_test = "text sentence",
#'   Conclusion = "text sentence"
#' )
#' print(summaryInfo)
#' }
#'
#' @export

Go_sumInfo <- function(Prep_overview=NA,
                       Taxa_overview=NA,
                       Alpha_div=NA,
                       Beta_div=NA,
                       Bacterial_div=NA,
                       DA_test=NA,
                       Summary=NA,
                       Conclusion=NA,
                       Useful_info=NA) {

  # Check if all arguments are missing and print options if they are
  if (missing(Prep_overview) && missing(Taxa_overview) && missing(Alpha_div) &&
      missing(Beta_div) && missing(Bacterial_div) && missing(DA_test) &&missing(Conclusion) && missing(Useful_info) && missing(Summary)  && missing(Useful_info) ) {
    cat(
      "Provide paths for the respective summary text. Each parameter expects a path or a vector of paths for its respective text.\n",
      "Available parameters are: Extraction_Summary, Taxa_overview, Alpha_div, Beta_div, DA_test, and Conclusion. \n\n",
      "Formatting options:\n",
      "Link text: [Link Text](http://example.com).\n",
      "*text* : italic. \n",
      "**text** : bold. \n",
      "***text*** : bold + italic. \n"
    )
    return(invisible())
  }

  # Return a list with all the inputs, treating them directly as paths or vectors of paths
  return(list(
    prep_overview = Prep_overview,
    taxa_overview = Taxa_overview,
    alpha_div = Alpha_div,
    beta_div = Beta_div,
    bacterial_div = Bacterial_div,
    da_test = DA_test,
    summary=Summary,
    conclusion = Conclusion,
    useful_info = Useful_info
  ))
}

