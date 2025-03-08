#' Retrieve and Organize Table Information
#'
#' This function loads specified data tables from given file paths and returns a list of these tables.
#' If no file paths are provided, it outputs information about what each parameter expects.
#'
#' @param ASVs_Tab File path for the ASVs (Amplicon Sequence Variants) table in CSV format.
#' @param Tract_Tab File path for the sequencing QC (Quality Control) track table in CSV format.
#' @param Alpha_divTab Already loaded or computed alpha diversity table.
#' @param Alpha_div_LmerTab Already loaded or computed linear mixed-effects model results for alpha diversity.
#' @param PermanovaTab Already loaded or computed PERMANOVA (Permutational Multivariate Analysis of Variance) results.
#'
#' @return A list containing the following elements:
#' \itemize{
#'   \item{asvs}{Data frame of ASVs loaded from \code{ASVs_Tab}, if available.}
#'   \item{track}{Data frame of sequencing QC track loaded from \code{Tract_Tab}, if available.}
#'   \item{adiv}{Alpha diversity table passed through \code{Alpha_divTab}.}
#'   \item{lmer.tab}{Results of the alpha diversity linear mixed-effects model passed through \code{Alpha_div_LmerTab}.}
#'   \item{Permanova}{PERMANOVA results passed through \code{PermanovaTab}.}
#' }
#'
#' @details If no file paths are provided, the function prints out what each parameter is intended for and returns nothing.
#' This is particularly useful for debugging or when setting up the function for the first time.
#'
#' @examples
#' \dontrun{
#' tabInfo <- Go_tabInfo(
#'   ASVs_Tab = "path/to/ASVs.csv",
#'   Tract_Tab = "path/to/tract.csv",
#'   Alpha_divTab = preloaded_alpha_diversity_table,
#'   Alpha_div_LmerTab = preloaded_lmer_results,
#'   PermanovaTab = preloaded_permanova_results
#' )
#' print(tabInfo)
#' }
#'
#' @export

Go_tabInfo <- function(Taxa_Tab=NA,
                       Tract_Tab=NA,
                       Func_Tab =NA,
                       Other_Tab=NA,
                       Alpha_divTab=NA,
                       Alpha_div_LmerTab=NA,
                       RNAseq=NA,
                       HumannTab=NA,
                       Tab1=NA,
                       Tab2=NA,
                       PermanovaTab=NA) {
  # Helper function to check if all elements are NA
  all_na <- function(x) {
    all(is.na(x)) && length(x) == 1
  }

  # Check if all arguments are missing and print options if they are
  if (all_na(Taxa_Tab) && all_na(Tract_Tab) && all_na(Alpha_divTab) && all_na(Alpha_div_LmerTab) &&  all_na(HumannTab) &&
      all_na(PermanovaTab) && all_na(RNAseq) && all_na(Tab1)) {
    cat(
      "Taxa_Tab: Add the location of the ASVs table. \n",
      "Tract_Tab: Add the location of the sequencing QC tract table.\n\n",
      "Alpha_divTab: Add the calculation of the alpha diversity table. \n",
      "Alpha_div_LmerTab: Add the calculation of the alpha diversity lmer table.\n",
      "PermanovaTab: Add the calculation of the PERMANOVA table.\n")
    return(invisible())
  }

  # Using safely read.csv to handle potential read errors or empty paths
  safely_read_table <- function(path) {
    if (!is.null(path) && !all(is.na(path)) && nzchar(path)) {
      ext <- tolower(tools::file_ext(path))
      cat("[DEBUG] ext = ", ext, "\n")

      result <- tryCatch({
        if (ext %in% c("txt", "tsv", "tab")) {
          cat("[DEBUG] -> read.delim() branch\n")
          read.delim(path, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
        } else {
          cat("[DEBUG] -> read.csv() branch\n")
          read.csv(path, header = TRUE, row.names = 1, check.names = FALSE)
        }
      }, error = function(e) {
        cat("[ERROR]", e$message, "\n")
        NULL
      })
      return(result)
    } else {
      cat("[DEBUG] path is invalid\n")
      return(NULL)
    }
  }


  return(list(
    asvs = lapply(Taxa_Tab, safely_read_table),
    func = lapply(Func_Tab, safely_read_table),
    otherTab = lapply(Other_Tab, safely_read_table),
    track = safely_read_table(Tract_Tab),
    rnaseq = safely_read_table(RNAseq),
    humanntab = safely_read_table(HumannTab),
    tab1 = safely_read_table(Tab1),
    tab2 = safely_read_table(Tab2),
    adiv = lapply(Alpha_divTab, safely_read_table),
    lmer.tab = Alpha_div_LmerTab,
    Permanova = PermanovaTab
  ))
}
