
#' Create Directories for Project Outputs
#'
#' This function creates directories for storing various outputs of a microbiome analysis project.
#'
#' @param project The name of the project, used to name the main directory.
#' @param pdf Indicator (yes/no) for creating a directory to store PDF files.
#' @param table Indicator (yes/no) for creating a directory to store table files.
#' @param path Additional specific path for creating another directory.
#'
#' @return
#' A list of paths to the created directories, with keys corresponding to the types of outputs (pdf, table, and other specified paths).
#'
#' @details
#' The function creates a main directory for the project and, based on the input parameters, can create additional subdirectories for:
#'   - PDF files (`pdf` parameter)
#'   - Table files (`table` parameter)
#'   - Any other specified type of files or outputs (`path` parameter)
#'
#' It is useful for organizing different types of outputs in a structured manner within the project.
#'
#' @examples
#' # Example usage

Go_path <- function(project, pdf, table, path) {
  # Validate the project name
  if (is.null(project) || project == "") {
    stop("Invalid input: 'project' cannot be NULL or empty.")
  }

  # Helper function to create a directory
  createDir <- function(dirPath, dirType) {
    if (!dir.exists(dirPath)) {
      dir.create(dirPath, recursive = TRUE)
      cat(sprintf("%s directory created at: %s\n", dirType, dirPath))
    }
  }

  # Main project directory
  out <- file.path(sprintf("%s_%s", project, format(Sys.Date(), "%y%m%d")))
  createDir(out, "Main")

  # Initialize the list of directories
  dirs <- list(main = out)

  # PDF directory
  if (!is.null(pdf) && tolower(pdf) == "yes") {
    out_pdf <- file.path(out, "pdf")
    createDir(out_pdf, "PDF")
    dirs$pdf <- out_pdf
  }

  # Table directory
  if (!is.null(table) && tolower(table) == "yes") {
    out_tab <- file.path(out, "table")
    createDir(out_tab, "Table")
    dirs$tab <- out_tab
  }

  # Custom path
  if (!is.null(path)) {
    out_path <- file.path(out, path)
    createDir(out_path, "Custom")
    dirs$path <- out_path
  }

  return(dirs)
}
