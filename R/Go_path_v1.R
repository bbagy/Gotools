
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

Go_path <- function(project, pdf, table, path){
  # main dir
  out <- file.path(sprintf("%s_%s",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out)) dir.create(out)
  

  # main pdf
  if(is.null(pdf)){
    print("No pdf dir.")
  } else if (pdf == "yes" | pdf == "Yes"|pdf == "YES"){
    out_pdf <- file.path(sprintf("%s/pdf",out)) 
    if(!file_test("-d", out_pdf)) dir.create(out_pdf)
    
    print("pdf is in your working dir. Use /dir$pdf/ for save.")
  }

  
  
  # main table
  if(is.null(table)){
    print("No table dir.")
  } else if (table == "yes" | table == "Yes"|table == "YES"){
    out_tab <- file.path(sprintf("%s/table",out)) 
    if(!file_test("-d", out_tab)) dir.create(out_tab)

    print("table is in your working dir.Use /dir$tab/ for save.")
  }
  
  if(is.null(path)){
    print("No another dir.")
  } else if(!is.null(path)){
    out_path <- file.path(sprintf("%s/%s",out,path)) 
    if(!file_test("-d", out_path)) dir.create(out_path)
    print("path is in your working dir. Use /dir$path/ for save.")
  }


  # 한개 이상 return 하기
  
  
  functionReturningTwoValues <- function() {
    dirs <- list()
    if(is.null(pdf)){
    } else if (pdf == "yes" | pdf == "Yes"|pdf == "YES"){
      dirs$pdf <- out_pdf
    }
    if(is.null(table)){
      next
    } else if (table == "yes" | table == "Yes"|table == "YES"){
      dirs$tab <- out_tab
    }
    if(is.null(path)){
    } else if(!is.null(path)){
      dirs$path <- out_path
    }

    return(dirs) 
  }

  functionReturningTwoValues ()
  
}
