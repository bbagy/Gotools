
#' Go_renderReport - Render the HTML summary report
#'
#' Renders the bundled 16S/MG/RNAseq Rmd report template into a single-file
#' HTML summary report, using the report objects produced by
#' \code{Go_expInfo()}, \code{Go_dirInfo()}, \code{Go_tabInfo()},
#' \code{Go_imgInfo()}, and \code{Go_sumInfo()}.
#'
#' @param family One of "16S", "MG", "RNAseq".
#' @param project Project name, used in the output HTML filename.
#' @param expInfo,dirInfo,tabInfo,imgInfo,sumInfo Report objects built with
#' the matching \code{Go_*Info()} functions.
#' @param output_dir Directory the rendered HTML is written into. Created if
#' it does not exist.
#' @param logo_path Optional path to a logo image shown at the bottom of the
#' report. \code{NULL} (the default) uses the Columbia/Uhlemann Lab logo
#' bundled with Gotools; pass a different path to use another logo.
#'
#' @return Invisibly, the path to the rendered HTML file.
#'
#' @export
Go_renderReport <- function(family,
                             project,
                             expInfo, dirInfo, tabInfo, imgInfo, sumInfo,
                             output_dir,
                             logo_path = NULL) {
  family <- match.arg(family, c("16S", "MG", "RNAseq"))
  if (is.null(logo_path)) {
    logo_path <- system.file("logo", "Columbia_CUIMC_Uhlemann4.png", package = "Gotools")
  }

  rmd_file <- switch(family,
    "16S"    = "16S/CUMC Microbiome Core Summary 16S Report_V2.Rmd",
    "MG"     = "MG/CUMC Microbiome Core Summary MG Report_V3.Rmd",
    "RNAseq" = "RNAseq/CUMC Microbiome Core Summary RNAseq Report_V1.Rmd"
  )
  input_path <- system.file("rmd", rmd_file, package = "Gotools")
  if (!nzchar(input_path)) stop("Could not find bundled report template for family: ", family)
  rmd_module_dir <- system.file("rmd", "2_modules", package = "Gotools")

  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  output_file <- sprintf("%s_%s_Summary_report.html", format(Sys.Date(), "%Y%m%d"), project)

  rmarkdown::render(
    input = input_path,
    output_format = "html_document",
    output_file = output_file,
    output_dir = output_dir,
    intermediates_dir = tempdir(),
    envir = environment()
  )

  invisible(file.path(output_dir, output_file))
}
