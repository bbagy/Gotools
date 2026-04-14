#' Assemble Multiple GoTools Plots with patchwork
#'
#' Takes plot objects returned by any GoTools function called with
#' \code{patchwork = TRUE} and assembles them into a single figure using the
#' \pkg{patchwork} package, then saves to PDF.
#'
#' @param plots A ggplot object, or a list of ggplot objects (optionally
#'   nested, as returned by \code{Go_barchart}).
#' @param project Project name used for the output folder and file naming.
#' @param name Optional name tag for the output filename.
#' @param ncol Number of columns in the assembled layout. Passed to
#'   \code{patchwork::wrap_plots()}.
#' @param nrow Number of rows in the assembled layout.
#' @param design A patchwork layout string (e.g. \code{"AB\nCD"}).  When
#'   supplied, \code{ncol} and \code{nrow} are ignored.
#' @param guides How to collect guides: \code{"collect"} (default),
#'   \code{"keep"}, or \code{"auto"}.
#' @param height PDF height in inches.
#' @param width PDF width in inches.
#' @param tag_levels Tag panels with letters (\code{"A"}) or numbers
#'   (\code{"1"}). \code{NULL} disables tagging (default).
#' @param title Optional overall figure title.
#' @param subtitle Optional overall subtitle.
#'
#' @return Invisibly returns the assembled patchwork object.
#'
#' @examples
#' \dontrun{
#' p1 <- Go_boxplot(df = dat, cate.vars = "Group", project = "P",
#'                  outcomes = "Shannon", patchwork = TRUE)
#' p2 <- Go_volcanoPlot(project = "P", result = da_dir,
#'                      fc = 1, font = 10, height = 5, width = 5,
#'                      patchwork = TRUE)
#'
#' Go_patchwork(plots = list(p1[[1]], p2[[1]]),
#'              project = "P", name = "fig1",
#'              ncol = 2, height = 5, width = 10)
#' }
#'
#' @export
Go_patchwork <- function(plots,
                          project,
                          name      = NULL,
                          ncol      = NULL,
                          nrow      = NULL,
                          design    = NULL,
                          guides    = "collect",
                          height    = 7,
                          width     = 10,
                          tag_levels = NULL,
                          title     = NULL,
                          subtitle  = NULL) {

  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop("Package 'patchwork' is required. Install with: install.packages('patchwork')")
  }

  # ── flatten nested lists ──────────────────────────────────────────────────
  flatten_plots <- function(x) {
    if (inherits(x, "gg") || inherits(x, "ggplot") || inherits(x, "patchwork")) {
      return(list(x))
    }
    if (is.list(x)) {
      result <- list()
      for (item in x) {
        result <- c(result, flatten_plots(item))
      }
      return(result)
    }
    list()
  }

  plot_list <- flatten_plots(plots)
  plot_list <- Filter(Negate(is.null), plot_list)

  if (length(plot_list) == 0) {
    message("[Gg_patchwork] No valid ggplot objects found in `plots`.")
    return(invisible(NULL))
  }

  # ── assemble ─────────────────────────────────────────────────────────────
  pw_args <- list()
  if (!is.null(ncol))   pw_args$ncol   <- ncol
  if (!is.null(nrow))   pw_args$nrow   <- nrow
  if (!is.null(design)) pw_args$design <- design
  pw_args$guides <- guides

  combined <- do.call(patchwork::wrap_plots, c(plot_list, pw_args))

  # ── annotation ───────────────────────────────────────────────────────────
  annot_args <- list()
  if (!is.null(title))      annot_args$title      <- title
  if (!is.null(subtitle))   annot_args$subtitle   <- subtitle
  if (!is.null(tag_levels)) annot_args$tag_levels <- tag_levels

  if (length(annot_args) > 0) {
    combined <- combined + do.call(patchwork::plot_annotation, annot_args)
  }

  # ── output directory ─────────────────────────────────────────────────────
  out_dir <- file.path(sprintf("%s_%s/pdf", project, format(Sys.Date(), "%y%m%d")))
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  file_name <- sprintf("patchwork.%s.%s%s.pdf",
                       project,
                       ifelse(is.null(name), "", paste0(name, ".")),
                       format(Sys.Date(), "%y%m%d"))

  pdf_path <- file.path(out_dir, file_name)
  ggplot2::ggsave(pdf_path, plot = combined, width = width, height = height)
  message("[Gg_patchwork] Saved: ", pdf_path)

  invisible(combined)
}
