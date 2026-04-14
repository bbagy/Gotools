#' Assemble Multiple GoTools Plots with patchwork
#'
#' Takes plot objects returned by any GoTools function called with
#' \code{patchwork = TRUE} and assembles them into a single figure using the
#' \pkg{patchwork} package, then saves to PDF.
#'
#' Plots can be passed as named arguments directly (preferred) or as a single
#' named list via the \code{plots} argument for backward compatibility.
#'
#' @section Layout with \code{design}:
#' Use \code{"/"} to separate rows. Each letter refers to the plot in that
#' position (assigned A, B, C ... in the order the plots are passed).
#' \code{"#"} is an empty cell. Repeat a letter to span cells.
#' Short rows are automatically padded by repeating their last character.
#'
#' \preformatted{
#'   "AB"        1 row   : [A][B]
#'   "A/B"       2 rows  : [A]
#'                         [B]
#'   "AB/CD"     2x2     : [A][B]
#'                         [C][D]
#'   "ABC/DEC"   right panel spans both rows:
#'               [A][B][C]
#'               [D][E][C]
#'   "ABC/DE#"   bottom-right empty:
#'               [A][B][C]
#'               [D][E][ ]
#'   "ABC/DEF/H" last row auto-pads H to full width:
#'               [A][B][C]
#'               [D][E][F]
#'               [H][H][H]
#'   "AB/CDE"    top row pads B to fill:
#'               [A][B][B]
#'               [C][D][E]
#' }
#' When \code{design} is \code{NULL}, patchwork arranges plots automatically.
#'
#' @section Panel tags (\code{tag_levels}):
#' Pass a single string encoding the base alphabet plus surrounding punctuation.
#' Valid bases: \code{A}, \code{a}, \code{1}, \code{I}, \code{i}.
#' \preformatted{
#'   "A"    ->  A   B   C  ...
#'   "A)"   ->  A)  B)  C) ...
#'   "(A)"  ->  (A) (B) (C) ...
#'   "(a)"  ->  (a) (b) (c) ...
#'   "(1)"  ->  (1) (2) (3) ...
#' }
#'
#' @param ... Named ggplot objects, e.g.
#'   \code{box = p_box[[1]], roc = p_mrs$roc}.
#'   Or pass \code{plots = list(...)} for backward compatibility.
#'   Nested lists are flattened automatically.
#' @param project Project name for output folder and file naming.
#' @param name Optional name tag appended to the filename.
#' @param figure_num Integer. When supplied, file is named
#'   \code{Figure<N>.<project>.<name>.<date>.pdf}.
#' @param sfigure_num Integer. When supplied, file is named
#'   \code{Supplemental_Figure<N>.<project>.<name>.<date>.pdf}.
#'   Takes precedence over \code{figure_num}.
#' @param design Layout string using \code{"/"} as row separator (see section
#'   above). When \code{NULL}, patchwork auto-arranges.
#' @param guides Guide collection: \code{"collect"} (default), \code{"keep"},
#'   \code{"auto"}.
#' @param height PDF height in inches (default \code{7}).
#' @param width  PDF width  in inches (default \code{10}).
#' @param tag_levels Panel tag string, e.g. \code{"A"}, \code{"A)"},
#'   \code{"(a)"}. \code{NULL} disables tagging (default).
#' @param title Optional overall figure title.
#' @param subtitle Optional overall subtitle.
#'
#' @return Invisibly returns the assembled patchwork object.
#'   \code{attr(result, "saved_path")} holds the PDF path.
#'
#' @examples
#' \dontrun{
#' p_box <- Go_boxplot(df, cate.vars = "Group", project = "P",
#'                     outcomes = "Shannon", patchwork = TRUE)
#' p_mrs <- Go_MRS_plot(fit, patchwork = TRUE)  # $score $roc $coef
#'
#' # 3 plots side by side — name args as A, B, C to match design string
#' Go_patchwork(
#'   A = p_box[[1]],
#'   B = p_mrs$score,
#'   C = p_mrs$roc,
#'   project = "P", name = "overview", figure_num = 1,
#'   design = "ABC", height = 5, width = 14,
#'   tag_levels = "A)"
#' )
#'
#' # score+coef on top, roc spanning full bottom
#' Go_patchwork(
#'   A = p_mrs$score,
#'   B = p_mrs$coef,
#'   C = p_mrs$roc,
#'   project = "P", sfigure_num = 1,
#'   design = "AB/C",
#'   tag_levels = "(A)"
#' )
#' }
#'
#' @export
Go_patchwork <- function(...,
                         project,
                         name        = NULL,
                         figure_num  = NULL,
                         sfigure_num = NULL,
                         design      = NULL,
                         guides      = "collect",
                         height      = 7,
                         width       = 10,
                         tag_levels  = NULL,
                         title       = NULL,
                         subtitle    = NULL) {

  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop("Package 'patchwork' is required. ",
         "Install with: install.packages('patchwork')")
  }
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required.")
  }

  `%||%` <- function(a, b) {
    if (is.null(a) || (length(a) == 1 && is.na(a))) b else a
  }

  # ── parse_design: "ABC/DE#" → "ABC\nDE#" with auto-padding ───────────────
  # Split on "/", find max row width, pad each short row by repeating its
  # last character so all rows reach the same width.
  #   "ABC/DEF/H"  →  "ABC\nDEF\nHHH"
  #   "ABC/DE#"    →  "ABC\nDE#"       (# pads to ###? No — # explicit stop)
  #   "AB/CDE"     →  "ABB\nCDE"
  parse_design <- function(d) {
    if (is.null(d)) return(NULL)
    rows    <- trimws(strsplit(d, "/", fixed = TRUE)[[1]])
    rows    <- rows[nzchar(rows)]
    max_len <- max(nchar(rows))
    padded  <- vapply(rows, function(r) {
      n <- nchar(r)
      if (n < max_len) {
        last <- substr(r, n, n)
        r    <- paste0(r, strrep(last, max_len - n))
      }
      r
    }, character(1), USE.NAMES = FALSE)
    paste(padded, collapse = "\n")
  }

  # ── parse_tag: "(A)" → level="A" prefix="(" suffix=")" ───────────────────
  parse_tag <- function(s) {
    if (is.null(s)) return(NULL)
    m <- regexpr("[Aa1Ii]", s)
    if (m < 0) return(list(level = s, prefix = "", suffix = ""))
    list(
      level  = substr(s, m, m),
      prefix = if (m > 1) substr(s, 1, m - 1) else "",
      suffix = if (m < nchar(s)) substr(s, m + 1, nchar(s)) else ""
    )
  }

  # ── resolve plot inputs ───────────────────────────────────────────────────
  dots <- list(...)
  if (length(dots) == 1 &&
      !is.null(names(dots)) && names(dots) == "plots" &&
      is.list(dots[[1]]) && !inherits(dots[[1]], "gg")) {
    input_list <- dots[[1]]
  } else {
    input_list <- dots
  }

  # ── is_gg ─────────────────────────────────────────────────────────────────
  is_gg <- function(x) {
    inherits(x, "gg") || inherits(x, "ggplot") || inherits(x, "patchwork")
  }

  # ── flatten: collect ggplot leaves with names ─────────────────────────────
  collect_plots <- function(x, name_prefix = "") {
    if (is.null(x)) return(list())
    if (is_gg(x)) {
      nm <- if (nzchar(name_prefix)) name_prefix else NA_character_
      return(list(list(name = nm, plot = x)))
    }
    if (is.list(x)) {
      result <- list()
      nms <- names(x)
      for (i in seq_along(x)) {
        child_nm <- if (!is.null(nms) && nzchar(nms[i])) {
          nms[i]
        } else {
          as.character(i)
        }
        key <- if (nzchar(name_prefix)) {
          paste(name_prefix, child_nm, sep = ".")
        } else {
          child_nm
        }
        result <- c(result, collect_plots(x[[i]], name_prefix = key))
      }
      return(result)
    }
    list()
  }

  entries    <- collect_plots(input_list)
  plot_list  <- lapply(entries, `[[`, "plot")
  plot_names <- vapply(entries, function(e) e$name %||% "", character(1))
  plot_list  <- Filter(Negate(is.null), plot_list)

  if (length(plot_list) == 0) {
    message("[Go_patchwork] No valid ggplot objects found in inputs.")
    return(invisible(NULL))
  }

  message(sprintf("[Go_patchwork] Assembling %d plot(s): %s",
                  length(plot_list),
                  paste(plot_names[nzchar(plot_names)], collapse = ", ")))

  # ── assemble ─────────────────────────────────────────────────────────────
  pw_design <- parse_design(design)

  pw_args          <- list()
  pw_args$guides   <- guides
  if (!is.null(pw_design)) pw_args$design <- pw_design

  combined <- tryCatch(
    do.call(patchwork::wrap_plots, c(plot_list, pw_args)),
    error = function(e) {
      message("[Go_patchwork] wrap_plots failed: ", conditionMessage(e),
              "\nRetrying with guides = 'keep'...")
      pw_args$guides <- "keep"
      tryCatch(
        do.call(patchwork::wrap_plots, c(plot_list, pw_args)),
        error = function(e2) {
          stop("[Go_patchwork] Assembly failed: ", conditionMessage(e2))
        }
      )
    }
  )

  # ── annotation ───────────────────────────────────────────────────────────
  tag_parsed <- parse_tag(tag_levels)

  annot_args <- list()
  if (!is.null(title))      annot_args$title      <- title
  if (!is.null(subtitle))   annot_args$subtitle   <- subtitle
  if (!is.null(tag_parsed)) {
    annot_args$tag_levels <- tag_parsed$level
    if (nzchar(tag_parsed$prefix)) annot_args$tag_prefix <- tag_parsed$prefix
    if (nzchar(tag_parsed$suffix)) annot_args$tag_suffix <- tag_parsed$suffix
  }

  if (length(annot_args) > 0) {
    combined <- combined + do.call(patchwork::plot_annotation, annot_args)
  }

  # ── filename ──────────────────────────────────────────────────────────────
  date_tag  <- format(Sys.Date(), "%y%m%d")
  name_part <- if (!is.null(name) && nzchar(name)) paste0(name, ".") else ""

  if (!is.null(sfigure_num)) {
    file_name <- sprintf("Supplemental_Figure%s.%s.%s%s.pdf",
                         sfigure_num, project, name_part, date_tag)
  } else if (!is.null(figure_num)) {
    file_name <- sprintf("Figure%s.%s.%s%s.pdf",
                         figure_num, project, name_part, date_tag)
  } else {
    file_name <- sprintf("patchwork.%s.%s%s.pdf",
                         project, name_part, date_tag)
  }

  # ── save ─────────────────────────────────────────────────────────────────
  out_dir <- file.path(sprintf("%s_%s/pdf", project, date_tag))
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  pdf_path <- file.path(out_dir, file_name)
  ggplot2::ggsave(pdf_path, plot = combined, width = width, height = height)
  message("[Go_patchwork] Saved: ", pdf_path)

  attr(combined, "saved_path") <- pdf_path
  invisible(combined)
}
