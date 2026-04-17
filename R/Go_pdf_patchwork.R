#' Merge pre-rendered PDF files into a multi-panel figure
#'
#' Takes PDF file paths as named arguments and assembles them into a single
#' PDF figure using the same \code{design} string syntax as
#' \code{Go_patchwork()}.  Each input PDF is converted to Cairo SVG and then
#' imported as a vector grob, so panel contents remain editable instead of
#' being rasterized during assembly.
#'
#' Unlike \code{Go_patchwork()}, \code{strip_labels} is not available because
#' the PDFs are already fully rendered.
#'
#' @section Layout with \code{design}:
#' Identical syntax to \code{Go_patchwork()}. Use \code{"/"} to separate rows,
#' repeat a letter to span cells, \code{"#"} for empty cells.
#' \preformatted{
#'   "AB/CD"        2x2 grid
#'   "ABB/ACC"      A spans left column; B top-right, C bottom-right
#'   "AABCC/DDEEF"  score(2) | roc(1) | coef(2) layout per row
#'   "AB#/CDE/FGH"  2 panels on top, 3+3 on rows 2-3
#' }
#'
#' @section Panel tags (\code{tag_levels}):
#' \preformatted{
#'   "A"    ->  A   B   C  ...
#'   "(a)"  ->  (a) (b) (c) ...
#'   "(1)"  ->  (1) (2) (3) ...
#' }
#'
#' @param ... Named PDF file paths, e.g.
#'   \code{A = "Figure1.pdf", B = "Figure2.pdf"}.
#' @param project   Project name for output folder and file naming.
#' @param name      Optional tag appended to the output filename.
#' @param figure_num   Integer. Output as \code{Figure<N>.<project>...pdf}.
#' @param sfigure_num  Integer. Output as \code{Supplemental_Figure<N>...pdf}.
#'   Takes precedence over \code{figure_num}.
#' @param design    Layout string (see above). \code{NULL} = auto grid.
#' @param height    Output PDF height in inches (default 7).
#' @param width     Output PDF width  in inches (default 10).
#' @param padding   Gap between panels in inches (default 0.05).
#' @param tag_levels Panel tag string, e.g. \code{"A"}, \code{"(a)"}.
#'   \code{NULL} disables tagging (default).
#' @param tag_size  Font size for panel tags in points (default 14).
#' @param tag_face  Font face: \code{"bold"} (default), \code{"plain"}, etc.
#' @param title     Optional overall figure title.
#' @param subtitle  Optional overall subtitle.
#'
#' @details
#' This function preserves vector content by converting each input PDF page to
#' Cairo SVG via the external \command{pdftocairo} utility, then importing that
#' SVG with \pkg{grImport2}.  Text may still be converted to outlines depending
#' on the source PDF and your downstream editor, but the panel is no longer
#' embedded as a bitmap.
#'
#' @return Invisibly returns the output PDF path.
#'
#' @examples
#' \dontrun{
#' Go_pdf_patchwork(
#'   A = "CCM2_16S_260416/pdf/Figure7.CCM2_16S.OR_forest.260416.pdf",
#'   B = "CCM2_16S_260416/pdf/Figure8.CCM2_16S.MRS_treated.260416.pdf",
#'   project = "CCM2_16S", name = "combined",
#'   figure_num = 10,
#'   design = "A/B",
#'   height = 16, width = 14,
#'   tag_levels = "(a)",
#'   title = "Figure 10. Prediction Summary"
#' )
#' }
#'
#' @export
Go_pdf_patchwork <- function(
    ...,
    project,
    name        = NULL,
    figure_num  = NULL,
    sfigure_num = NULL,
    design      = NULL,
    height      = 7,
    width       = 10,
    padding     = 0.05,
    tag_levels  = NULL,
    tag_size    = 14,
    tag_face    = "bold",
    title       = NULL,
    subtitle    = NULL
) {

  # ── dependencies ────────────────────────────────────────────────────────────
  for (pkg in c("grid", "grImport2")) {
    if (!requireNamespace(pkg, quietly = TRUE))
      stop(sprintf(
        "Package '%s' is required. Install: install.packages('%s')", pkg, pkg
      ))
  }
  if (!nzchar(Sys.which("pdftocairo"))) {
    stop(
      "[Go_pdf_patchwork] External command 'pdftocairo' was not found on PATH.\n",
      "Install poppler so PDF pages can be converted to Cairo SVG without rasterization."
    )
  }

  # ── collect PDF paths ────────────────────────────────────────────────────────
  dots <- list(...)
  if (length(dots) == 0) stop("[Go_pdf_patchwork] No PDF paths provided.")
  if (is.null(names(dots)) || any(!nzchar(names(dots))))
    stop("[Go_pdf_patchwork] All arguments must be named (e.g. A = 'path.pdf').")

  panel_paths <- unlist(dots)  # named: A="p1.pdf", B="p2.pdf", ...

  # expand ~ in paths
  panel_paths <- vapply(panel_paths, path.expand, character(1))

  missing_files <- panel_paths[!file.exists(panel_paths)]
  if (length(missing_files) > 0)
    stop("[Go_pdf_patchwork] PDF files not found:\n  ",
         paste(sprintf("%s -> %s", names(missing_files), missing_files),
               collapse = "\n  "))

  # ── parse_design ─────────────────────────────────────────────────────────────
  parse_design_mat <- function(d, n) {
    if (is.null(d)) {
      nc  <- ceiling(sqrt(n))
      nr  <- ceiling(n / nc)
      lts <- c(LETTERS[seq_len(n)], rep("#", nr * nc - n))
      rows <- vapply(seq_len(nr), function(i)
        paste(lts[((i - 1) * nc + 1):(i * nc)], collapse = ""), character(1))
      d <- paste(rows, collapse = "/")
    }
    rows    <- trimws(strsplit(d, "/", fixed = TRUE)[[1]])
    rows    <- rows[nzchar(rows)]
    max_len <- max(nchar(rows))
    padded  <- vapply(rows, function(r) {
      if (nchar(r) < max_len) paste0(r, strrep("#", max_len - nchar(r))) else r
    }, character(1), USE.NAMES = FALSE)
    nr <- length(padded); nc <- max_len
    mat <- matrix(NA_character_, nrow = nr, ncol = nc)
    for (i in seq_len(nr))
      for (j in seq_len(nc)) mat[i, j] <- substr(padded[i], j, j)
    mat
  }

  # ── parse_tag ────────────────────────────────────────────────────────────────
  parse_tag <- function(s) {
    if (is.null(s)) return(NULL)
    m <- regexpr("[Aa1Ii]", s)
    if (m < 0) return(list(level = s, prefix = "", suffix = ""))
    list(level  = substr(s, m, m),
         prefix = if (m > 1) substr(s, 1, m - 1) else "",
         suffix = if (m < nchar(s)) substr(s, m + 1, nchar(s)) else "")
  }

  make_tag_labels <- function(tp, n) {
    if (is.null(tp)) return(NULL)
    base <- switch(tp$level,
      "A" = LETTERS[seq_len(n)],
      "a" = letters[seq_len(n)],
      "1" = as.character(seq_len(n)),
      "I" = as.character(utils::as.roman(seq_len(n))),
      "i" = tolower(as.character(utils::as.roman(seq_len(n)))),
      LETTERS[seq_len(n)]
    )
    paste0(tp$prefix, base, tp$suffix)
  }

  # ── build layout matrix ──────────────────────────────────────────────────────
  mat         <- parse_design_mat(design, length(panel_paths))
  nr_cells    <- nrow(mat)
  nc_cells    <- ncol(mat)
  unique_lets <- setdiff(unique(as.vector(mat)), "#")

  first_pos <- vapply(unique_lets, function(l) {
    idx <- which(mat == l, arr.ind = TRUE)
    (idx[1L, 1L] - 1L) * nc_cells + idx[1L, 2L]
  }, integer(1))
  letter_order <- unique_lets[order(first_pos)]

  if (length(letter_order) != length(panel_paths)) {
    stop(sprintf(
      "[Go_pdf_patchwork] design has %d slot(s) but %d PDF(s) given.",
      length(letter_order), length(panel_paths)
    ))
  }

  letter_to_path <- stats::setNames(as.list(panel_paths), letter_order)

  # ── read PDFs as vector grobs ───────────────────────────────────────────────
  message(sprintf("[Go_pdf_patchwork] Importing %d PDF(s) as vector graphics...",
                  length(panel_paths)))

  read_pdf_as_grob <- function(pdf_path) {
    svg_path <- tempfile(pattern = "go_pdf_patchwork_")
    on.exit(unlink(svg_path, force = TRUE), add = TRUE)
    cmd <- paste(
      shQuote(Sys.which("pdftocairo")),
      "-svg",
      shQuote(pdf_path),
      shQuote(svg_path)
    )
    exit_status <- system(cmd, intern = TRUE, ignore.stderr = FALSE)
    status_code <- attr(exit_status, "status")
    if (!is.null(status_code) && status_code != 0) {
      stop(
        "[Go_pdf_patchwork] pdftocairo failed for '", pdf_path, "':\n",
        paste(exit_status, collapse = "\n")
      )
    }
    if (!file.exists(svg_path)) {
      stop(
        "[Go_pdf_patchwork] SVG conversion did not produce an output file for '",
        pdf_path, "'."
      )
    }

    pic <- tryCatch(
      grImport2::readPicture(svg_path, warn = FALSE),
      error = function(e) {
        stop(
          "[Go_pdf_patchwork] Could not import converted SVG for '", pdf_path,
          "': ", conditionMessage(e)
        )
      }
    )
    grImport2::pictureGrob(
      pic,
      x = grid::unit(0.5, "npc"),
      y = grid::unit(0.5, "npc"),
      width = grid::unit(1, "npc"),
      height = grid::unit(1, "npc"),
      just = "centre",
      expansion = 0
    )
  }

  letter_to_grob <- lapply(letter_to_path, read_pdf_as_grob)

  # ── output path ──────────────────────────────────────────────────────────────
  date_tag  <- format(Sys.Date(), "%y%m%d")
  name_part <- if (!is.null(name) && nzchar(name)) paste0(name, ".") else ""

  file_name <- if (!is.null(sfigure_num)) {
    sprintf("Supplemental_Figure%s.%s.%s%s.pdf",
            sfigure_num, project, name_part, date_tag)
  } else if (!is.null(figure_num)) {
    sprintf("Figure%s.%s.%s%s.pdf",
            figure_num, project, name_part, date_tag)
  } else {
    sprintf("merged.%s.%s%s.pdf", project, name_part, date_tag)
  }

  out_dir  <- file.path(sprintf("%s_%s/pdf", project, date_tag))
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  pdf_path <- file.path(out_dir, file_name)

  # ── draw ─────────────────────────────────────────────────────────────────────
  tag_parsed <- parse_tag(tag_levels)
  tag_labels <- make_tag_labels(tag_parsed, length(letter_order))

  n_title_lines <- (!is.null(title)) + (!is.null(subtitle))
  title_h_frac  <- n_title_lines * 0.04

  panel_top  <- 1 - title_h_frac
  cell_h_npc <- panel_top / nr_cells
  cell_w_npc <- 1 / nc_cells
  pad        <- padding / height

  # cairo_pdf for clean vector-compatible output (Inkscape-editable)
  grDevices::cairo_pdf(pdf_path, width = width, height = height)
  on.exit(grDevices::dev.off(), add = TRUE)

  grid::grid.newpage()

  if (!is.null(title)) {
    grid::grid.text(title,
      x  = 0.5, y = 1 - title_h_frac * 0.30,
      gp = grid::gpar(fontsize = 12, fontface = "bold"))
  }
  if (!is.null(subtitle)) {
    grid::grid.text(subtitle,
      x  = 0.5, y = 1 - title_h_frac * 0.75,
      gp = grid::gpar(fontsize = 10))
  }

  tag_idx <- 1L

  for (l in letter_order) {
    pos   <- which(mat == l, arr.ind = TRUE)
    r_min <- min(pos[, 1]); r_max <- max(pos[, 1])
    c_min <- min(pos[, 2]); c_max <- max(pos[, 2])

    vp_x <- (c_min - 1) * cell_w_npc + pad
    vp_y <- panel_top - r_max * cell_h_npc + pad
    vp_w <- (c_max - c_min + 1) * cell_w_npc - 2 * pad
    vp_h <- (r_max - r_min + 1) * cell_h_npc - 2 * pad

    vp <- grid::viewport(
      x = vp_x + vp_w / 2, y = vp_y + vp_h / 2,
      width = vp_w, height = vp_h,
      just  = "centre", default.units = "npc"
    )
    grid::pushViewport(vp)

    grid::grid.draw(letter_to_grob[[l]])

    if (!is.null(tag_labels)) {
      grid::grid.text(
        tag_labels[[tag_idx]],
        x    = grid::unit(5, "points"),
        y    = grid::unit(1, "npc") - grid::unit(5, "points"),
        just = c("left", "top"),
        gp   = grid::gpar(fontsize = tag_size, fontface = tag_face)
      )
      tag_idx <- tag_idx + 1L
    }

    grid::popViewport()
  }

  message("[Go_pdf_patchwork] Saved: ", pdf_path)
  invisible(pdf_path)
}
