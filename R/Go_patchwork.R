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
#' Rows without \code{"#"} are automatically stretched to use the full width.
#' Use explicit \code{"#"} only when you want blank space.
#'
#' \preformatted{
#'   "AB"         1 row  : [A][B]
#'   "A/B"        2 rows : [A][ ]
#'                         [B][B]    (each row uses full width)
#'   "AB/CD"      2x2   : [A][B]
#'                         [C][D]
#'   "A/BCD"      top A spans full width automatically:
#'                [A][A][A]
#'                [B][C][D]
#'   "AAA/BCD"    top A spans full width:
#'                [A][A][A]
#'                [B][C][D]
#'   "AAB/CCD"    left 2/3 vs right 1/3:
#'                [A][A][B]
#'                [C][C][D]
#'   "ABC/DEC"    right panel spans both rows:
#'                [A][B][C]
#'                [D][E][C]
#'   "ABC/DE#"    bottom-right explicit empty cell:
#'                [A][B][C]
#'                [D][E][ ]
#'   "ABCD/EFG/HIJ" lower rows use wider 3-panel layout automatically
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
#' @param col_widths Optional numeric vector of relative column widths for the
#'   parsed \code{design} grid, e.g. \code{c(1, 2, 2)}.
#' @param row_heights Optional numeric vector of relative row heights for the
#'   parsed \code{design} grid, e.g. \code{c(1.3, 0.8)}.
#' @param free_panel Logical. If \code{TRUE}, wrap each ggplot with
#'   \code{patchwork::free(type = "panel")} before assembly so panel sizing is
#'   less affected by axis/label alignment across neighbors. Default
#'   \code{FALSE}.
#' @param guides Guide collection: \code{"collect"} (default), \code{"keep"},
#'   \code{"auto"}.
#' @param height PDF height in inches (default \code{7}).
#' @param width  PDF width  in inches (default \code{10}).
#' @param pdf Logical. If \code{TRUE} (default), save the assembled figure to a
#'   PDF file. If \code{FALSE}, return the assembled object without saving.
#' @param tag_levels Panel tag string, e.g. \code{"A"}, \code{"A)"},
#'   \code{"(a)"}. \code{NULL} disables tagging (default).
#' @param strip_labels Logical. If \code{TRUE}, call
#'   \code{labs(title = NULL, subtitle = NULL)} on every panel before assembly,
#'   useful for clean patchwork figures where each panel's title/subtitle is
#'   redundant. Default \code{FALSE}.
#' @param strip_legend Logical. If \code{TRUE}, call
#'   \code{theme(legend.position = "none")} on each ggplot panel before
#'   assembly. Default \code{FALSE}.
#' @param title Optional overall figure title.
#' @param subtitle Optional overall subtitle.
#'
#' @section Nesting patchwork objects:
#' The return value is itself a patchwork object and can be passed as a panel
#' to a second \code{Go_patchwork()} call.  Patchwork assigns tag letters to
#' the top-level "slots" only — a nested patchwork counts as one slot.  If you
#' want to re-tag a composite, pass it as a single letter and let the outer
#' call handle tagging.
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
                         name         = NULL,
                         figure_num   = NULL,
                         sfigure_num  = NULL,
                         design       = NULL,
                         col_widths   = NULL,
                         row_heights  = NULL,
                         free_panel   = FALSE,
                         guides       = "collect",
                         height       = 7,
                         width        = 10,
                         pdf          = TRUE,
                         tag_levels   = NULL,
                         strip_labels = FALSE,
                         strip_legend = FALSE,
                         title        = NULL,
                         subtitle     = NULL) {

  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop("Package 'patchwork' is required. ",
         "Install with: install.packages('patchwork')")
  }
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required.")
  }

  # Raise the R expression-node stack to avoid "node stack overflow" when
  # patchwork traverses complex ggplot objects during guide collection.
  .old_expr <- getOption("expressions")
  options(expressions = max(.old_expr, 5e4L))
  on.exit(options(expressions = .old_expr), add = TRUE)

  `%||%` <- function(a, b) {
    if (is.null(a) || (length(a) == 1 && is.na(a))) b else a
  }

  gcd_int <- function(a, b) {
    a <- as.integer(abs(a))
    b <- as.integer(abs(b))
    while (b != 0L) {
      tmp <- b
      b <- a %% b
      a <- tmp
    }
    a
  }

  lcm_int <- function(x) {
    x <- as.integer(x[x > 0])
    if (!length(x)) return(1L)
    Reduce(function(a, b) as.integer(a / gcd_int(a, b) * b), x)
  }

  original_design_rows <- function(d) {
    if (is.null(d)) return(character(0))
    rows <- trimws(strsplit(d, "/", fixed = TRUE)[[1]])
    rows[nzchar(rows)]
  }

  expand_design_rows <- function(d) {
    rows <- original_design_rows(d)
    if (!length(rows)) return(character(0))

    row_lens <- nchar(rows)
    common_len <- lcm_int(row_lens)

    vapply(rows, function(r) {
      chars <- strsplit(r, "", fixed = TRUE)[[1]]
      rep_each <- common_len %/% length(chars)
      paste(rep(chars, each = rep_each), collapse = "")
    }, character(1), USE.NAMES = FALSE)
  }

  # ── parse_design: row-wise expansion unless "#" is explicit ──────────────
  # Rows without "#" are stretched to the common width automatically.
  # Explicit "#" remains a real blank area.
  parse_design <- function(d) {
    if (is.null(d)) return(NULL)
    rows <- expand_design_rows(d)
    if (!length(rows)) return(NULL)
    paste(rows, collapse = "\n")
  }

  parse_design_mat <- function(d) {
    if (is.null(d)) return(NULL)
    rows <- expand_design_rows(d)
    if (!length(rows)) return(NULL)
    mat <- matrix("#", nrow = length(rows), ncol = nchar(rows[1]))
    for (i in seq_along(rows)) {
      for (j in seq_len(nchar(rows[i]))) {
        mat[i, j] <- substr(rows[i], j, j)
      }
    }
    mat
  }

  design_letters_in_order <- function(d) {
    mat <- parse_design_mat(d)
    if (is.null(mat)) return(character(0))
    unique_lets <- setdiff(unique(as.vector(mat)), "#")
    first_pos <- vapply(unique_lets, function(l) {
      idx <- which(mat == l, arr.ind = TRUE)
      (idx[1L, 1L] - 1L) * ncol(mat) + idx[1L, 2L]
    }, integer(1))
    unique_lets[order(first_pos)]
  }

  resolve_explicit_layout <- function(design, col_widths, row_heights) {
    if (is.null(design) || (is.null(col_widths) && is.null(row_heights))) {
      return(NULL)
    }
    mat <- parse_design_mat(design)
    if (is.null(mat)) return(NULL)
    raw_rows <- original_design_rows(design)
    raw_max_cols <- if (length(raw_rows)) max(nchar(raw_rows)) else ncol(mat)
    out <- list()
    if (!is.null(col_widths)) {
      if (!is.numeric(col_widths) || any(!is.finite(col_widths)) || any(col_widths <= 0)) {
        warning("[Go_patchwork] col_widths must be positive numeric. Ignoring.")
      } else if (length(col_widths) == ncol(mat)) {
        out$widths <- as.numeric(col_widths)
      } else if (length(col_widths) == raw_max_cols && (ncol(mat) %% raw_max_cols) == 0) {
        rep_each <- ncol(mat) %/% raw_max_cols
        out$widths <- rep(as.numeric(col_widths), each = rep_each)
      } else {
        warning(sprintf(
          "[Go_patchwork] col_widths length (%d) must match design columns (%d). Ignoring.",
          length(col_widths), ncol(mat)))
      }
    }
    if (!is.null(row_heights)) {
      if (!is.numeric(row_heights) || any(!is.finite(row_heights)) || any(row_heights <= 0)) {
        warning("[Go_patchwork] row_heights must be positive numeric. Ignoring.")
      } else if (length(row_heights) != nrow(mat)) {
        warning(sprintf(
          "[Go_patchwork] row_heights length (%d) must match design rows (%d). Ignoring.",
          length(row_heights), nrow(mat)))
      } else {
        out$heights <- as.numeric(row_heights)
      }
    }
    if (!length(out)) return(NULL)
    out
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
  # Also catches ggExtraPlot (ggMarginal output) and raw gtable/grob objects.
  is_gg <- function(x) {
    inherits(x, "gg")         ||
    inherits(x, "ggplot")     ||
    inherits(x, "patchwork")  ||
    inherits(x, "ggExtraPlot")||
    inherits(x, "gtable")     ||
    inherits(x, "grob")
  }

  # ggExtraPlot / gtable / grob are not ggplot objects — wrap them so patchwork
  # can place them in a layout alongside regular ggplots.
  wrap_if_needed <- function(x) {
    if (inherits(x, "patchwork")) {
      patchwork::wrap_elements(full = grid::grid.grabExpr(print(x)))
    } else if (inherits(x, "gg") || inherits(x, "ggplot")) {
      x
    } else {
      patchwork::wrap_elements(full = x)
    }
  }

  # ── flatten: collect ggplot leaves with names ─────────────────────────────
  # .depth guard prevents infinite recursion into deeply-nested ggplot
  # internals (causes "node stack overflow" for complex PCoA/bdiv plots).
  collect_plots <- function(x, name_prefix = "", .depth = 0L) {
    if (.depth > 15L) return(list())   # safety: never recurse past ggplot internals
    if (is.null(x)) return(list())
    if (is_gg(x)) {
      nm <- if (nzchar(name_prefix)) name_prefix else NA_character_
      return(list(list(name = nm, plot = wrap_if_needed(x))))
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
        el <- tryCatch(x[[i]], error = function(e) NULL)
        result <- c(result, collect_plots(el, name_prefix = key,
                                          .depth = .depth + 1L))
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

  # ── strip title/subtitle from all panels if requested ─────────────────────
  if (isTRUE(strip_labels)) {
    plot_list <- lapply(plot_list, function(p) {
      if (inherits(p, "ggplot") && !inherits(p, "patchwork")) {
        p + ggplot2::labs(title = NULL, subtitle = NULL)
      } else {
        p
      }
    })
  }

  if (isTRUE(strip_legend)) {
    plot_list <- lapply(plot_list, function(p) {
      if (inherits(p, "ggplot") && !inherits(p, "patchwork")) {
        p + ggplot2::theme(legend.position = "none")
      } else {
        p
      }
    })
  }

  if (isTRUE(free_panel)) {
    plot_list <- lapply(plot_list, function(p) {
      if (inherits(p, "ggplot") && !inherits(p, "patchwork")) {
        patchwork::free(p, type = "panel")
      } else {
        p
      }
    })
  }

  snapshot_panel <- function(p) {
    if (inherits(p, "patchwork") || inherits(p, "ggplot") || inherits(p, "gg") ||
        inherits(p, "free_plot") || inherits(p, "wrapped_patch") || inherits(p, "patch")) {
      patchwork::wrap_elements(full = grid::grid.grabExpr(print(p)))
    } else {
      patchwork::wrap_elements(full = p)
    }
  }

  # ── assemble ─────────────────────────────────────────────────────────────
  pw_design <- parse_design(design)
  explicit_layout <- resolve_explicit_layout(design, col_widths, row_heights)
  design_letters <- design_letters_in_order(design)

  pw_args          <- list()
  pw_args$guides   <- guides
  if (!is.null(pw_design)) pw_args$design <- pw_design
  if (!is.null(explicit_layout)) {
    if (!is.null(explicit_layout$widths))  pw_args$widths  <- explicit_layout$widths
    if (!is.null(explicit_layout$heights)) pw_args$heights <- explicit_layout$heights
  }

  if (length(design_letters) > 0L && length(design_letters) == length(plot_list)) {
    names(plot_list) <- design_letters
  }

  combined_base <- tryCatch(
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
    combined <- combined_base + do.call(patchwork::plot_annotation, annot_args)
  } else {
    combined <- combined_base
  }

  date_tag  <- format(Sys.Date(), "%y%m%d")
  pdf_path <- NULL
  save_plot <- function(plot_obj) plot_obj

  if (isTRUE(pdf)) {
    # ── filename ────────────────────────────────────────────────────────────
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

    # ── save ───────────────────────────────────────────────────────────────
    out_dir <- file.path(sprintf("%s_%s/pdf", project, date_tag))
    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

    pdf_path <- file.path(out_dir, file_name)
    save_plot <- function(plot_obj) {
      ggplot2::ggsave(pdf_path, plot = plot_obj, width = width, height = height)
    }
  }

  save_snapshot_fallback <- function() {
    message("[Go_patchwork] Retrying with snapshot panels...")
    plot_list_snapshot <- lapply(plot_list, snapshot_panel)
    if (length(design_letters) > 0L && length(design_letters) == length(plot_list_snapshot)) {
      names(plot_list_snapshot) <- design_letters
    }
    snapshot_args <- pw_args
    snapshot_args$guides <- "keep"
    combined_snapshot <- do.call(patchwork::wrap_plots, c(plot_list_snapshot, snapshot_args))

    if (length(annot_args) > 0) {
      combined_snapshot_full <- combined_snapshot +
        do.call(patchwork::plot_annotation, annot_args)

      tryCatch(
        save_plot(combined_snapshot_full),
        error = function(e) {
          message("[Go_patchwork] Snapshot save with tags failed: ",
                  conditionMessage(e),
                  "\nRetrying snapshot without panel tags...")
          annot_args_fallback <- annot_args
          annot_args_fallback$tag_levels <- NULL
          annot_args_fallback$tag_prefix <- NULL
          annot_args_fallback$tag_suffix <- NULL
          annot_args_fallback <- Filter(Negate(is.null), annot_args_fallback)

          if (length(annot_args_fallback) > 0) {
            combined_snapshot_notags <- combined_snapshot +
              do.call(patchwork::plot_annotation, annot_args_fallback)
            tryCatch(
              save_plot(combined_snapshot_notags),
              error = function(e2) {
                message("[Go_patchwork] Snapshot save without tags failed: ",
                        conditionMessage(e2),
                        "\nRetrying snapshot without plot_annotation...")
                save_plot(combined_snapshot)
              }
            )
          } else {
            save_plot(combined_snapshot)
          }
        }
      )
    } else {
      save_plot(combined_snapshot)
    }
  }

  if (isTRUE(pdf)) {
    tryCatch(
      save_plot(combined),
      error = function(e) {
        if (!is.null(tag_parsed)) {
          message("[Go_patchwork] Save with tags failed: ", conditionMessage(e),
                  "\nRetrying without panel tags...")
          annot_args_fallback <- annot_args
          annot_args_fallback$tag_levels <- NULL
          annot_args_fallback$tag_prefix <- NULL
          annot_args_fallback$tag_suffix <- NULL

          combined_fallback <- combined_base
          annot_args_fallback <- Filter(Negate(is.null), annot_args_fallback)
          if (length(annot_args_fallback) > 0) {
            combined_fallback <- combined_base +
              do.call(patchwork::plot_annotation, annot_args_fallback)
          }

          tryCatch(
            save_plot(combined_fallback),
            error = function(e2) {
              message("[Go_patchwork] Save without tags failed: ", conditionMessage(e2),
                      "\nRetrying without plot_annotation...")
              tryCatch(
                save_plot(combined_base),
                error = function(e3) {
                  message("[Go_patchwork] Save without plot_annotation failed: ",
                          conditionMessage(e3))
                  save_snapshot_fallback()
                }
              )
            }
          )
          return(invisible(NULL))
        }

        if (length(annot_args) > 0) {
          message("[Go_patchwork] Save with plot_annotation failed: ", conditionMessage(e),
                  "\nRetrying without plot_annotation...")
          tryCatch(
            save_plot(combined_base),
            error = function(e2) {
              message("[Go_patchwork] Save without plot_annotation failed: ",
                      conditionMessage(e2))
              save_snapshot_fallback()
            }
          )
          return(invisible(NULL))
        }

        message("[Go_patchwork] Save failed: ", conditionMessage(e))
        save_snapshot_fallback()
      }
    )
    message("[Go_patchwork] Saved: ", pdf_path)
  }

  attr(combined, "saved_path") <- pdf_path
  invisible(combined)
}
