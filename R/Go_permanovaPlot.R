#' PERMANOVA R² bar chart with optional inset summary
#'
#' Runs one-variable PERMANOVA (adonis2) for each metadata variable against a
#' beta-diversity distance matrix, then produces a publication-ready figure
#' consisting of:
#'
#' * **Main panel** – horizontal bar chart of R² (%) per variable, faceted by
#'   category and sorted within each facet.
#' * **Inset panel** (optional) – lollipop chart of mean R² per category,
#'   embedded in the upper-right corner of the main panel via
#'   `patchwork::inset_element()`.
#'
#' A CSV table of all PERMANOVA results is also written to disk alongside the
#' PDF figure.
#'
#' @param psIN Phyloseq object.
#' @param project Project name used by `Go_path()` to build output directories.
#' @param vars Character vector of metadata variable names to screen. If `NULL`
#'   (default), variables are detected automatically from `sample_data(psIN)`.
#' @param category_map Optional named list mapping category labels (character)
#'   to vectors of variable names, e.g.
#'   `list("Clinical" = c("BMI", "Age"), "Treatment" = c("Group"))`.
#'   Variables not present in the map are labelled `"Metadata"`.
#' @param distance_metric Distance metric passed to `phyloseq::distance()`.
#'   Default `"bray"`.
#' @param name Optional extra suffix appended to output file names.
#' @param permutations Number of permutations for `vegan::adonis2()`.
#'   Default `999`.
#' @param p_adjust P-value adjustment method passed to `stats::p.adjust()`.
#'   Default `"BH"`.
#' @param min_n Minimum number of complete-case samples required to run
#'   PERMANOVA for a variable. Default `3`.
#' @param max_levels Maximum number of factor levels allowed when
#'   auto-detecting categorical variables. Default `12`.
#' @param include_numeric Include numeric metadata columns in auto-detection.
#'   Default `TRUE`.
#' @param include_character Include factor/character/logical metadata columns
#'   in auto-detection. Default `TRUE`.
#' @param cat_colors Optional named character vector of hex colours, one per
#'   category. If `NULL` (default), a built-in palette is used.
#' @param show_inset Logical. Whether to embed the per-category lollipop inset.
#'   Default `TRUE`.
#' @param inset_left,inset_bottom,inset_right,inset_top Fractional coordinates
#'   (0–1) that define the inset bounding box within the main panel.
#'   Defaults place the inset at roughly the upper-right quadrant
#'   (`left=0.58`, `bottom=0.34`, `right=0.99`, `top=0.73`).
#' @param show_strip_text Logical. Whether to display facet strip labels.
#'   Default `FALSE` (strips are hidden; category identity is conveyed by fill
#'   colour and legend).
#' @param base_size Base font size (pt) for `ggplot2::theme_classic()`.
#'   Default `10`.
#' @param height Height of the saved PDF in inches. If `NULL`, height is
#'   determined automatically from the number of variables.
#' @param width Width of the saved PDF in inches. If `NULL`, width is
#'   determined automatically from the longest variable name.
#'
#' @return Invisibly returns a named list:
#' \describe{
#'   \item{`results`}{`data.frame` with one row per variable containing
#'     `variable`, `category`, `type`, `n`, `levels`, `R2`, `pval`,
#'     `p_adj`, `sig`.}
#'   \item{`plot_main`}{`ggplot` object – the main bar chart.}
#'   \item{`plot_inset`}{`ggplot` object – the lollipop inset (`NULL` when
#'     `show_inset = FALSE`).}
#'   \item{`plot_combined`}{Final combined plot written to PDF.}
#'   \item{`csv`}{Path to the CSV results file.}
#'   \item{`pdf`}{Path to the PDF figure file.}
#' }
#'
#' @examples
#' \dontrun{
#' Go_permanovaPlot(
#'   psIN = ps,
#'   project = "DemoProject",
#'   distance_metric = "bray",
#'   category_map = list(
#'     "Demographic & Lifestyle" = c("age", "sex", "BMI", "smoking"),
#'     "Periodontal Status"      = c("periodontitis_severity", "Mean_PD"),
#'     "Cognitive Function"      = c("cognitive_status", "MMSE"),
#'     "Medical History"         = c("hypertension", "diabetes")
#'   )
#' )
#' }
#'
#' @export
Go_permanovaPlot <- function(psIN,
                             project,
                             vars             = NULL,
                             category_map     = NULL,
                             distance_metric  = "bray",
                             name             = NULL,
                             permutations     = 999,
                             p_adjust         = "BH",
                             min_n            = 3,
                             max_levels       = 12,
                             include_numeric  = TRUE,
                             include_character = TRUE,
                             cat_colors       = NULL,
                             show_inset       = TRUE,
                             inset_left       = 0.58,
                             inset_bottom     = 0.34,
                             inset_right      = 0.99,
                             inset_top        = 0.73,
                             show_strip_text  = FALSE,
                             base_size        = 10,
                             height           = NULL,
                             width            = NULL) {

  # ── input validation ────────────────────────────────────────────────────────
  if (!inherits(psIN, "phyloseq")) {
    stop("`psIN` must be a phyloseq object.")
  }

  # ── output directories ──────────────────────────────────────────────────────
  output_dirs <- Go_path(project = project, pdf = "yes", table = "yes", path = NULL)
  out_pdf  <- output_dirs$pdf
  out_perm <- file.path(output_dirs$main, "table", "perm")
  if (!dir.exists(out_perm)) dir.create(out_perm, recursive = TRUE)

  # ── metadata ─────────────────────────────────────────────────────────────────
  metadata <- as.data.frame(phyloseq::sample_data(psIN), stringsAsFactors = FALSE)
  if (nrow(metadata) == 0) stop("`sample_data(psIN)` is empty.")

  # ── variable auto-detection ──────────────────────────────────────────────────
  detect_variables <- function(df) {
    is_id_like <- function(x, nm) {
      non_na <- x[!is.na(x)]
      if (length(non_na) == 0) return(TRUE)
      nm_id    <- grepl("(^|_)(sample|subject|study|patient|participant).*id$",
                        nm, ignore.case = TRUE)
      unique_all <- length(unique(non_na)) == length(non_na)
      nm_id || unique_all
    }
    keep <- vapply(names(df), function(nm) {
      x      <- df[[nm]]
      non_na <- x[!is.na(x)]
      if (length(non_na) < min_n)       return(FALSE)
      nunique <- length(unique(non_na))
      if (nunique < 2)                   return(FALSE)
      if (is.numeric(x))
        return(isTRUE(include_numeric) && !is_id_like(x, nm))
      if (is.logical(x) || is.factor(x) || is.character(x)) {
        if (!isTRUE(include_character))  return(FALSE)
        if (nunique > max_levels)        return(FALSE)
        return(!is_id_like(as.character(x), nm))
      }
      FALSE
    }, logical(1))
    names(df)[keep]
  }

  if (is.null(vars)) {
    vars <- if (!is.null(category_map)) {
      unique(unlist(category_map, use.names = FALSE))
    } else {
      detect_variables(metadata)
    }
  }
  vars <- intersect(vars, colnames(metadata))
  if (length(vars) == 0) stop("No usable metadata variables found.")

  # ── category lookup ──────────────────────────────────────────────────────────
  category_lookup <- stats::setNames(rep("Metadata", length(vars)), vars)
  if (!is.null(category_map)) {
    if (is.null(names(category_map)) || any(!nzchar(names(category_map))))
      stop("`category_map` must be a named list.")
    for (cat_name in names(category_map)) {
      hit <- intersect(category_map[[cat_name]], vars)
      category_lookup[hit] <- cat_name
    }
  }

  # ── distance matrix ──────────────────────────────────────────────────────────
  full_dist <- phyloseq::distance(psIN, method = distance_metric)
  full_mat  <- as.matrix(full_dist)

  # ── per-variable PERMANOVA ───────────────────────────────────────────────────
  analyze_one <- function(var_name) {
    df_var <- metadata[, var_name, drop = FALSE]
    keep   <- stats::complete.cases(df_var)
    df_sub <- metadata[keep, , drop = FALSE]
    if (nrow(df_sub) < min_n) return(NULL)

    x <- df_sub[[var_name]]
    if (is.logical(x) || is.factor(x) || is.character(x)) {
      x <- droplevels(factor(x))
      if (nlevels(x) < 2 || nlevels(x) > max_levels) return(NULL)
      df_sub[[var_name]] <- x
      var_type <- "categorical"
      n_levels <- nlevels(x)
    } else if (is.numeric(x)) {
      if (length(unique(x)) < 2) return(NULL)
      var_type <- "numeric"
      n_levels <- NA_integer_
    } else {
      return(NULL)
    }

    sids     <- rownames(df_sub)
    dist_sub <- stats::as.dist(full_mat[sids, sids, drop = FALSE])
    model_df <- data.frame(value = df_sub[[var_name]], stringsAsFactors = FALSE)

    res <- tryCatch(
      vegan::adonis2(dist_sub ~ value, data = model_df,
                     permutations = permutations),
      error = function(e) NULL
    )
    if (is.null(res)) return(NULL)

    data.frame(
      variable = var_name,
      category = unname(category_lookup[[var_name]]),
      type     = var_type,
      n        = nrow(df_sub),
      levels   = n_levels,
      R2       = as.numeric(res$R2[1]) * 100,
      pval     = as.numeric(res$`Pr(>F)`[1]),
      stringsAsFactors = FALSE
    )
  }

  results <- dplyr::bind_rows(lapply(vars, analyze_one))
  if (nrow(results) == 0) stop("No variables produced valid PERMANOVA results.")

  # ── p-value adjustment & significance labels ─────────────────────────────────
  results$p_adj <- stats::p.adjust(results$pval, method = p_adjust)
  results$sig   <- dplyr::case_when(
    results$p_adj < 0.01 ~ "**",
    results$p_adj < 0.05 ~ "*",
    TRUE                 ~ ""
  )

  # ── category factor ordering ─────────────────────────────────────────────────
  if (!is.null(category_map)) {
    cat_levels <- c(names(category_map),
                    setdiff(unique(results$category), names(category_map)))
  } else {
    cat_levels <- unique(results$category)
  }
  results$category      <- factor(results$category, levels = cat_levels)
  results               <- results[order(results$category, results$R2), ]
  results$variable_plot <- stats::reorder(results$variable, results$R2)

  # ── auto figure sizing ──────────────────────────────────────────────────────
  if (is.null(height)) {
    n_vars <- nrow(results)
    panel_height <- 0.8
    extra_height <- if (n_vars <= 3) 0 else (n_vars - 3) * 0.45
    outer_height <- 1.1
    height <- panel_height + extra_height + outer_height
    if (isTRUE(show_inset)) height <- height + 0.3
    height <- min(height, 12)
  }

  if (is.null(width)) {
    max_label_chars <- max(nchar(as.character(results$variable)), na.rm = TRUE)
    panel_width <- 2.5
    label_width <- 1.0 + max(0, max_label_chars - 12) * 0.06
    outer_width <- 0.8
    width <- panel_width + label_width + outer_width
    if (isTRUE(show_inset)) width <- width + 0.4
    width <- max(3, min(width, 10))
  }

  # ── colour palette ───────────────────────────────────────────────────────────
  default_palette <- c(
    "#2196F3", "#F44336", "#4CAF50", "#FF9800", "#9C27B0",
    "#009688", "#795548", "#607D8B", "#E91E63", "#8BC34A"
  )
  if (is.null(cat_colors)) {
    cat_colors <- stats::setNames(
      rep(default_palette, length.out = length(cat_levels)),
      cat_levels
    )
  } else {
    missing_cats <- setdiff(cat_levels, names(cat_colors))
    if (length(missing_cats)) {
      extra <- stats::setNames(
        rep(default_palette, length.out = length(missing_cats)),
        missing_cats
      )
      cat_colors <- c(cat_colors, extra)
    }
  }

  # ── main bar chart ───────────────────────────────────────────────────────────
  strip_text_theme <- if (isTRUE(show_strip_text)) {
    ggplot2::element_text(face = "bold", size = base_size * 0.85)
  } else {
    ggplot2::element_blank()
  }

  p_main <- ggplot2::ggplot(
    results,
    ggplot2::aes(x = R2, y = variable_plot, fill = category)
  ) +
    ggplot2::geom_col(width = 0.7) +
    ggplot2::geom_text(ggplot2::aes(label = sig), hjust = -0.2, size = 3) +
    ggplot2::scale_fill_manual(values = cat_colors, name = NULL, drop = FALSE) +
    ggplot2::scale_x_continuous(
      labels = function(x) paste0(round(x, 1), "%"),
      expand = ggplot2::expansion(mult = c(0, 0.15))
    ) +
    ggplot2::facet_grid(category ~ ., scales = "free_y", space = "free_y") +
    ggplot2::labs(
      x     = sprintf("R\u00b2 in PERMANOVA (%%) [%s]", distance_metric),
      y     = NULL,
      title = "Beta-diversity variance explained"
    ) +
    ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      strip.background   = ggplot2::element_blank(),
      strip.text         = strip_text_theme,
      panel.grid.major.x = ggplot2::element_line(color = "grey90", linewidth = 0.3),
      legend.position    = if (isTRUE(show_strip_text)) "none" else "right"
    )

  # ── inset lollipop ───────────────────────────────────────────────────────────
  p_inset  <- NULL
  p_combined <- p_main

  if (isTRUE(show_inset)) {
    overall <- results |>
      dplyr::group_by(category) |>
      dplyr::summarise(mean_R2 = mean(R2, na.rm = TRUE), .groups = "drop") |>
      dplyr::mutate(category = factor(category, levels = rev(cat_levels)))

    p_inset <- ggplot2::ggplot(
      overall,
      ggplot2::aes(x = mean_R2, y = category, color = category)
    ) +
      ggplot2::geom_segment(
        ggplot2::aes(x = 0, xend = mean_R2, y = category, yend = category),
        linewidth = 0.8, alpha = 0.9
      ) +
      ggplot2::geom_point(size = 3.5) +
      ggplot2::scale_color_manual(values = cat_colors, guide = "none") +
      ggplot2::scale_x_continuous(
        labels = function(x) paste0(round(x, 1), "%"),
        expand = ggplot2::expansion(mult = c(0.02, 0.10))
      ) +
      ggplot2::labs(x = NULL, y = NULL, title = "mean R\u00b2") +
      ggplot2::theme_classic(base_size = base_size * 0.8) +
      ggplot2::theme(
        plot.title       = ggplot2::element_text(
                             size = base_size * 0.9, hjust = 0.5, face = "bold"),
        axis.text.y      = ggplot2::element_text(size = base_size * 0.7),
        axis.text.x      = ggplot2::element_text(size = base_size * 0.7),
        axis.ticks.y     = ggplot2::element_blank(),
        panel.grid.major.x = ggplot2::element_line(color = "grey92", linewidth = 0.3),
        panel.border     = ggplot2::element_rect(color = "grey60", fill = NA,
                                                  linewidth = 0.5)
      )

    p_combined <- p_main +
      patchwork::inset_element(
        p_inset,
        left   = inset_left,
        bottom = inset_bottom,
        right  = inset_right,
        top    = inset_top
      )
  }

  # ── save outputs ─────────────────────────────────────────────────────────────
  file_stub <- sprintf(
    "permanovaPlot.%s.%s%s.%s",
    distance_metric,
    project,
    ifelse(is.null(name), "", paste0(".", name)),
    format(Sys.Date(), "%y%m%d")
  )
  csv_path <- file.path(out_perm, paste0(file_stub, ".csv"))
  pdf_path <- file.path(out_pdf,  paste0(file_stub, ".pdf"))

  utils::write.csv(results, csv_path, row.names = FALSE)
  ggplot2::ggsave(pdf_path, p_combined, width = width, height = height)

  invisible(list(
    results       = results,
    plot_main     = p_main,
    plot_inset    = p_inset,
    plot_combined = p_combined,
    csv           = csv_path,
    pdf           = pdf_path
  ))
}
