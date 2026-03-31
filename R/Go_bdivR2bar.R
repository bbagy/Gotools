#' Screen metadata variables with PERMANOVA R2 bar plots
#'
#' Runs one-variable PERMANOVA models against a beta-diversity distance matrix
#' and summarizes the explained variance (`R2`) as a faceted bar chart.
#'
#' The function is designed as a companion to `Go_bdivPM()` when there are many
#' metadata variations and a single ordination plot is no longer sufficient.
#'
#' @param psIN Phyloseq object.
#' @param project Project name used for output directories.
#' @param vars Optional metadata variables to screen. If `NULL`, variables are
#'   detected automatically from `sample_data(psIN)`.
#' @param distance_metric Distance metric passed to `phyloseq::distance()`.
#' @param category_map Optional named list mapping category labels to vectors of
#'   variable names.
#' @param name Optional extra suffix for output file names.
#' @param permutations Number of permutations used in `adonis2()`.
#' @param p_adjust P-value adjustment method.
#' @param min_n Minimum number of complete samples required per variable.
#' @param max_levels Maximum number of levels allowed for categorical variables
#'   in auto-detection mode.
#' @param include_numeric Include numeric metadata columns in auto-detection.
#' @param include_character Include factor/character/logical metadata columns in
#'   auto-detection.
#' @param height Height of the output PDF.
#' @param width Width of the output PDF.
#'
#' @return Invisibly returns a list with `results`, `plot`, and output paths.
#'   A CSV table and PDF figure are written to disk.
#'
#' @examples
#' \dontrun{
#' Go_bdivR2bar(
#'   psIN = ps,
#'   project = "DemoProject",
#'   distance_metric = "bray",
#'   category_map = list(
#'     "Clinical" = c("BMI", "TreatmentGroup"),
#'     "Demographic" = c("Age", "Sex")
#'   )
#' )
#' }
#'
#' @export
Go_bdivR2bar <- function(psIN,
                         project,
                         vars = NULL,
                         distance_metric = "bray",
                         category_map = NULL,
                         name = NULL,
                         permutations = 999,
                         p_adjust = "BH",
                         min_n = 3,
                         max_levels = 12,
                         include_numeric = TRUE,
                         include_character = TRUE,
                         height = 8,
                         width = 7) {

  if (!inherits(psIN, "phyloseq")) {
    stop("`psIN` must be a phyloseq object.")
  }

  output_dirs <- Go_path(project = project, pdf = "yes", table = "yes", path = NULL)
  out_pdf <- output_dirs$pdf
  out_tab <- output_dirs$tab
  out_perm <- file.path(output_dirs$main, "table", "perm")
  if (!dir.exists(out_perm)) dir.create(out_perm, recursive = TRUE)

  metadata <- as.data.frame(phyloseq::sample_data(psIN), stringsAsFactors = FALSE)
  if (nrow(metadata) == 0) {
    stop("`sample_data(psIN)` is empty.")
  }

  detect_variables <- function(df) {
    is_id_like <- function(x, nm) {
      non_na <- x[!is.na(x)]
      if (length(non_na) == 0) return(TRUE)
      nm_id <- grepl("(^|_)(sample|subject|study|patient|participant).*id$", nm, ignore.case = TRUE)
      unique_all <- length(unique(non_na)) == length(non_na)
      nm_id || unique_all
    }

    keep <- vapply(names(df), function(nm) {
      x <- df[[nm]]
      non_na <- x[!is.na(x)]
      if (length(non_na) < min_n) return(FALSE)
      nunique <- length(unique(non_na))
      if (nunique < 2) return(FALSE)
      if (is.numeric(x)) return(isTRUE(include_numeric) && !is_id_like(x, nm))
      if (is.logical(x) || is.factor(x) || is.character(x)) {
        if (!isTRUE(include_character)) return(FALSE)
        if (nunique > max_levels) return(FALSE)
        return(!is_id_like(as.character(x), nm))
      }
      FALSE
    }, logical(1))

    names(df)[keep]
  }

  if (is.null(vars)) {
    if (!is.null(category_map)) {
      vars <- unique(unlist(category_map, use.names = FALSE))
    } else {
      vars <- detect_variables(metadata)
    }
  }
  vars <- intersect(vars, colnames(metadata))
  if (length(vars) == 0) {
    stop("No usable metadata variables were available for screening.")
  }

  category_lookup <- stats::setNames(rep("Metadata", length(vars)), vars)
  if (!is.null(category_map)) {
    if (is.null(names(category_map)) || any(!nzchar(names(category_map)))) {
      stop("`category_map` must be a named list.")
    }
    for (cat_name in names(category_map)) {
      hit_vars <- intersect(category_map[[cat_name]], vars)
      category_lookup[hit_vars] <- cat_name
    }
  }

  full_dist <- phyloseq::distance(psIN, method = distance_metric)
  full_mat <- as.matrix(full_dist)

  analyze_one <- function(var_name) {
    df_var <- metadata[, var_name, drop = FALSE]
    keep <- stats::complete.cases(df_var)
    df_sub <- base::data.frame(metadata[keep, , drop = FALSE], stringsAsFactors = FALSE)
    if (nrow(df_sub) < min_n) return(NULL)

    x <- df_sub[[var_name]]
    if (is.logical(x) || is.factor(x) || is.character(x)) {
      x <- factor(x)
      x <- droplevels(x)
      if (nlevels(x) < 2) return(NULL)
      if (nlevels(x) > max_levels) return(NULL)
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

    sample_ids <- rownames(df_sub)
    dist_sub <- stats::as.dist(full_mat[sample_ids, sample_ids, drop = FALSE])
    model_df <- base::data.frame(value = df_sub[[var_name]], stringsAsFactors = FALSE)
    form <- stats::as.formula("dist_sub ~ value")

    res <- tryCatch(
      vegan::adonis2(form, data = model_df, permutations = permutations),
      error = function(e) NULL
    )
    if (is.null(res)) return(NULL)

    data.frame(
      variable = var_name,
      category = unname(category_lookup[[var_name]]),
      type = var_type,
      n = nrow(df_sub),
      levels = n_levels,
      R2 = as.numeric(res$R2[1]) * 100,
      pval = as.numeric(res$`Pr(>F)`[1]),
      stringsAsFactors = FALSE
    )
  }

  results <- dplyr::bind_rows(lapply(vars, analyze_one))
  if (nrow(results) == 0) {
    stop("No variables produced valid PERMANOVA results.")
  }

  results$p_adj <- stats::p.adjust(results$pval, method = p_adjust)
  results$sig <- dplyr::case_when(
    results$p_adj < 0.01 ~ "**",
    results$p_adj < 0.05 ~ "*",
    TRUE ~ ""
  )

  if (!is.null(category_map)) {
    category_levels <- names(category_map)
    other_levels <- setdiff(unique(results$category), category_levels)
    results$category <- factor(results$category, levels = c(category_levels, other_levels))
  } else {
    results$category <- factor(results$category, levels = unique(results$category))
  }

  results <- results[order(results$category, results$R2), , drop = FALSE]
  results$variable_plot <- stats::reorder(results$variable, results$R2)

  default_colors <- c(
    "#2196F3", "#F44336", "#4CAF50", "#FF9800", "#9C27B0",
    "#009688", "#795548", "#607D8B", "#E91E63", "#8BC34A"
  )
  cat_levels <- levels(results$category)
  cat_colors <- stats::setNames(rep(default_colors, length.out = length(cat_levels)), cat_levels)

  p <- ggplot2::ggplot(results, ggplot2::aes(x = R2, y = variable_plot, fill = category)) +
    ggplot2::geom_col(width = 0.7) +
    ggplot2::geom_text(ggplot2::aes(label = sig), hjust = -0.2, size = 3) +
    ggplot2::scale_fill_manual(values = cat_colors, name = NULL, drop = FALSE) +
    ggplot2::scale_x_continuous(
      labels = function(x) paste0(round(x, 1), "%"),
      expand = ggplot2::expansion(mult = c(0, 0.15))
    ) +
    ggplot2::facet_grid(category ~ ., scales = "free_y", space = "free_y") +
    ggplot2::labs(
      x = "R2 in PERMANOVA (%)",
      y = NULL,
      title = sprintf("Beta-diversity variance explained (%s)", distance_metric)
    ) +
    ggplot2::theme_classic(base_size = 10) +
    ggplot2::theme(
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold"),
      panel.grid.major.x = ggplot2::element_line(color = "grey90", linewidth = 0.3),
      legend.position = "none"
    )

  file_stub <- sprintf(
    "bdivR2bar.%s.%s%s.%s",
    distance_metric,
    project,
    ifelse(is.null(name), "", paste0(".", name)),
    format(Sys.Date(), "%y%m%d")
  )
  csv_path <- file.path(out_perm, paste0(file_stub, ".csv"))
  pdf_path <- file.path(out_pdf, paste0(file_stub, ".pdf"))

  utils::write.csv(results, csv_path, row.names = FALSE)
  ggplot2::ggsave(pdf_path, p, width = width, height = height)

  invisible(list(
    results = results,
    plot = p,
    csv = csv_path,
    pdf = pdf_path
  ))
}
