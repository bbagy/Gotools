#' Screen metadata variables with MiRKAT association bar plots
#'
#' Runs one-variable MiRKAT or MiRKAT-LMM models against a beta-diversity
#' distance matrix and summarizes the association signal as a faceted bar chart.
#'
#' The function is designed as a companion to `Go_bdivMK()` when there are many
#' metadata variations and a single ordination plot is no longer sufficient.
#'
#' @param psIN Phyloseq object.
#' @param project Project name used for output directories.
#' @param vars Optional metadata variables to screen. If `NULL`, variables are
#'   detected automatically from `sample_data(psIN)` or from `category_map`.
#' @param distance_metric Distance metric passed to `phyloseq::distance()`.
#' @param category_map Optional named list mapping category labels to vectors of
#'   variable names.
#' @param cate.conf Optional covariate columns added to each screened model.
#' @param strata_var Optional subject/block ID used for MiRKAT-LMM.
#' @param name Optional extra suffix for output file names.
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
#' Go_bdivMKbar(
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
Go_bdivMKbar <- function(psIN,
                         project,
                         vars = NULL,
                         distance_metric = "bray",
                         category_map = NULL,
                         cate.conf = NULL,
                         strata_var = NULL,
                         name = NULL,
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
  if (!requireNamespace("MiRKAT", quietly = TRUE)) {
    stop("MiRKAT is required for Go_bdivMKbar().")
  }

  output_dirs <- Go_path(project = project, pdf = "yes", table = "yes", path = NULL)
  out_pdf <- output_dirs$pdf
  out_mk <- file.path(output_dirs$main, "table", "mirkat")
  if (!dir.exists(out_mk)) dir.create(out_mk, recursive = TRUE)

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
      if (nm %in% cate.conf || identical(nm, strata_var)) return(FALSE)
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
  vars <- setdiff(vars, c(cate.conf, strata_var))
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

  .is_strata_confounded <- function(id_vec, mvar_vec) {
    if (is.null(id_vec) || is.null(mvar_vec)) return(FALSE)
    id_vec <- as.character(id_vec)
    mvar_vec <- as.character(mvar_vec)
    if (length(unique(stats::na.omit(id_vec))) < 2) return(FALSE)
    all(tapply(mvar_vec, id_vec, function(x) length(unique(stats::na.omit(x)))) <= 1)
  }

  .mk_y <- function(vec) {
    f <- factor(vec)
    if (nlevels(f) < 2) return(NULL)
    as.numeric(f) - 1L
  }

  .mk_X <- function(map_df, cate.conf, mvar) {
    if (is.null(cate.conf) || length(cate.conf) == 0) return(NULL)
    cvars <- intersect(setdiff(cate.conf, c("SampleType", mvar)), names(map_df))
    if (length(cvars) == 0) return(NULL)
    out <- do.call(cbind, lapply(cvars, function(cv) {
      x <- map_df[[cv]]
      if (is.numeric(x)) x else as.numeric(factor(x))
    }))
    if (is.null(dim(out))) matrix(out, ncol = 1) else out
  }

  .run_mirkat <- function(map_df, mvar, dist_mat, strata_var, cate.conf) {
    y <- .mk_y(map_df[[mvar]])
    if (is.null(y)) return(list(pval = NA_real_, krv = NA_real_, method_tag = "MiRKAT"))

    K <- try(MiRKAT::D2K(as.matrix(dist_mat)), silent = TRUE)
    if (inherits(K, "try-error")) {
      return(list(pval = NA_real_, krv = NA_real_, method_tag = "MiRKAT"))
    }

    X_mat <- .mk_X(map_df, cate.conf, mvar)
    out_type <- if (length(unique(y)) == 2) "D" else "C"

    use_lmm <- !is.null(strata_var) && strata_var %in% names(map_df) &&
      !.is_strata_confounded(map_df[[strata_var]], map_df[[mvar]])

    if (use_lmm && !("MiRKAT_LMM" %in% getNamespaceExports("MiRKAT"))) {
      use_lmm <- FALSE
    }

    if (use_lmm) {
      tt <- try(
        res <- MiRKAT::MiRKAT_LMM(
          y = y,
          X = X_mat,
          id = as.character(map_df[[strata_var]]),
          Ks = list(K = K),
          method = "davies"
        ),
        silent = TRUE
      )
      if (inherits(tt, "try-error")) {
        return(list(pval = NA_real_, krv = NA_real_, method_tag = "MiRKAT-LMM"))
      }
      pval <- tryCatch(as.numeric(res$p_values[1]), error = function(e) NA_real_)
      return(list(pval = pval, krv = NA_real_, method_tag = "MiRKAT-LMM"))
    }

    tt <- try(
      res <- MiRKAT::MiRKAT(
        y = y,
        Ks = list(K = K),
        X = X_mat,
        out_type = out_type,
        method = "davies",
        omnibus = "permutation",
        returnKRV = TRUE,
        returnR2 = FALSE
      ),
      silent = TRUE
    )
    if (inherits(tt, "try-error")) {
      return(list(pval = NA_real_, krv = NA_real_, method_tag = "MiRKAT"))
    }

    pval <- tryCatch(as.numeric(res$p_values[1]), error = function(e) NA_real_)
    krv <- tryCatch(as.numeric(res$KRV[1]), error = function(e) NA_real_)
    list(pval = pval, krv = krv, method_tag = "MiRKAT")
  }

  full_dist <- phyloseq::distance(psIN, method = distance_metric)
  full_mat <- as.matrix(full_dist)

  analyze_one <- function(var_name) {
    model_vars <- unique(c(var_name, cate.conf, strata_var))
    model_vars <- intersect(model_vars, colnames(metadata))
    df_sub <- base::data.frame(metadata[, model_vars, drop = FALSE], stringsAsFactors = FALSE)
    keep <- stats::complete.cases(df_sub)
    df_sub <- df_sub[keep, , drop = FALSE]
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
    mk_res <- .run_mirkat(df_sub, var_name, dist_sub, strata_var, cate.conf)
    if (!is.finite(mk_res$pval)) return(NULL)

    data.frame(
      variable = var_name,
      category = unname(category_lookup[[var_name]]),
      type = var_type,
      n = nrow(df_sub),
      levels = n_levels,
      KRV = mk_res$krv,
      pval = mk_res$pval,
      method = mk_res$method_tag,
      stringsAsFactors = FALSE
    )
  }

  results <- dplyr::bind_rows(lapply(vars, analyze_one))
  if (nrow(results) == 0) {
    stop("No variables produced valid MiRKAT results.")
  }

  results$p_adj <- stats::p.adjust(results$pval, method = p_adjust)
  results$sig <- dplyr::case_when(
    results$p_adj < 0.01 ~ "**",
    results$p_adj < 0.05 ~ "*",
    TRUE ~ ""
  )

  results$score <- ifelse(is.finite(results$KRV), results$KRV, -log10(pmax(results$pval, 1e-12)))
  results$score_label <- ifelse(
    is.finite(results$KRV),
    sprintf("KRV = %.3f", results$KRV),
    sprintf("-log10(p) = %.2f", -log10(pmax(results$pval, 1e-12)))
  )

  if (!is.null(category_map)) {
    category_levels <- names(category_map)
    other_levels <- setdiff(unique(results$category), category_levels)
    results$category <- factor(results$category, levels = c(category_levels, other_levels))
  } else {
    results$category <- factor(results$category, levels = unique(results$category))
  }

  results <- results[order(results$category, results$score), , drop = FALSE]
  results$variable_plot <- stats::reorder(results$variable, results$score)

  default_colors <- c(
    "#2196F3", "#F44336", "#4CAF50", "#FF9800", "#9C27B0",
    "#009688", "#795548", "#607D8B", "#E91E63", "#8BC34A"
  )
  cat_levels <- levels(results$category)
  cat_colors <- stats::setNames(rep(default_colors, length.out = length(cat_levels)), cat_levels)

  x_title <- if (all(is.finite(results$KRV))) "KRV from MiRKAT" else "MiRKAT association score"

  p <- ggplot2::ggplot(results, ggplot2::aes(x = score, y = variable_plot, fill = category)) +
    ggplot2::geom_col(width = 0.7) +
    ggplot2::geom_text(ggplot2::aes(label = sig), hjust = -0.2, size = 3) +
    ggplot2::scale_fill_manual(values = cat_colors, name = NULL, drop = FALSE) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0, 0.15))) +
    ggplot2::facet_grid(category ~ ., scales = "free_y", space = "free_y") +
    ggplot2::labs(
      x = x_title,
      y = NULL,
      title = sprintf("Beta-diversity association screen (%s)", distance_metric)
    ) +
    ggplot2::theme_classic(base_size = 10) +
    ggplot2::theme(
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold"),
      panel.grid.major.x = ggplot2::element_line(color = "grey90", linewidth = 0.3),
      legend.position = "none"
    )

  file_stub <- sprintf(
    "bdivMKbar.%s.%s%s.%s",
    distance_metric,
    project,
    ifelse(is.null(name), "", paste0(".", name)),
    format(Sys.Date(), "%y%m%d")
  )
  csv_path <- file.path(out_mk, paste0(file_stub, ".csv"))
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
