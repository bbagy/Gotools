#' Go_MRS_plot
#'
#' Plot helper for \code{Go_MRS_fit()} results.
#'
#' @description
#' Draws standardized plots from a \code{Go_MRS_fit()} object. The plot subtitle
#' automatically includes the recorded model engine, family, outcome type,
#' random-effect status, feature count, and summary metric.
#'
#' @param fit A \code{Go_MRS_fit()} result object.
#' @param plot_type One of \code{"score"}, \code{"roc"}, or \code{"coef"}.
#' @param top_n Integer; number of coefficients to display for
#'   \code{plot_type = "coef"}.
#' @param title Optional custom title.
#' @param style One of \code{"auto"}, \code{"default"}, or \code{"paper"}.
#' @param project Optional project name used for automatic output directory
#'   creation when saving.
#' @param name Optional filename tag appended to saved PDFs.
#' @param order Optional character vector controlling group order in binary
#'   score plots. Must have length 2 when provided.
#' @param mycol Optional color vector. The first color is used for the first
#'   group (or negative direction), and the second for the second group (or
#'   positive direction).
#'
#' @return A ggplot object. The returned plot may carry attributes such as
#'   \code{recommended_width}, \code{recommended_height}, \code{n_grp}, and
#'   \code{max_lbl_chars} for downstream PDF sizing. The output path is stored
#'   in \code{attr(plot, "saved_path")}.
#'   the output path is stored in \code{attr(plot, "saved_path")}.
#'
#' @param patchwork Logical. If \code{TRUE}, skip saving and return the plot object(s) for use with \code{Gg_patchwork()} or the \pkg{patchwork} package. Default \code{FALSE}.
#' @export
Go_MRS_plot <- function(fit,
                        plot_type = c("all", "score", "roc", "coef"),
                        top_n = 20,
                        title = NULL,
                        style = c("auto", "default", "paper"),
                        project = NULL,
                        name = NULL,
                        order = NULL,
                        mycol = c("#1f77b4", "#d62728"),
                        patchwork = FALSE) {

  `%||%` <- function(x, y) if (is.null(x) || length(x) == 0) y else x

  # Deduplicate repeated taxonomic name segments.
  # Handles both single-word and multi-word repetition:
  #   "Lactobacillus Lactobacillus crispatus"         → "Lactobacillus crispatus"
  #   "Rikenellaceae RC9 gut group Rikenellaceae RC9 gut group" → "Rikenellaceae RC9 gut group"
  .dedup_taxon <- function(x) {
    x <- trimws(x)
    w <- strsplit(x, "\\s+")[[1]]
    n <- length(w)
    for (k in seq_len(n %/% 2)) {
      if (identical(w[seq_len(k)], w[seq(k + 1L, 2L * k)])) {
        rest <- if (n > 2L * k) w[seq(2L * k + 1L, n)] else character(0)
        return(paste(c(w[seq_len(k)], rest), collapse = " "))
      }
    }
    x
  }

  plot_type <- if (missing(plot_type)) "all" else plot_type
  style <- match.arg(style)
  if (!inherits(fit, "Go_MRS_fit")) stop("`fit` must be a Go_MRS_fit object.")
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Package `ggplot2` is required.")
  if (length(mycol) < 2) stop("`mycol` must contain at least 2 colors.")
  mycol <- rep(mycol, length.out = 2)

  clean_tag <- function(x) {
    x <- gsub("[^A-Za-z0-9._-]+", "_", x)
    x <- gsub("_+", "_", x)
    x <- gsub("^_|_$", "", x)
    x
  }

  resolve_outdir <- function(project) {
    project_use <- project %||% "MRS"
    today <- format(Sys.Date(), "%y%m%d")
    root <- file.path(sprintf("%s_%s", clean_tag(project_use), today), "pdf", "MRS_plot")
    dir.create(root, recursive = TRUE, showWarnings = FALSE)
    root
  }

  build_filename <- function(plot_type, engine, name) {
    tag <- paste(c("MRS", plot_type, engine, clean_tag(name %||% "")), collapse = "_")
    tag <- gsub("_+", "_", tag)
    tag <- gsub("_$", "", tag)
    paste0(tag, ".pdf")
  }

  maybe_save_plot <- function(p, plot_type, engine, project, name) {
    if (isTRUE(patchwork) || is.null(project) || is.null(name)) return(p)
    out_path <- resolve_outdir(project = project)
    file_name <- build_filename(plot_type = plot_type, engine = engine, name = name)
    pdf_w <- attr(p, "recommended_width") %||% 5
    pdf_h <- attr(p, "recommended_height") %||% 5
    message(sprintf("[Go_MRS_plot] PDF: %.1f x %.1f in", pdf_w, pdf_h))
    ggplot2::ggsave(
      filename = file.path(out_path, file_name),
      plot = p,
      width = pdf_w,
      height = pdf_h
    )
    attr(p, "saved_path") <- file.path(out_path, file_name)
    p
  }

  annotate_wilcox <- function(grp_vec, val_vec) {
    sub_df <- data.frame(g = grp_vec, v = val_vec)
    sub_df <- sub_df[stats::complete.cases(sub_df), , drop = FALSE]
    if (length(unique(sub_df$g)) != 2 || nrow(sub_df) < 3) return(NULL)
    wt <- tryCatch(stats::wilcox.test(v ~ g, data = sub_df), error = function(e) NULL)
    if (is.null(wt)) return(NULL)
    p <- wt$p.value
    if (is.na(p)) return(NULL)
    lbl <- if (p < 0.001) {
      "p < 0.001"
    } else if (p < 0.05) {
      sprintf("p = %.3f", p)
    } else {
      sprintf("p = %.2f (NS)", p)
    }
    v_max  <- max(val_vec, na.rm = TRUE)
    v_min  <- min(val_vec, na.rm = TRUE)
    v_rng  <- v_max - v_min
    y_pos  <- v_max + v_rng * 0.05
    ggplot2::annotate(
      "text",
      x = 1.5,
      y = y_pos,
      label = lbl,
      size = 4,
      fontface = "italic"
    )
  }

  attach_plot_size_info <- function(p, group_labels, panel_width = 2.5, panel_height = 0.8, min_width = 3.8) {
    n_grp <- length(group_labels)
    max_lbl_chars <- max(nchar(as.character(group_labels)), na.rm = TRUE)
    label_width <- 0.9 + max(0, max_lbl_chars - 10) * 0.06
    outer_width <- 0.8
    extra_height <- if (n_grp <= 3) 0 else (n_grp - 3) * 0.22
    outer_height <- 0.9
    pdf_w <- max(min_width, min(9, panel_width + label_width + outer_width))
    pdf_h <- max(3.2, min(8, panel_height + extra_height + outer_height))
    attr(p, "recommended_width") <- pdf_w
    attr(p, "recommended_height") <- pdf_h
    attr(p, "n_grp") <- n_grp
    attr(p, "max_lbl_chars") <- max_lbl_chars
    p
  }

  subtitle <- fit$subtitle
  info <- fit$model_info
  pred <- fit$predictions
  coef_df <- fit$coefficients

  format_plot_subtitle <- function(subtitle, metrics, info) {
    subtitle <- subtitle %||% ""
    # AUC is shown as in-panel annotation — strip it from the subtitle
    subtitle_clean <- gsub("\\s*\\|\\s*AUC = [^|]+", "", subtitle)
    trimws(subtitle_clean)
  }

  # Helper: build in-panel AUC label (returns NULL when AUC unavailable)
  auc_panel_label <- function(metrics) {
    if (is.null(metrics$auc) || is.na(metrics$auc)) return(NULL)
    ci <- metrics$auc_ci
    if (!is.null(ci) && length(ci) == 3 && all(is.finite(ci))) {
      sprintf("AUC = %.3f\n(95%% CI: %.3f\u2013%.3f)", metrics$auc, ci[1], ci[3])
    } else {
      sprintf("AUC = %.3f", metrics$auc)
    }
  }

  wrap_plot_subtitle <- function(x, width = 55) {
    if (!nzchar(x %||% "")) return(x)
    parts <- strsplit(x, "\n", fixed = TRUE)[[1]]
    wrapped <- unlist(lapply(parts, function(part) {
      part <- trimws(part)
      if (!nzchar(part)) return("")
      segs <- strsplit(part, " \\| ", fixed = FALSE)[[1]]
      segs <- trimws(segs)
      if (length(segs) <= 1) return(strwrap(part, width = width))

      core_idx <- which(grepl("^(outcome =|validation =|features =|AUC =)", segs))
      if (!length(core_idx)) return(strwrap(part, width = width))

      out <- character(0)
      first_core <- min(core_idx)
      if (first_core > 1) {
        out <- c(out, paste(segs[seq_len(first_core - 1)], collapse = " | "))
      }

      outcome_idx <- which(grepl("^outcome =", segs))
      validation_idx <- which(grepl("^validation =", segs))
      features_idx <- which(grepl("^features =", segs))
      auc_idx <- which(grepl("^AUC =", segs))

      if (length(outcome_idx)) {
        end_idx <- c(validation_idx, features_idx, auc_idx)
        end_idx <- end_idx[end_idx > outcome_idx[1]]
        next_idx <- if (length(end_idx)) min(end_idx) - 1 else length(segs)
        out <- c(out, paste(segs[outcome_idx[1]:next_idx], collapse = " | "))
      }

      if (length(features_idx)) {
        end_idx <- c(auc_idx)
        end_idx <- end_idx[end_idx > features_idx[1]]
        next_idx <- if (length(end_idx)) min(end_idx) - 1 else length(segs)
        out <- c(out, paste(segs[features_idx[1]:next_idx], collapse = " | "))
      }

      if (length(auc_idx)) {
        out <- c(out, paste(segs[auc_idx[1]:length(segs)], collapse = " | "))
      }

      out[nzchar(out)]
    }))
    paste(wrapped[nzchar(wrapped)], collapse = "\n")
  }

  format_plot_title <- function(x, default) {
    out <- x %||% default
    out <- gsub("_", " ", out, fixed = TRUE)
    trimws(out)
  }

  subtitle <- format_plot_subtitle(subtitle, fit$metrics, info)
  subtitle <- wrap_plot_subtitle(subtitle, width = 52)

  if (length(plot_type) > 1 || identical(plot_type, "all")) {
    plot_types <- unique(if (identical(plot_type, "all")) c("score", "roc", "coef") else plot_type)
    out_list <- lapply(plot_types, function(pt) {
      Go_MRS_plot(
        fit = fit,
        plot_type = pt,
        top_n = top_n,
        title = title,
        style = style,
        project = project,
        name = name,
        order = order,
        mycol = mycol,
        patchwork = patchwork
      )
    })
    names(out_list) <- plot_types
    return(invisible(out_list))
  }

  plot_type <- match.arg(plot_type, c("score", "roc", "coef"))

  if (identical(style, "auto")) {
    style <- if (info$engine %in% c("cv.glmnet", "weighted_sum_glmnet", "mrs_glmnet") && identical(plot_type, "score")) "paper" else "default"
  }

  if (plot_type == "roc") {
    if (info$outcome_type != "binary") {
      stop("ROC plot is available only for binary outcomes.")
    }
    if (is.null(fit$metrics$roc)) {
      stop("ROC object not available. Install `pROC` and refit if needed.")
    }

    roc_obj <- fit$metrics$roc
    roc_df <- data.frame(
      specificity = roc_obj$specificities,
      sensitivity = roc_obj$sensitivities
    )

    p <- ggplot2::ggplot(roc_df, ggplot2::aes(x = 1 - specificity, y = sensitivity)) +
      ggplot2::geom_path(linewidth = 1, colour = "#1B9E77") +
      ggplot2::geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey60") +
      ggplot2::coord_equal() +
      ggplot2::labs(
        title = title %||% "ROC Curve",
        subtitle = subtitle,
        x = "1 - Specificity",
        y = "Sensitivity"
      ) +
      ggplot2::theme_bw(base_size = 12) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
        plot.subtitle = ggplot2::element_text(size = 9.5, hjust = 0.5, lineheight = 1.05, margin = ggplot2::margin(b = 8))
      )
    # AUC annotation inside ROC panel (bottom-right)
    auc_lbl <- auc_panel_label(fit$metrics)
    if (!is.null(auc_lbl)) {
      p <- p + ggplot2::annotate("text",
        x = 0.60, y = 0.18,
        label = auc_lbl,
        size = 3.5, hjust = 0,
        lineheight = 0.9,
        fontface = "bold"
      )
    }
    p <- attach_plot_size_info(p, group_labels = c("ROC"), panel_width = 3, panel_height = 3, min_width = 4.2)
    p <- maybe_save_plot(p, plot_type = plot_type, engine = info$engine, project = project, name = name)
    return(invisible(p))
  }

  if (plot_type == "coef") {
    if (!nrow(coef_df)) stop("No coefficients available to plot.")
    coef_sub <- utils::head(coef_df, top_n)
    coef_sub$feature <- vapply(as.character(coef_sub$feature), .dedup_taxon, character(1))
    coef_sub$feature <- factor(coef_sub$feature, levels = rev(unique(coef_sub$feature)))

    p <- ggplot2::ggplot(
      coef_sub,
      ggplot2::aes(x = estimate, y = feature, fill = estimate > 0)
    ) +
      ggplot2::geom_col(width = 0.75) +
      ggplot2::scale_fill_manual(values = c("TRUE" = mycol[2], "FALSE" = mycol[1]), guide = "none") +
      ggplot2::labs(
        title = title %||% "Feature Weights",
        subtitle = subtitle,
        x = "Estimate",
        y = NULL
      ) +
      ggplot2::theme_bw(base_size = 12) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold", hjust = 0.2),
        plot.subtitle = ggplot2::element_text(size = 9.5, hjust = 1, lineheight = 1.05, margin = ggplot2::margin(b = 8)),
        axis.text.y = ggplot2::element_text(face = "italic")
      )
    p <- attach_plot_size_info(p, group_labels = coef_sub$feature, panel_width = 2.5, panel_height = 0.8, min_width = 5.4)
    p <- maybe_save_plot(p, plot_type = plot_type, engine = info$engine, project = project, name = name)
    return(invisible(p))
  }

  if (info$outcome_type == "binary") {
    neg_lab <- info$negative_class %||% "0"
    pos_lab <- info$positive_class %||% "1"
    grp_order <- order %||% c(neg_lab, pos_lab)
    if (length(grp_order) != 2) stop("`order` must have length 2 for binary score plots.")
    pred$group <- factor(
      pred$observed,
      levels = c(0, 1),
      labels = c(neg_lab, pos_lab)
    )
    pred$group <- factor(as.character(pred$group), levels = grp_order)

    fill_vals <- stats::setNames(mycol[seq_len(2)], grp_order)
    grp_vec <- pred$group

    if (identical(style, "paper")) {
      p <- ggplot2::ggplot(
        pred,
        ggplot2::aes(x = group, y = score, fill = group)
      ) +
        ggplot2::geom_boxplot(alpha = 0.7, outlier.shape = NA, width = 0.55) +
        ggplot2::geom_jitter(width = 0.15, size = 1.5, alpha = 0.5) +
        ggplot2::scale_fill_manual(values = fill_vals, labels = grp_order, name = NULL) +
        ggplot2::labs(
          title = format_plot_title(title, "Score Distribution"),
          subtitle = subtitle,
          x = "Trajectory",
          y = "Microbial Risk Score (MRS)"
        ) +
        ggplot2::theme_classic(base_size = 12) +
        ggplot2::theme(
          plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
          plot.subtitle = ggplot2::element_text(size = 9.5, hjust = 0.1, lineheight = 1.05, margin = ggplot2::margin(b = 8)),
          axis.text.x = ggplot2::element_text(angle = 0, hjust = 0.5)
        )
    } else {
      p <- ggplot2::ggplot(
        pred,
        ggplot2::aes(x = group, y = score, fill = group)
      ) +
        ggplot2::geom_violin(alpha = 0.7, colour = NA) +
        ggplot2::geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white") +
        ggplot2::geom_jitter(width = 0.08, size = 1.5, alpha = 0.6) +
        ggplot2::scale_fill_manual(values = fill_vals) +
        ggplot2::labs(
          title = format_plot_title(title, "Score Distribution"),
          subtitle = subtitle,
          x = "Trajectory",
          y = "Model score"
        ) +
        ggplot2::theme_bw(base_size = 12)
    }

    # AUC annotation inside score panel (top-left corner)
    auc_lbl <- auc_panel_label(fit$metrics)
    if (!is.null(auc_lbl)) {
      p <- p + ggplot2::annotate("text",
        x = -Inf, y = Inf,
        label = auc_lbl,
        hjust = -0.1, vjust = 1.2,
        size = 3.2,
        lineheight = 0.9,
        fontface = "italic"
      )
    }
    p <- attach_plot_size_info(
      p,
      group_labels = grp_order,
      panel_width = 3,
      panel_height = 3,
      min_width = 3.0
    )
    p <- maybe_save_plot(p, plot_type = plot_type, engine = info$engine, project = project, name = name)
    return(invisible(p))
  }

  p <- ggplot2::ggplot(
    pred,
    ggplot2::aes(x = observed, y = score)
  ) +
    ggplot2::geom_point(alpha = 0.7, colour = mycol[1]) +
    ggplot2::geom_smooth(method = "lm", se = FALSE, colour = mycol[2], linewidth = 0.9) +
    ggplot2::labs(
      title = title %||% "Observed vs Predicted",
      subtitle = subtitle,
      x = "Observed",
      y = "Predicted score"
    ) +
    ggplot2::theme_bw(base_size = 12)

  p <- attach_plot_size_info(p, group_labels = c("Observed"), panel_width = 2.5, panel_height = 0.8, min_width = 4.4)
  p <- maybe_save_plot(p, plot_type = plot_type, engine = info$engine, project = project, name = name)
  invisible(p)
}
