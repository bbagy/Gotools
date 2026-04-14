#' Go_OR_plot
#'
#' Plot helper for \code{Go_OR_fit()} results.
#'
#' @param fit A \code{Go_OR_fit} result object.
#' @param title Optional custom title.
#' @param order Optional character vector giving the preferred display order of
#'   classes in text annotations. Only classes present in the fitted comparison
#'   are used.
#' @param style One of \code{"default"} or \code{"paper"}.
#' @param project Optional project name used for automatic PDF output.
#' @param name Optional filename tag for the saved PDF.
#'
#' @return A ggplot object with optional saved-path attributes.
#'
#' @param patchwork Logical. If \code{TRUE}, skip saving and return the plot object(s) for use with \code{Gg_patchwork()} or the \pkg{patchwork} package. Default \code{FALSE}.
#' @export
Go_OR_plot <- function(fit,
                       title = NULL,
                       order = NULL,
                       style = c("default", "paper"),
                       project = NULL,
                       name = NULL,
                       patchwork = FALSE) {

  `%||%` <- function(x, y) if (is.null(x) || length(x) == 0) y else x

  style <- match.arg(style)
  if (!inherits(fit, "Go_OR_fit")) stop("`fit` must be a Go_OR_fit object.")
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Package `ggplot2` is required.")
  if (!is.null(order)) {
    if (!is.character(order)) {
      stop("`order` must be a character vector.")
    }
  }

  clean_tag <- function(x) {
    x <- gsub("[^A-Za-z0-9._-]+", "_", x)
    x <- gsub("_+", "_", x)
    x <- gsub("^_|_$", "", x)
    x
  }

  attach_plot_size_info <- function(p, feature_labels, panel_width = 2.5, panel_height = 0.8) {
    n_feat <- length(feature_labels)
    max_lbl_chars <- max(nchar(as.character(feature_labels)), na.rm = TRUE)
    label_width <- 1.8 + max(0, max_lbl_chars - 12) * 0.05
    outer_width <- 1.6
    extra_height <- if (n_feat <= 2) 0 else (n_feat - 2) * 0.22
    outer_height <- 1.6
    pdf_w <- max(5.8, min(11, panel_width + label_width + outer_width))
    pdf_h <- min(11, panel_height + extra_height + outer_height)
    attr(p, "recommended_width") <- pdf_w
    attr(p, "recommended_height") <- pdf_h
    p
  }

  maybe_save_plot <- function(p, project, name) {
    if (isTRUE(patchwork) || is.null(project) || is.null(name)) return(p)
    out_dirs <- Go_path(project = project, pdf = "yes", table = "no", path = NULL)
    file_name <- paste0("OR_forest_", clean_tag(name), ".pdf")
    pdf_w <- attr(p, "recommended_width") %||% 7
    pdf_h <- attr(p, "recommended_height") %||% 6
    message(sprintf("[Go_OR_plot] PDF: %.1f x %.1f in", pdf_w, pdf_h))
    ggplot2::ggsave(
      filename = file.path(out_dirs$pdf, file_name),
      plot = p,
      width = pdf_w,
      height = pdf_h
    )
    attr(p, "saved_path") <- file.path(out_dirs$pdf, file_name)
    p
  }

  format_or_subtitle <- function(subtitle_raw, info, order = NULL) {
    subtitle_raw <- subtitle_raw %||% ""
    pos_class <- info$positive_class %||% "positive"
    neg_class <- info$negative_class %||% "negative"
    available_classes <- c(neg_class, pos_class)

    if (is.null(order)) {
      display_order <- c(pos_class, neg_class)
    } else {
      display_order <- intersect(order, available_classes)
      display_order <- c(display_order, setdiff(available_classes, display_order))
    }

    outcome_txt <- paste(display_order, collapse = " vs ")
    subtitle_fmt <- gsub(
      "outcome = [^|]+",
      paste0("outcome = ", outcome_txt),
      subtitle_raw
    )
    direction_txt <- paste0("Left of 1: ", neg_class, " | Right of 1: ", pos_class)
    paste(subtitle_fmt, direction_txt, sep = " | ")
  }

  res <- fit$results
  if (is.null(res) || !nrow(res)) stop("No OR results available to plot.")
  res <- res[order(res$OR, decreasing = TRUE), , drop = FALSE]
  res$feature_label <- factor(res$feature_label, levels = rev(res$feature_label))

  has_fdr     <- "sig_fdr"     %in% colnames(res)
  has_nominal <- "sig_nominal" %in% colnames(res)
  has_padj    <- "p.adj"       %in% colnames(res)

  direction <- ifelse(res$OR >= 1, "up", "down")

  res$sig_cat <- if (has_fdr && has_nominal) {
    ifelse(res$sig_fdr     & direction == "up",   "FDR_up",
    ifelse(res$sig_nominal & direction == "up",   "Nominal_up",
    ifelse(res$sig_fdr     & direction == "down", "FDR_down",
    ifelse(res$sig_nominal & direction == "down", "Nominal_down",
    "NS"))))
  } else {
    ifelse(res$sig & direction == "up",   "Nominal_up",
    ifelse(res$sig & direction == "down", "Nominal_down", "NS"))
  }
  res$sig_cat <- factor(res$sig_cat,
    levels = c("FDR_up", "Nominal_up", "FDR_down", "Nominal_down", "NS"))

  res_plot <- res[res$sig_cat != "NS", , drop = FALSE]
  if (!nrow(res_plot)) {
    stop("No significant OR results available to plot.")
  }

  res_plot$or_label <- ifelse(
    res_plot$sig_cat %in% c("FDR_up", "FDR_down") & has_padj,
    sprintf("%.2f\n(p=%.3f, q=%.3f)", res_plot$OR, res_plot$p.value, res_plot$p.adj),
    ifelse(
      res_plot$sig_cat %in% c("Nominal_up", "Nominal_down"),
      sprintf("%.2f\n(p=%.3f)", res_plot$OR, res_plot$p.value),
      ""
    )
  )

  res_plot$effect_size <- abs(log(pmax(res_plot$OR, 1e-6)))

  x_lo <- min(res_plot$LCI[is.finite(res_plot$LCI)], na.rm = TRUE)
  x_hi <- max(res_plot$UCI[is.finite(res_plot$UCI)], na.rm = TRUE)
  x_lo <- max(x_lo * 0.7, 1e-3)
  x_hi <- x_hi * 1.5

  p <- ggplot2::ggplot(
    res_plot,
    ggplot2::aes(y = feature_label, x = OR, color = sig_cat)
  ) +
    ggplot2::geom_vline(xintercept = 1, linetype = "dashed", color = "grey60") +
    ggplot2::geom_errorbarh(
      ggplot2::aes(xmin = pmax(LCI, 1e-3), xmax = UCI),
      height = 0.22,
      linewidth = 0.7
    ) +
    ggplot2::geom_point(ggplot2::aes(size = effect_size), shape = 16) +
    ggplot2::geom_text(
      ggplot2::aes(label = or_label),
      hjust = -0.15,
      size = 2.8,
      lineheight = 0.9
    ) +
    ggplot2::scale_x_log10(
      limits = c(x_lo, x_hi),
      breaks  = c(0.1, 0.25, 0.5, 1, 2, 4, 10),
      labels  = c("0.1", "0.25", "0.5", "1", "2", "4", "10")
    ) +
    ggplot2::scale_color_manual(
      values = c(
        "FDR_up"      = "#D73027",
        "Nominal_up"  = "#F4A582",
        "FDR_down"    = "#2166AC",
        "Nominal_down"= "#92C5DE",
        "NS"          = "grey60"
      ),
      labels = c(
        "FDR_up"      = "Increased, q < 0.05",
        "Nominal_up"  = "Increased, p < 0.05",
        "FDR_down"    = "Decreased, q < 0.05",
        "Nominal_down"= "Decreased, p < 0.05"
      ),
      breaks = c("FDR_up", "Nominal_up", "FDR_down", "Nominal_down"),
      drop = TRUE,
      name = NULL
    ) +
    ggplot2::scale_size_continuous(
      range = c(2, 6),
      name  = "|log(OR)|"
    )

  subtitle_raw <- format_or_subtitle(fit$subtitle, fit$model_info, order = order)
  subtitle_2ln <- if (nzchar(subtitle_raw)) {
    parts <- strsplit(subtitle_raw, " \\| ")[[1]]
    outcome_idx <- which(grepl("^outcome =", parts))
    feature_idx <- which(grepl("^features =", parts))
    split_candidates <- c(outcome_idx, feature_idx)
    split_candidates <- split_candidates[split_candidates > 1]
    if (length(split_candidates) > 0) {
      split_idx <- min(split_candidates) - 1
    } else {
      split_idx <- ceiling(length(parts) / 2)
    }
    paste0(
      paste(parts[seq_len(split_idx)], collapse = " | "),
      "\n",
      paste(parts[seq(split_idx + 1, length(parts))], collapse = " | ")
    )
  } else { "" }

  p <- p +
    ggplot2::labs(
      title    = title %||% "Adjusted Odds Ratio",
      subtitle = subtitle_2ln,
      x        = "Adjusted OR (95% CI, log scale)",
      y        = NULL
    )

  if (identical(style, "paper")) {
    p <- p +
      ggplot2::theme_classic(base_size = 12) +
      ggplot2::theme(
        plot.title    = ggplot2::element_text(face = "bold", hjust = 0),
        plot.subtitle = ggplot2::element_text(size = 8.5, hjust = 0.2, lineheight = 1.3,
                                              margin = ggplot2::margin(b = 8)),
        axis.text.y   = ggplot2::element_text(face = "italic", size = 9),
        legend.position      = "right",
        legend.justification = "top",
        legend.text          = ggplot2::element_text(size = 8),
        legend.key.size      = ggplot2::unit(0.8, "lines")
      )
  } else {
    p <- p +
      ggplot2::theme_bw(base_size = 12) +
      ggplot2::theme(
        plot.subtitle       = ggplot2::element_text(size = 8.5, hjust = 0.2, lineheight = 1.3,
                                                    margin = ggplot2::margin(b = 8)),
        axis.text.y          = ggplot2::element_text(face = "italic", size = 9),
        legend.position      = "right",
        legend.justification = "top",
        legend.text          = ggplot2::element_text(size = 8),
        legend.key.size      = ggplot2::unit(0.8, "lines")
      )
  }

  p <- attach_plot_size_info(p, feature_labels = res_plot$feature_label)
  p <- maybe_save_plot(p, project = project, name = name)
  p
}
