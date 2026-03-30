#' Generate alluvial plots for selected taxa
#'
#' Builds alluvial plots from a significant-feature table and sample metadata.
#' Each target taxon is rendered as one panel showing how feature-level signal
#' flows across the requested metadata axes.
#'
#' @param project Project name used for output directory and file naming.
#' @param SigASVs Either a CSV path or a data frame containing significant
#'   features and taxonomy columns.
#' @param map Either a CSV path or a data frame containing sample metadata.
#' @param targets.bac Character vector of taxa to plot.
#' @param target.rank Taxonomy column used to match `targets.bac`, for example
#'   `"Genus"` or `"Species"`.
#' @param outcome Primary metadata variable shown on the last axis.
#' @param column1 Optional intermediate metadata variable.
#' @param column2 Optional second intermediate metadata variable.
#' @param orders Optional level order used for categorical axes.
#' @param name Optional label appended to the output file name.
#' @param height Height of the output PDF.
#' @param width Width of the output PDF.
#' @param plotCols Number of columns for multi-panel layout.
#' @param plotRows Number of rows for multi-panel layout.
#' @param size_by Whether ribbon size should reflect summed abundance or sample
#'   count per path. One of `"abundance"` or `"count"`.
#' @param alpha Ribbon transparency.
#' @param palette Optional named color vector. Names should match feature
#'   labels. Unnamed vectors are recycled in feature order.
#'
#' @return Invisibly returns a list of ggplot objects, one per requested target.
#'   A PDF is written to the dated project `pdf/` directory.
#'
#' @examples
#' \dontrun{
#' Go_alluvialplot(
#'   project = "MyProject",
#'   SigASVs = "path_to_sig_ASVs.csv",
#'   map = "path_to_metadata.csv",
#'   targets.bac = c("Bacteroides", "Escherichia"),
#'   target.rank = "Genus",
#'   outcome = "TreatmentResponse",
#'   column1 = "Diet",
#'   column2 = "BMI_group",
#'   orders = c("Control", "Treatment"),
#'   size_by = "count"
#' )
#' }
#'
#' @export
Go_alluvialplot <- function(project,
                            SigASVs,
                            map,
                            targets.bac,
                            target.rank,
                            outcome,
                            column1 = NULL,
                            column2 = NULL,
                            orders = NULL,
                            name = NULL,
                            height = 2,
                            width = 3,
                            plotCols = 2,
                            plotRows = 1,
                            size_by = c("abundance", "count"),
                            alpha = 0.85,
                            palette = NULL) {

  if (!requireNamespace("ggalluvial", quietly = TRUE)) {
    stop("ggalluvial is required for Go_alluvialplot().")
  }
  if (!requireNamespace("reshape2", quietly = TRUE)) {
    stop("reshape2 is required for Go_alluvialplot().")
  }

  if (!is.null(dev.list())) {
    grDevices::dev.off()
  }

  read_input <- function(x) {
    if (is.character(x) && length(x) == 1) {
      return(utils::read.csv(x, row.names = 1, check.names = FALSE))
    }
    if (is.data.frame(x)) {
      return(x)
    }
    stop("Inputs must be either a CSV file path or a data frame.")
  }

  sanitize_feature_label <- function(x) {
    x <- as.character(x)
    x[is.na(x) | !nzchar(trimws(x))] <- "Unknown"
    x <- gsub("\\s+", "_", x)
    x
  }

  build_plot_title <- function(target, rank_name) {
    paste0(rank_name, ": ", target)
  }

  tab <- read_input(SigASVs)
  sampledata <- read_input(map)

  if (is.null(project) || !nzchar(project)) {
    stop("`project` must be provided.")
  }
  if (is.null(targets.bac) || length(targets.bac) == 0) {
    stop("`targets.bac` must contain at least one target taxon.")
  }
  size_by <- match.arg(size_by)
  if (is.null(target.rank) || !target.rank %in% colnames(tab)) {
    stop("`target.rank` must be a column present in `SigASVs`.")
  }

  axis_vars <- c(column1, column2, outcome)
  axis_vars <- axis_vars[!is.na(axis_vars) & nzchar(axis_vars)]
  missing_axis_cols <- setdiff(axis_vars, colnames(sampledata))
  if (length(missing_axis_cols) > 0) {
    stop("Metadata columns not found in `map`: ", paste(missing_axis_cols, collapse = ", "))
  }
  if (!outcome %in% axis_vars) {
    stop("`outcome` must be present in `map`.")
  }

  numeric_cols <- names(tab)[vapply(tab, is.numeric, logical(1))]
  sample_cols <- intersect(numeric_cols, rownames(sampledata))
  if (length(sample_cols) == 0) {
    stop("No overlapping sample columns were found between `SigASVs` and `map`.")
  }

  tab <- tab[, c(setdiff(colnames(tab), sample_cols), sample_cols), drop = FALSE]
  sampledata <- sampledata[sample_cols, , drop = FALSE]
  sampledata$SampleID <- rownames(sampledata)

  if (!"Species" %in% colnames(tab)) {
    tab$Species <- tab[[target.rank]]
  }

  tab$RowSum <- rowSums(tab[, sample_cols, drop = FALSE], na.rm = TRUE)
  tab$.feature_label <- paste0(sanitize_feature_label(tab$Species), "_", tab$RowSum)
  rownames(tab) <- make.unique(tab$.feature_label)

  output_dirs <- Go_path(project = project, pdf = "yes", table = "no", path = NULL)
  out_path <- output_dirs$pdf

  resolve_palette <- function(features, user_palette = NULL) {
    if (length(features) == 0) {
      return(stats::setNames(character(0), character(0)))
    }

    if (!is.null(user_palette)) {
      if (!is.null(names(user_palette)) && all(features %in% names(user_palette))) {
        return(user_palette[features])
      }
      if (length(user_palette) >= length(features)) {
        vals <- user_palette[seq_along(features)]
        return(stats::setNames(vals, features))
      }
    }

    if (requireNamespace("scales", quietly = TRUE)) {
      vals <- scales::hue_pal(l = 65, c = 100)(length(features))
    } else {
      vals <- grDevices::rainbow(length(features))
    }
    stats::setNames(vals, features)
  }

  plots <- list()
  rendered_targets <- character(0)

  for (target in targets.bac) {
    tab.sel <- tab[as.character(tab[[target.rank]]) == target, , drop = FALSE]
    if (nrow(tab.sel) == 0) {
      message("[Go_alluvialplot] No rows found for target: ", target)
      next
    }

    feature_labels <- rownames(tab.sel)
    taxa_tab <- data.frame(
      SampleID = sample_cols,
      t(as.matrix(tab.sel[, sample_cols, drop = FALSE])),
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
    colnames(taxa_tab)[-1] <- feature_labels
    taxa_tab <- merge(sampledata, taxa_tab, by = "SampleID", all.x = FALSE, sort = FALSE)

    taxa_melt <- reshape2::melt(
      taxa_tab,
      id.vars = c("SampleID", axis_vars),
      measure.vars = feature_labels,
      variable.name = "feature",
      value.name = "value"
    )

    taxa_melt$value <- suppressWarnings(as.numeric(taxa_melt$value))
    taxa_melt <- taxa_melt[is.finite(taxa_melt$value) & taxa_melt$value > 0, , drop = FALSE]
    if (nrow(taxa_melt) == 0) {
      message("[Go_alluvialplot] No positive values remained for target: ", target)
      next
    }

    taxa_melt <- dplyr::group_by(taxa_melt, feature)
    taxa_melt <- dplyr::mutate(taxa_melt, feature_count = dplyr::n())
    taxa_melt <- dplyr::ungroup(taxa_melt)
    taxa_melt$label <- paste0(taxa_melt$feature, " (n=", taxa_melt$feature_count, ")")

    for (var_name in axis_vars) {
      taxa_melt[[var_name]] <- as.character(taxa_melt[[var_name]])
      taxa_melt[[var_name]][is.na(taxa_melt[[var_name]]) | !nzchar(taxa_melt[[var_name]])] <- "NA"

      if (!is.null(orders)) {
        observed_levels <- unique(taxa_melt[[var_name]])
        var_levels <- c(intersect(orders, observed_levels), setdiff(observed_levels, orders))
        taxa_melt[[var_name]] <- factor(taxa_melt[[var_name]], levels = var_levels)
      } else {
        taxa_melt[[var_name]] <- factor(taxa_melt[[var_name]])
      }
    }

    group_cols <- c("label", "feature", axis_vars)
    df_plot <- dplyr::summarise(
      dplyr::group_by(taxa_melt, dplyr::across(dplyr::all_of(group_cols))),
      sample_n = dplyr::n(),
      abundance_sum = sum(value, na.rm = TRUE),
      .groups = "drop"
    )
    df_plot$weight <- if (identical(size_by, "count")) df_plot$sample_n else log1p(df_plot$abundance_sum)
    df_plot <- df_plot[df_plot$weight > 0 & is.finite(df_plot$weight), , drop = FALSE]
    if (nrow(df_plot) == 0) {
      message("[Go_alluvialplot] No aggregated paths remained for target: ", target)
      next
    }

    feature_levels <- unique(df_plot$feature)
    feature_palette <- resolve_palette(feature_levels, palette)

    axis_mapping <- list(axis1 = quote(label))
    for (i in seq_along(axis_vars)) {
      axis_mapping[[paste0("axis", i + 1)]] <- as.name(axis_vars[i])
    }
    axis_mapping$y <- quote(weight)
    axis_mapping$fill <- quote(feature)

    p <- ggplot2::ggplot(
      df_plot,
      do.call(ggplot2::aes, axis_mapping)
    ) +
      ggalluvial::geom_alluvium(width = 5 / 12, alpha = alpha) +
      ggalluvial::geom_stratum(width = 5 / 12, fill = "aliceblue", color = "black") +
      ggalluvial::stat_stratum(
        ggplot2::aes(label = after_stat(stratum)),
        geom = "text",
        size = 3
      ) +
      ggplot2::scale_fill_manual(values = feature_palette, drop = FALSE) +
      ggplot2::scale_x_discrete(
        limits = c("label", axis_vars),
        labels = c(target.rank, axis_vars),
        expand = c(0.12, 0.12)
      ) +
      ggplot2::labs(
        title = build_plot_title(target, target.rank),
        x = NULL,
        y = if (identical(size_by, "count")) "sample count" else "log1p summed abundance"
      ) +
      ggplot2::theme_void() +
      ggplot2::theme(
        legend.position = "none",
        axis.title.y = ggplot2::element_blank(),
        axis.text.y = ggplot2::element_blank(),
        axis.ticks.y = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_text(size = 9, colour = "black"),
        plot.title = ggplot2::element_text(hjust = 0.1, size = 10),
        text = ggplot2::element_text(size = 8)
      )

    plots[[length(plots) + 1]] <- p
    rendered_targets <- c(rendered_targets, target)
  }

  if (length(plots) == 0) {
    stop("No alluvial plots were generated. Check `targets.bac`, `target.rank`, and sample overlap.")
  }

  file_stub <- sprintf(
    "Alluvial.%s.%s.%s%s%s.pdf",
    project,
    outcome,
    ifelse(is.null(column1), "", paste0(column1, ".")),
    ifelse(is.null(column2), "", paste0(column2, ".")),
    ifelse(is.null(name), format(Sys.Date(), "%y%m%d"), paste0(name, ".", format(Sys.Date(), "%y%m%d")))
  )
  pdf_path <- file.path(out_path, file_stub)

  grDevices::pdf(pdf_path, height = height, width = width)
  on.exit(grDevices::dev.off(), add = TRUE)

  if (length(plots) == 1) {
    print(plots[[1]])
  } else {
    multiplot(plotlist = plots, cols = plotCols, rows = plotRows)
  }

  invisible(plots)
}
