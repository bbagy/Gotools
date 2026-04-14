#' Generate alluvial plots for feature or transition data
#'
#' Builds alluvial plots in two modes with package-aware fallbacks:
#' \itemize{
#'   \item \code{mode = "feature"} uses the original SigASVs + metadata
#'   workflow for taxa-centered alluvial plots. When the \code{alluvial}
#'   package is available, a classic alluvial layout is used; otherwise the
#'   function falls back to \code{ggalluvial}.
#'   \item \code{mode = "transition"} draws state-transition alluvial plots
#'   directly from a wide-format table such as \code{V1 -> V2 -> V3}. When the
#'   \code{alluvial} package is available, a classic alluvial layout is used;
#'   otherwise the function falls back to \code{ggalluvial}.
#' }
#'
#' @param project Project name used for output directory and file naming.
#' @param SigASVs Either a CSV path or a data frame containing significant
#'   features and taxonomy columns. Used in \code{mode = "feature"}.
#' @param map Either a CSV path or a data frame containing sample metadata.
#'   Used in \code{mode = "feature"}.
#' @param targets.bac Character vector of taxa to plot. Used in
#'   \code{mode = "feature"}.
#' @param target.rank Taxonomy column used to match \code{targets.bac}, for
#'   example \code{"Genus"} or \code{"Species"}. Used in
#'   \code{mode = "feature"}.
#' @param outcome Primary metadata variable shown on the last axis in
#'   \code{mode = "feature"}.
#' @param column1 Optional intermediate metadata variable in
#'   \code{mode = "feature"}.
#' @param column2 Optional second intermediate metadata variable in
#'   \code{mode = "feature"}.
#' @param Addcolumn Optional extra metadata columns. When provided in
#'   \code{mode = "feature"}, these are inserted before \code{outcome}.
#' @param data Either a CSV path or a data frame containing wide-format
#'   transition data. Used in \code{mode = "transition"}.
#' @param mode One of \code{"feature"} or \code{"transition"}.
#' @param axes Character vector of ordered axis columns for
#'   \code{mode = "transition"}, for example \code{c("V1","V2","V3")}.
#' @param fill_var Optional column used for ribbon fill in
#'   \code{mode = "transition"}.
#' @param facet_var Optional faceting column in \code{mode = "transition"}.
#' @param weight_var Optional weight column in \code{mode = "transition"}.
#'   When NULL, paths are counted.
#' @param orders Optional level order vector applied to categorical axes.
#' @param name Optional label appended to the output file name.
#' @param height Height of the output PDF.
#' @param width Width of the output PDF.
#' @param plotCols Number of columns for multi-panel layout.
#' @param plotRows Number of rows for multi-panel layout.
#' @param size_by Whether ribbon size should reflect summed abundance or sample
#'   count per path in \code{mode = "feature"}. One of \code{"abundance"} or
#'   \code{"count"}.
#' @param alpha Ribbon transparency.
#' @param mycol Optional color vector. In \code{mode = "feature"} it is mapped
#'   to features. In \code{mode = "transition"} it is mapped to fill groups.
#'
#' @return Invisibly returns a list of plot objects. A PDF is written to the
#'   dated project \code{pdf/} directory.
#'
#' @param patchwork Logical. If \code{TRUE}, skip saving and return the plot object(s) for use with \code{Gg_patchwork()} or the \pkg{patchwork} package. Default \code{FALSE}.
#' @export
Go_alluvialplot <- function(project,
                            SigASVs = NULL,
                            map = NULL,
                            targets.bac = NULL,
                            target.rank = NULL,
                            outcome = NULL,
                            column1 = NULL,
                            column2 = NULL,
                            Addcolumn = NULL,
                            data = NULL,
                            mode = c("feature", "transition"),
                            axes = NULL,
                            fill_var = NULL,
                            facet_var = NULL,
                            weight_var = NULL,
                            orders = NULL,
                            name = NULL,
                            height = 2,
                            width = 3,
                            plotCols = 2,
                            plotRows = 1,
                            size_by = c("abundance", "count"),
                            alpha = 0.72,
                            mycol = NULL,
                            patchwork = FALSE) {

  mode <- match.arg(mode)
  size_by <- match.arg(size_by)

  if (identical(mode, "feature")) {
    if (!requireNamespace("reshape2", quietly = TRUE)) {
      stop("reshape2 is required for Go_alluvialplot(mode = 'feature').")
    }
    if (!requireNamespace("alluvial", quietly = TRUE) &&
        !requireNamespace("ggalluvial", quietly = TRUE)) {
      stop("Feature mode requires either the 'alluvial' or 'ggalluvial' package.")
    }
  }
  if (identical(mode, "transition")) {
    if (!requireNamespace("ggalluvial", quietly = TRUE) &&
        !requireNamespace("alluvial", quietly = TRUE)) {
      stop("Transition mode requires either the 'alluvial' or 'ggalluvial' package.")
    }
    if (!requireNamespace("dplyr", quietly = TRUE)) {
      stop("dplyr is required for Go_alluvialplot(mode = 'transition').")
    }
  }

  if (!is.null(dev.list())) {
    grDevices::dev.off()
  }

  read_input <- function(x, rownames_col = FALSE) {
    if (is.character(x) && length(x) == 1) {
      if (isTRUE(rownames_col)) {
        return(utils::read.csv(x, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE))
      }
      return(utils::read.csv(x, check.names = FALSE, stringsAsFactors = FALSE))
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

  clean_tag <- function(x) {
    x <- gsub("[^A-Za-z0-9._-]+", "_", x)
    x <- gsub("_+", "_", x)
    x <- gsub("^_|_$", "", x)
    x
  }

  default_cols <- c(
    "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",
    "#A65628", "#F781BF", "#1B9E77", "#D95F02", "#7570B3",
    "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666"
  )

  build_plot_title <- function(target, rank_name = NULL, width = 28) {
    title_txt <- as.character(target)
    paste(strwrap(title_txt, width = width), collapse = "\n")
  }

  resolve_palette <- function(features, user_cols = NULL) {
    features <- unique(as.character(features))
    if (length(features) == 0) {
      return(stats::setNames(character(0), character(0)))
    }

    if (!is.null(user_cols)) {
      if (!is.null(names(user_cols)) && all(features %in% names(user_cols))) {
        return(user_cols[features])
      }
      vals <- rep(user_cols, length.out = length(features))
      return(stats::setNames(vals, features))
    }

    if (length(default_cols) > 0) {
      vals <- rep(default_cols, length.out = length(features))
    } else if ("hcl.colors" %in% getNamespaceExports("grDevices")) {
      vals <- grDevices::hcl.colors(length(features), palette = "Dark 3")
    } else if (requireNamespace("scales", quietly = TRUE)) {
      vals <- scales::hue_pal(l = 55, c = 65)(length(features))
    } else {
      vals <- grDevices::rainbow(length(features), s = 0.65, v = 0.75)
    }
    stats::setNames(vals, features)
  }

  render_plots <- function(plots, pdf_path, width, height, plotCols, plotRows) {
    grDevices::pdf(pdf_path, height = height, width = width)
    on.exit(grDevices::dev.off(), add = TRUE)

    if (length(plots) == 1) {
      print(plots[[1]])
    } else if (exists("multiplot", mode = "function")) {
      multiplot(plotlist = plots, cols = plotCols, rows = plotRows)
    } else {
      for (p in plots) print(p)
    }
  }

  if (is.null(project) || !nzchar(project)) {
    stop("`project` must be provided.")
  }

  output_dirs <- Go_path(project = project, pdf = "yes", table = "no", path = NULL)
  out_path <- output_dirs$pdf

  if (identical(mode, "feature")) {
    tab <- read_input(SigASVs, rownames_col = TRUE)
    sampledata <- read_input(map, rownames_col = TRUE)

    if (is.null(targets.bac) || length(targets.bac) == 0) {
      stop("`targets.bac` must contain at least one target taxon in feature mode.")
    }
    if (is.null(target.rank) || !target.rank %in% colnames(tab)) {
      stop("`target.rank` must be a column present in `SigASVs`.")
    }

    axis_vars <- c(column1, column2, Addcolumn, outcome)
    axis_vars <- axis_vars[!is.na(axis_vars) & nzchar(axis_vars)]
    axis_vars <- unique(axis_vars)
    if (!length(axis_vars)) {
      stop("At least one metadata axis is required in feature mode.")
    }

    missing_axis_cols <- setdiff(axis_vars, colnames(sampledata))
    if (length(missing_axis_cols) > 0) {
      stop("Metadata columns not found in `map`: ", paste(missing_axis_cols, collapse = ", "))
    }

    if ((is.null(rownames(sampledata)) || !length(rownames(sampledata)) || all(!nzchar(rownames(sampledata)))) &&
        "SampleID" %in% colnames(sampledata)) {
      rownames(sampledata) <- sampledata$SampleID
    }

    numeric_cols <- names(tab)[vapply(tab, is.numeric, logical(1))]
    sample_ids <- unique(c(rownames(sampledata), sampledata$SampleID))
    sample_ids <- sample_ids[!is.na(sample_ids) & nzchar(sample_ids)]
    sample_cols <- intersect(numeric_cols, sample_ids)
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

    plots <- list()
    rendered_targets <- character(0)
    use_classic_alluvial <- requireNamespace("alluvial", quietly = TRUE)

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

      feature_counts <- table(taxa_melt$feature)
      taxa_melt$feature_count <- as.integer(feature_counts[match(taxa_melt$feature, names(feature_counts))])
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
      key_df <- taxa_melt[, group_cols, drop = FALSE]
      key_id <- do.call(paste, c(lapply(key_df, as.character), sep = "\r"))

      count_tab <- stats::aggregate(
        x = list(sample_n = rep(1, length(key_id))),
        by = c(key_df, list(.key = key_id)),
        FUN = sum
      )
      sum_tab <- stats::aggregate(
        x = list(abundance_sum = as.numeric(taxa_melt$value)),
        by = c(key_df, list(.key = key_id)),
        FUN = sum,
        na.rm = TRUE
      )

      df_plot <- merge(count_tab, sum_tab[, c(".key", "abundance_sum"), drop = FALSE], by = ".key", all = TRUE, sort = FALSE)
      df_plot <- df_plot[, c(group_cols, "sample_n", "abundance_sum"), drop = FALSE]
      df_plot <- as.data.frame(df_plot, stringsAsFactors = FALSE)
      df_plot[["sample_n"]] <- as.numeric(unlist(df_plot[["sample_n"]], use.names = FALSE))
      df_plot[["abundance_sum"]] <- as.numeric(unlist(df_plot[["abundance_sum"]], use.names = FALSE))
      if (identical(size_by, "count")) {
        df_plot[["weight"]] <- base::as.numeric(df_plot[["sample_n"]])
      } else {
        df_plot[["weight"]] <- base::log1p(base::as.numeric(df_plot[["abundance_sum"]]))
      }
      df_plot <- df_plot[df_plot[["weight"]] > 0 & is.finite(df_plot[["weight"]]), , drop = FALSE]
      df_plot[["plot_weight"]] <- sqrt(df_plot[["weight"]])
      if (nrow(df_plot) == 0) {
        message("[Go_alluvialplot] No aggregated paths remained for target: ", target)
        next
      }

      feature_levels <- unique(df_plot$feature)
      feature_palette <- resolve_palette(feature_levels, mycol)
      plot_cols <- grDevices::adjustcolor(
        feature_palette[match(df_plot$feature, names(feature_palette))],
        alpha.f = alpha
      )
      axis_names <- c(target.rank, axis_vars)
      flow_df <- df_plot[, c("label", axis_vars), drop = FALSE]

      if (use_classic_alluvial) {
        plots[[length(plots) + 1]] <- list(
          engine = "alluvial",
          data = flow_df,
          freq = df_plot$plot_weight,
          col = plot_cols,
          axis_labels = axis_names,
          title = build_plot_title(target, target.rank)
        )
      } else {
        axis_positions <- seq_along(c("label", axis_vars))
        names(axis_positions) <- c("label", axis_vars)
        df_plot$alluvium <- seq_len(nrow(df_plot))
        lodes_plot <- ggalluvial::to_lodes_form(
          data = df_plot,
          axes = c("label", axis_vars),
          key = "x_axis",
          value = "stratum",
          id = "alluvium"
        )
        lodes_plot$x_num <- unname(axis_positions[as.character(lodes_plot$x_axis)])

        p <- ggplot2::ggplot(
          lodes_plot,
          ggplot2::aes(
            x = x_num,
            stratum = stratum,
            alluvium = alluvium,
            y = plot_weight,
            fill = feature
          )
        ) +
          ggalluvial::geom_alluvium(
            ggplot2::aes(fill = feature),
            alpha = alpha,
            width = 0.30,
            knot.pos = 0.35,
            color = NA
          ) +
          ggalluvial::geom_stratum(
            width = 0.30, fill = "grey96",
            color = "grey55", linewidth = 0.55
          ) +
          ggalluvial::stat_stratum(
            ggplot2::aes(label = ggplot2::after_stat(stratum)),
            geom = "text",
            size = 3.4,
            fontface = "bold",
            color = "grey15"
          ) +
          ggplot2::scale_fill_manual(values = feature_palette, drop = FALSE, name = target.rank) +
          ggplot2::scale_x_continuous(
            breaks = axis_positions,
            labels = axis_names,
            expand = ggplot2::expansion(mult = c(0.06, 0.10))
          ) +
          ggplot2::labs(
            title = build_plot_title(target, target.rank),
            x = NULL,
            y = if (identical(size_by, "count")) "Sample count" else "log1p(Abundance)"
          ) +
          ggplot2::theme_minimal(base_size = 11) +
          ggplot2::theme(
            legend.position = "none",
            axis.title = ggplot2::element_blank(),
            axis.text.y = ggplot2::element_blank(),
            axis.ticks = ggplot2::element_blank(),
            axis.text.x = ggplot2::element_text(size = 10, colour = "grey20", face = "bold"),
            plot.title = ggplot2::element_text(hjust = 0.5, size = 11, face = "bold"),
            panel.grid = ggplot2::element_blank(),
            plot.background = ggplot2::element_rect(fill = "white", color = NA),
            panel.background = ggplot2::element_rect(fill = "white", color = NA)
          )
        plots[[length(plots) + 1]] <- p
      }

      rendered_targets <- c(rendered_targets, target)
    }

    if (length(plots) == 0) {
      stop("No alluvial plots were generated. Check `targets.bac`, `target.rank`, and sample overlap.")
    }

    file_stub <- sprintf(
      "Alluvial.%s.feature.%s%s.pdf",
      project,
      ifelse(is.null(name), "", paste0(clean_tag(name), ".")),
      format(Sys.Date(), "%y%m%d")
    )
    if (isTRUE(patchwork)) return(invisible(plots))
    pdf_path <- file.path(out_path, file_stub)
    grDevices::pdf(pdf_path, height = height, width = width)
    on.exit(grDevices::dev.off(), add = TRUE)
    for (p in plots) {
      if (is.list(p) && identical(p$engine, "alluvial")) {
        old_par <- graphics::par(no.readonly = TRUE)
        on.exit(graphics::par(old_par), add = TRUE)
        graphics::par(
          mar = c(4.2, 1.2, 1.6, 1.2),
          oma = c(0, 0, 2.2, 0),
          xpd = NA,
          bg = "white"
        )
        alluvial::alluvial(
          p$data,
          freq = p$freq,
          border = NA,
          alpha = alpha,
          col = p$col,
          cex = 0.75,
          axis_labels = p$axis_labels
        )
        graphics::mtext(
          text = p$title,
          side = 3,
          outer = TRUE,
          line = 0.3,
          cex = 0.95,
          font = 2
        )
      } else {
        print(p)
      }
    }
    return(invisible(plots))
  }

  ## transition mode ---------------------------------------------------------
  trans_df <- read_input(data, rownames_col = FALSE)
  if (is.null(axes) || length(axes) < 2) {
    stop("`axes` must contain at least two columns in transition mode.")
  }
  missing_axes <- setdiff(axes, names(trans_df))
  if (length(missing_axes) > 0) {
    stop("Transition axes not found in `data`: ", paste(missing_axes, collapse = ", "))
  }
  if (!is.null(fill_var) && !fill_var %in% names(trans_df)) {
    stop("`fill_var` not found in `data`.")
  }
  if (!is.null(facet_var) && !facet_var %in% names(trans_df)) {
    stop("`facet_var` not found in `data`.")
  }
  if (!is.null(weight_var) && !weight_var %in% names(trans_df)) {
    stop("`weight_var` not found in `data`.")
  }

  trans_df <- as.data.frame(trans_df, stringsAsFactors = FALSE)
  for (ax in axes) {
    trans_df[[ax]] <- as.character(trans_df[[ax]])
    trans_df[[ax]][is.na(trans_df[[ax]]) | !nzchar(trans_df[[ax]])] <- "NA"
    if (!is.null(orders)) {
      observed_levels <- unique(trans_df[[ax]])
      levs <- c(intersect(orders, observed_levels), setdiff(observed_levels, orders))
      trans_df[[ax]] <- factor(trans_df[[ax]], levels = levs)
    } else {
      trans_df[[ax]] <- factor(trans_df[[ax]])
    }
  }

  if (is.null(fill_var)) {
    trans_df$.fill_group <- "All"
    fill_var <- ".fill_group"
  }
  trans_df[[fill_var]] <- as.character(trans_df[[fill_var]])
  trans_df[[fill_var]][is.na(trans_df[[fill_var]]) | !nzchar(trans_df[[fill_var]])] <- "NA"

  group_cols <- c(fill_var, facet_var, axes)
  group_cols <- group_cols[!is.na(group_cols) & nzchar(group_cols)]
  df_plot <- dplyr::summarise(
    dplyr::group_by(trans_df, dplyr::across(dplyr::all_of(group_cols))),
    weight = if (is.null(weight_var)) dplyr::n() else sum(rlang::.data[[weight_var]], na.rm = TRUE),
    .groups = "drop"
  )
  df_plot <- as.data.frame(df_plot, stringsAsFactors = FALSE)
  df_plot[["weight"]] <- as.numeric(df_plot[["weight"]])
  df_plot <- df_plot[is.finite(df_plot[["weight"]]) & df_plot[["weight"]] > 0, , drop = FALSE]
  if (nrow(df_plot) == 0) {
    stop("No transition paths remained after aggregation.")
  }

  file_stub <- sprintf(
    "Alluvial.%s.transition.%s%s.pdf",
    project,
    ifelse(is.null(name), "", paste0(clean_tag(name), ".")),
    format(Sys.Date(), "%y%m%d")
  )
  pdf_path <- file.path(out_path, file_stub)
  fill_levels <- unique(df_plot[[fill_var]])
  fill_levels <- fill_levels[!is.na(fill_levels) & nzchar(as.character(fill_levels)) & as.character(fill_levels) != "NA"]
  fill_palette <- resolve_palette(fill_levels, mycol)
  use_classic_alluvial <- FALSE

  if (use_classic_alluvial) {
    df_plot[["plot_weight"]] <- sqrt(df_plot[["weight"]])
    split_list <- if (is.null(facet_var)) {
      list(.all = df_plot)
    } else {
      split(df_plot, df_plot[[facet_var]], drop = TRUE)
    }

    grDevices::pdf(pdf_path, height = height, width = width)
    on.exit(grDevices::dev.off(), add = TRUE)
    for (nm in names(split_list)) {
      sub_df <- split_list[[nm]]
      sub_df <- sub_df[is.finite(sub_df[["plot_weight"]]) & sub_df[["plot_weight"]] > 0, , drop = FALSE]
      if (!nrow(sub_df)) next

      old_par <- graphics::par(no.readonly = TRUE)
      on.exit(graphics::par(old_par), add = TRUE)
      graphics::par(
        mar = c(4.2, 1.2, 3.2, 8.8),
        oma = c(0, 0, 2.2, 0),
        xpd = NA,
        bg = "white"
      )

      flow_cols <- grDevices::adjustcolor(
        fill_palette[match(sub_df[[fill_var]], names(fill_palette))],
        alpha.f = alpha
      )
      plot_title <- if (is.null(name)) "Transition alluvial plot" else name
      if (!is.null(facet_var)) {
        plot_title <- paste0(plot_title, "\n", nm)
      }

      alluvial::alluvial(
        sub_df[, axes, drop = FALSE],
        freq = sub_df[["plot_weight"]],
        border = NA,
        alpha = alpha,
        col = flow_cols,
        cex = 0.85,
        axis_labels = rep("", length(axes))
      )
      graphics::mtext(
        text = axes,
        side = 3,
        at = seq_along(axes),
        line = -2.1,
        cex = 0.95,
        font = 2
      )
      graphics::legend(
        "right",
        legend = names(fill_palette),
        fill = fill_palette,
        border = NA,
        bty = "n",
        cex = 0.72,
        title = NULL,
        xpd = NA,
        inset = c(-0.18, 0)
      )
      graphics::mtext(
        text = build_plot_title(plot_title, width = 36),
        side = 3,
        outer = TRUE,
        line = 0.3,
        cex = 0.95,
        font = 2
      )
    }
    return(invisible(split_list))
  }

  axis_mapping <- list()
  for (i in seq_along(axes)) {
    axis_mapping[[paste0("axis", i)]] <- as.name(axes[i])
  }
  axis_mapping$y <- quote(weight)
  axis_mapping$fill <- as.name(fill_var)

  p <- ggplot2::ggplot(
    df_plot,
    do.call(ggplot2::aes, axis_mapping)
  ) +
    ggalluvial::geom_alluvium(
      ggplot2::aes(fill = !!as.name(fill_var)),
      width = 0.32, alpha = alpha, knot.pos = 0.4
    ) +
    ggalluvial::geom_stratum(
      width = 0.32, fill = "grey96",
      color = "grey55", linewidth = 0.55
    ) +
    ggalluvial::stat_stratum(
      ggplot2::aes(label = ggplot2::after_stat(stratum)),
      geom = "text",
      size = 3.5,
      fontface = "bold",
      color = "grey15"
    ) +
    ggplot2::scale_fill_manual(
      values = fill_palette,
      name = if (identical(fill_var, "trajectory")) "Trajectory" else gsub("_", " ", tools::toTitleCase(fill_var)),
      na.translate = FALSE
    ) +
    ggplot2::scale_x_discrete(
      limits = axes,
      labels = axes,
      expand = c(0.06, 0.06)
    ) +
    ggplot2::labs(
      title = ifelse(is.null(name), "Transition alluvial plot", name),
      x = NULL,
      y = ifelse(is.null(weight_var), "Number of rows", weight_var)
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5, size = 13),
      axis.text.x = ggplot2::element_text(size = 11, face = "bold", colour = "grey20"),
      axis.text.y = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank(),
      legend.position = "right",
      plot.background = ggplot2::element_rect(fill = "white", color = NA),
      panel.background = ggplot2::element_rect(fill = "white", color = NA)
    )

  if (!is.null(facet_var)) {
    p <- p + ggplot2::facet_wrap(stats::as.formula(paste("~", facet_var)))
  }

  if (isTRUE(patchwork)) return(invisible(list(p)))
  render_plots(
    list(p), pdf_path = pdf_path,
    width = width, height = height, plotCols = 1, plotRows = 1
  )
  invisible(list(p))
}
