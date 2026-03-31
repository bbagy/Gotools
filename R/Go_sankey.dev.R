#' Generate sankey plots for feature or transition data
#'
#' Uses \code{ggsankey::geom_sankey()} to draw sankey plots from either
#' feature tables plus metadata or from wide-format transition tables.
#' \itemize{
#'   \item \code{mode = "feature"} keeps the original SigASVs + metadata
#'   workflow for taxa-centered sankey plots.
#'   \item \code{mode = "transition"} draws transition sankey plots directly
#'   from a wide-format table such as \code{V1 -> V2 -> V3} states.
#' }
#'
#' @param project Project name used for output directory and file naming.
#' @param SigASVs Either a CSV path or a data frame containing significant
#'   features and taxonomy columns. Used in \code{mode = "feature"}.
#' @param map Either a CSV path or a data frame containing sample metadata.
#'   Used in \code{mode = "feature"}.
#' @param targets.bac Character vector of taxa to plot. Used in
#'   \code{mode = "feature"}.
#' @param target.rank Taxonomy column used to match \code{targets.bac}.
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
#' @param fill_var Optional column used for flow fill in
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
#' @param size_by Whether ribbon size should reflect summed abundance or
#'   sample count per path in \code{mode = "feature"}. One of
#'   \code{"abundance"} or \code{"count"}.
#' @param alpha Flow transparency.
#' @param mycol Optional color vector. In \code{mode = "feature"} it is
#'   mapped to features. In \code{mode = "transition"} it is mapped to
#'   fill groups.
#'
#' @return Invisibly returns a list of ggplot objects. A PDF is written to
#'   the dated project \code{pdf/} directory.
#'
#' @export
Go_sankey <- function(project,
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
                      height = 4,
                      width = 6,
                      plotCols = 2,
                      plotRows = 1,
                      size_by = c("abundance", "count"),
                      alpha = 0.72,
                      mycol = NULL) {

  mode <- match.arg(mode)
  size_by <- match.arg(size_by)

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 is required for Go_sankey().")
  }
  if (!requireNamespace("ggsankey", quietly = TRUE)) {
    stop(paste0(
      "ggsankey is required for Go_sankey(). ",
      "Install with: remotes::install_github('davidsjoberg/ggsankey')"
    ))
  }
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("dplyr is required for Go_sankey().")
  }
  if (!requireNamespace("rlang", quietly = TRUE)) {
    stop("rlang is required for Go_sankey().")
  }
  if (identical(mode, "feature")) {
    if (!requireNamespace("reshape2", quietly = TRUE)) {
      stop("reshape2 is required for Go_sankey(mode = 'feature').")
    }
  }

  if (!is.null(dev.list())) {
    grDevices::dev.off()
  }

  # --- helpers ---------------------------------------------------------------

  read_input <- function(x, rownames_col = FALSE) {
    if (is.character(x) && length(x) == 1) {
      if (isTRUE(rownames_col)) {
        return(utils::read.csv(
          x, row.names = 1,
          check.names = FALSE, stringsAsFactors = FALSE
        ))
      }
      return(utils::read.csv(
        x, check.names = FALSE, stringsAsFactors = FALSE
      ))
    }
    if (is.data.frame(x)) return(x)
    stop("Inputs must be either a CSV file path or a data frame.")
  }

  sanitize_feature_label <- function(x) {
    x <- as.character(x)
    x[is.na(x) | !nzchar(trimws(x))] <- "Unknown"
    x <- gsub("\\s+", "_", x)
    x
  }

  compact_feature_labels <- function(x) {
    x_chr <- as.character(x)
    looks_like_asv <- grepl("^[ACGTNacgtn-]+$", x_chr) & nchar(x_chr) >= 20
    out <- x_chr
    if (any(looks_like_asv, na.rm = TRUE)) {
      out[looks_like_asv] <- paste0("ASV_", seq_len(sum(looks_like_asv)))
    }
    out
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
    paste(strwrap(as.character(target), width = width), collapse = "\n")
  }

  resolve_palette <- function(features, user_cols = NULL) {
    features <- unique(as.character(features))
    if (length(features) == 0) {
      return(stats::setNames(character(0), character(0)))
    }
    if (!is.null(user_cols)) {
      if (!is.null(names(user_cols)) &&
          all(features %in% names(user_cols))) {
        return(user_cols[features])
      }
      return(stats::setNames(
        rep(user_cols, length.out = length(features)), features
      ))
    }
    stats::setNames(
      rep(default_cols, length.out = length(features)), features
    )
  }

  render_plots <- function(plots, pdf_path, width, height) {
    grDevices::pdf(pdf_path, width = width, height = height)
    on.exit(grDevices::dev.off(), add = TRUE)
    for (p in plots) print(p)
  }

  build_sankey_plot <- function(
      flow_df, axis_cols, fill_col, plot_title, palette_vals
  ) {
    sankey_long <- rlang::exec(
      ggsankey::make_long,
      .df = flow_df,
      !!!rlang::syms(axis_cols),
      value = rlang::sym("plot_weight")
    )
    row_idx <- rep(seq_len(nrow(flow_df)), each = length(axis_cols))
    sankey_long[[".fill_group"]] <- flow_df[[fill_col]][row_idx]

    node_df <- stats::aggregate(
      value ~ x + node,
      data = sankey_long[!is.na(sankey_long$node), c("x", "node", "value")],
      FUN = sum
    )
    node_key <- unique(sankey_long[!is.na(sankey_long$node), c("x", "node")])
    node_key[[".ord"]] <- seq_len(nrow(node_key))
    node_df <- merge(node_df, node_key, by = c("x", "node"), all.x = TRUE, sort = FALSE)
    node_df <- node_df[order(node_df$x, node_df$.ord), , drop = FALSE]
    node_df <- do.call(rbind, lapply(split(node_df, node_df$x), function(d) {
      d$ymax <- cumsum(d$value)
      d$ymin <- d$ymax - d$value
      d$y <- (d$ymin + d$ymax) / 2
      if (nrow(d) > 1) {
        total_h <- max(d$ymax, na.rm = TRUE)
        min_gap <- max(total_h * 0.08, 0.12)
        ord_y <- order(d$y)
        y_adj <- d$y[ord_y]
        for (i in 2:length(y_adj)) {
          if ((y_adj[i] - y_adj[i - 1]) < min_gap) {
            y_adj[i] <- y_adj[i - 1] + min_gap
          }
        }
        overflow <- max(y_adj, na.rm = TRUE) - total_h
        if (is.finite(overflow) && overflow > 0) {
          y_adj <- y_adj - overflow
        }
        underflow <- min(y_adj, na.rm = TRUE)
        if (is.finite(underflow) && underflow < 0) {
          y_adj <- y_adj - underflow
        }
        d$y[ord_y] <- y_adj
      }
      d
    }))

    ggplot2::ggplot(
      sankey_long,
      ggplot2::aes(
        x = .data$x,
        next_x = .data$next_x,
        node = .data$node,
        next_node = .data$next_node,
        fill = .data$.fill_group,
        value = .data$value
      )
    ) +
      ggsankey::geom_sankey(
        flow.alpha = alpha,
        node.color = "grey35",
        width = 0.08,
        space = 0.02
      ) +
      ggplot2::geom_label(
        data = node_df,
        ggplot2::aes(x = .data$x, y = .data$y, label = .data$node),
        inherit.aes = FALSE,
        size = 3.0,
        fontface = "bold",
        color = "grey15",
        fill = "white",
        alpha = 0.9,
        label.size = 0.2
      ) +
      ggplot2::scale_fill_manual(values = palette_vals, drop = FALSE) +
      ggplot2::labs(title = plot_title, x = NULL, y = NULL) +
      ggsankey::theme_sankey(base_size = 11) +
      ggplot2::theme(
        legend.position = "none",
        plot.title = ggplot2::element_text(
          face = "bold", hjust = 0.5, size = 11
        ),
        axis.text.x = ggplot2::element_text(
          size = 10, face = "bold", colour = "grey20"
        ),
        axis.text.y = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank(),
        panel.grid = ggplot2::element_blank(),
        plot.background = ggplot2::element_rect(fill = "white", color = NA),
        panel.background = ggplot2::element_rect(fill = "white", color = NA)
      )
  }

  # --- validation ------------------------------------------------------------

  if (is.null(project) || !nzchar(project)) {
    stop("`project` must be provided.")
  }

  out_path <- Go_path(
    project = project, pdf = "yes", table = "no", path = NULL
  )$pdf

  # ===========================================================================
  # feature mode
  # ===========================================================================
  if (identical(mode, "feature")) {
    tab <- read_input(SigASVs, rownames_col = TRUE)
    sampledata <- read_input(map, rownames_col = TRUE)

    if (is.null(targets.bac) || length(targets.bac) == 0) {
      stop("`targets.bac` must contain at least one target taxon.")
    }
    if (is.null(target.rank) || !target.rank %in% colnames(tab)) {
      stop("`target.rank` must be a column present in `SigASVs`.")
    }

    axis_vars <- c(column1, column2, Addcolumn, outcome)
    axis_vars <- unique(axis_vars[!is.na(axis_vars) & nzchar(axis_vars)])
    if (!length(axis_vars)) {
      stop("At least one metadata axis is required in feature mode.")
    }

    missing_axis_cols <- setdiff(axis_vars, colnames(sampledata))
    if (length(missing_axis_cols) > 0) {
      stop(
        "Metadata columns not found in `map`: ",
        paste(missing_axis_cols, collapse = ", ")
      )
    }

    if ((is.null(rownames(sampledata)) ||
         !length(rownames(sampledata)) ||
         all(!nzchar(rownames(sampledata)))) &&
        "SampleID" %in% colnames(sampledata)) {
      rownames(sampledata) <- sampledata$SampleID
    }

    numeric_cols <- names(tab)[vapply(tab, is.numeric, logical(1))]
    sample_ids <- unique(c(rownames(sampledata), sampledata$SampleID))
    sample_ids <- sample_ids[!is.na(sample_ids) & nzchar(sample_ids)]
    sample_cols <- intersect(numeric_cols, sample_ids)
    if (length(sample_cols) == 0) {
      stop("No overlapping sample columns found between `SigASVs` and `map`.")
    }

    tab <- tab[, c(setdiff(colnames(tab), sample_cols), sample_cols),
               drop = FALSE]
    sampledata <- sampledata[sample_cols, , drop = FALSE]
    sampledata$SampleID <- rownames(sampledata)

    tab$RowSum <- rowSums(tab[, sample_cols, drop = FALSE], na.rm = TRUE)
    raw_display <- if (!is.null(target.rank) && target.rank %in% colnames(tab)) {
      tab[[target.rank]]
    } else {
      rownames(tab)
    }
    raw_display <- sanitize_feature_label(raw_display)
    if (identical(target.rank, "ASV")) {
      raw_display <- compact_feature_labels(raw_display)
    }
    tab$.display_label <- raw_display
    tab$.feature_label <- paste0(raw_display, "_", tab$RowSum)
    rownames(tab) <- make.unique(tab$.feature_label)

    plots <- list()

    for (target in targets.bac) {
      tab.sel <- tab[
        as.character(tab[[target.rank]]) == target, , drop = FALSE
      ]
      if (nrow(tab.sel) == 0) {
        message("[Go_sankey] No rows found for target: ", target)
        next
      }

      feature_labels <- rownames(tab.sel)
      display_map <- stats::setNames(tab.sel$.display_label, feature_labels)
      taxa_tab <- data.frame(
        SampleID = sample_cols,
        t(as.matrix(tab.sel[, sample_cols, drop = FALSE])),
        check.names = FALSE,
        stringsAsFactors = FALSE
      )
      colnames(taxa_tab)[-1] <- feature_labels
      taxa_tab <- merge(
        sampledata, taxa_tab,
        by = "SampleID", all.x = FALSE, sort = FALSE
      )

      taxa_melt <- reshape2::melt(
        taxa_tab,
        id.vars = c("SampleID", axis_vars),
        measure.vars = feature_labels,
        variable.name = "feature",
        value.name = "value"
      )
      taxa_melt$value <- suppressWarnings(as.numeric(taxa_melt$value))
      taxa_melt <- taxa_melt[
        is.finite(taxa_melt$value) & taxa_melt$value > 0, , drop = FALSE
      ]
      if (nrow(taxa_melt) == 0) {
        message(
          "[Go_sankey] No positive values remained for target: ", target
        )
        next
      }

      feature_counts <- table(taxa_melt$feature)
      taxa_melt$feature_count <- as.integer(
        feature_counts[match(taxa_melt$feature, names(feature_counts))]
      )
      taxa_melt$label <- unname(display_map[as.character(taxa_melt$feature)])
      taxa_melt$label[is.na(taxa_melt$label) | !nzchar(taxa_melt$label)] <- "Unknown"

      for (var_name in axis_vars) {
        taxa_melt[[var_name]] <- as.character(taxa_melt[[var_name]])
        taxa_melt[[var_name]][
          is.na(taxa_melt[[var_name]]) | !nzchar(taxa_melt[[var_name]])
        ] <- "NA"
        if (!is.null(orders)) {
          obs <- unique(taxa_melt[[var_name]])
          var_levels <- c(intersect(orders, obs), setdiff(obs, orders))
          taxa_melt[[var_name]] <- factor(taxa_melt[[var_name]], var_levels)
        } else {
          taxa_melt[[var_name]] <- factor(taxa_melt[[var_name]])
        }
      }

      group_cols <- c("label", "feature", axis_vars)
      df_plot <- dplyr::summarise(
        dplyr::group_by(
          taxa_melt,
          dplyr::across(dplyr::all_of(group_cols))
        ),
        sample_n = dplyr::n(),
        abundance_sum = sum(.data$value, na.rm = TRUE),
        .groups = "drop"
      )
      df_plot <- as.data.frame(df_plot, stringsAsFactors = FALSE)
      sample_n_vec <- suppressWarnings(as.numeric(unlist(df_plot[, "sample_n", drop = TRUE], use.names = FALSE)))
      abundance_vec <- suppressWarnings(as.numeric(unlist(df_plot[, "abundance_sum", drop = TRUE], use.names = FALSE)))
      if (identical(size_by, "count")) {
        weight_vec <- sample_n_vec
      } else {
        weight_vec <- base::log1p(abundance_vec)
      }
      df_plot <- data.frame(
        df_plot[, group_cols, drop = FALSE],
        sample_n = sample_n_vec,
        abundance_sum = abundance_vec,
        weight = weight_vec,
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
      df_plot <- df_plot[
        is.finite(df_plot[["weight"]]) & df_plot[["weight"]] > 0, ,
        drop = FALSE
      ]
      df_plot[["plot_weight"]] <- sqrt(df_plot[["weight"]])
      if (nrow(df_plot) == 0) {
        message(
          "[Go_sankey] No aggregated paths remained for target: ", target
        )
        next
      }

      feature_palette <- resolve_palette(unique(df_plot$feature), mycol)
      plots[[length(plots) + 1]] <- build_sankey_plot(
        flow_df = df_plot,
        axis_cols = c("label", axis_vars),
        fill_col = "feature",
        plot_title = build_plot_title(target, target.rank),
        palette_vals = feature_palette
      )
    }

    if (length(plots) == 0) {
      stop("No sankey plots were generated in feature mode.")
    }

    file_stub <- sprintf(
      "Sankey.%s.feature.%s%s.pdf",
      project,
      ifelse(is.null(name), "", paste0(clean_tag(name), ".")),
      format(Sys.Date(), "%y%m%d")
    )
    render_plots(plots, file.path(out_path, file_stub), width, height)
    return(invisible(plots))
  }

  # ===========================================================================
  # transition mode
  # ===========================================================================
  trans_df <- read_input(data, rownames_col = FALSE)
  if (is.null(axes) || length(axes) < 2) {
    stop("`axes` must contain at least two columns in transition mode.")
  }
  missing_axes <- setdiff(axes, names(trans_df))
  if (length(missing_axes) > 0) {
    stop(
      "Transition axes not found in `data`: ",
      paste(missing_axes, collapse = ", ")
    )
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
    trans_df[[ax]][
      is.na(trans_df[[ax]]) | !nzchar(trans_df[[ax]])
    ] <- "NA"
    if (!is.null(orders)) {
      obs <- unique(trans_df[[ax]])
      levs <- c(intersect(orders, obs), setdiff(obs, orders))
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
  trans_df[[fill_var]][
    is.na(trans_df[[fill_var]]) | !nzchar(trans_df[[fill_var]])
  ] <- "NA"

  group_cols <- c(fill_var, facet_var, axes)
  group_cols <- group_cols[!is.na(group_cols) & nzchar(group_cols)]
  df_plot <- dplyr::summarise(
    dplyr::group_by(trans_df, dplyr::across(dplyr::all_of(group_cols))),
    weight = if (is.null(weight_var)) {
      dplyr::n()
    } else {
      sum(.data[[weight_var]], na.rm = TRUE)
    },
    .groups = "drop"
  )
  df_plot <- as.data.frame(df_plot, stringsAsFactors = FALSE)
  weight_vec <- suppressWarnings(as.numeric(unlist(df_plot[, "weight", drop = TRUE], use.names = FALSE)))
  df_plot <- data.frame(
    df_plot[, setdiff(colnames(df_plot), "weight"), drop = FALSE],
    weight = weight_vec,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  df_plot <- df_plot[
    is.finite(df_plot[["weight"]]) & df_plot[["weight"]] > 0, ,
    drop = FALSE
  ]
  if (nrow(df_plot) == 0) {
    stop("No transition paths remained after aggregation.")
  }
  df_plot[["plot_weight"]] <- sqrt(df_plot[["weight"]])

  fill_palette <- resolve_palette(unique(df_plot[[fill_var]]), mycol)
  split_list <- if (is.null(facet_var)) {
    list(.all = df_plot)
  } else {
    split(df_plot, df_plot[[facet_var]], drop = TRUE)
  }

  plots <- lapply(names(split_list), function(nm) {
    sub_df <- split_list[[nm]]
    p_title <- if (!is.null(name)) {
      if (is.null(facet_var)) name else paste(name, nm, sep = "\n")
    } else {
      if (is.null(facet_var)) "Transition sankey plot" else nm
    }
    build_sankey_plot(
      flow_df = sub_df,
      axis_cols = axes,
      fill_col = fill_var,
      plot_title = p_title,
      palette_vals = fill_palette
    )
  })

  file_stub <- sprintf(
    "Sankey.%s.transition.%s%s.pdf",
    project,
    ifelse(is.null(name), "", paste0(clean_tag(name), ".")),
    format(Sys.Date(), "%y%m%d")
  )
  render_plots(plots, file.path(out_path, file_stub), width, height)
  invisible(plots)
}
