#' Generate Rarefaction Curves for Microbiome Data
#'
#' This function creates rarefaction curves for microbiome data using a Phyloseq object.
#' It provides an assessment of species richness across varying sequencing depths.
#'
#' @param physeq_object A Phyloseq object containing OTU/ASV counts and associated sample data.
#' @param step Step size for rarefaction curve plotting.
#' @param label Optional sample label mapping. Supply either a metadata column
#'   name or a vector with one value per sample.
#' @param color Optional color mapping. Supply either a metadata column name, a
#'   vector with one value per sample, or a single constant color value.
#' @param xlimit Optional maximum limit for the x-axis (sequence sample size).
#' @param plot Logical; if TRUE, the function plots the rarefaction curve.
#' @param parallel Logical; if TRUE, uses parallel processing for faster computation.
#' @param se Logical; if TRUE, includes standard error in the rarefaction curves.
#'
#' @return A ggplot object representing the rarefaction curve.
#'         This object is returned invisibly and plotted if 'plot' is TRUE.
#'
#' @examples
#' # Assuming 'ps' is a Phyloseq object
#' Go_rare(ps, step = 10, label = "SampleLabel", color = "Group", xlimit = 1000, plot = TRUE)
#'
#' @export
#' @importFrom phyloseq phyloseq otu_table taxa_are_rows sample_data
#' @importFrom vegan rarefy
#' @importFrom parallel mclapply
#' @importFrom ggplot2 ggplot aes_string labs geom_text geom_line geom_ribbon theme_classic


Go_rare <- function(physeq_object, step = 10, label = NULL, color = NULL, xlimit = NULL, plot = TRUE, parallel = FALSE, se = TRUE) {
  if (!inherits(physeq_object, "phyloseq")) {
    stop("`physeq_object` must be a phyloseq object.")
  }
  
  x <- methods::as(phyloseq::otu_table(physeq_object), "matrix")
  if (phyloseq::taxa_are_rows(physeq_object)) { x <- t(x) }
  
  ## This script is adapted from vegan `rarecurve` function
  tot <- rowSums(x)
  S <- rowSums(x > 0)
  nr <- nrow(x)
  
  rarefun <- function(i) {
    cat(paste("rarefying sample", rownames(x)[i]), sep = "\n")
    n <- seq(1, tot[i], by = step)
    if (n[length(n)] != tot[i]) {
      n <- c(n, tot[i])
    }
    y <- vegan::rarefy(x[i, ,drop = FALSE], n, se = se)
    if (nrow(y) != 1) {
      rownames(y) <- c(".S", ".se")
      return(data.frame(t(y), Size = n, Sample = rownames(x)[i]))
    } else {
      return(data.frame(.S = y[1, ], Size = n, Sample = rownames(x)[i]))
    }
  }
  if (parallel) {
    out <- parallel::mclapply(seq_len(nr), rarefun, mc.preschedule = FALSE)
  } else {
    # out <- lapply(seq_len(nr), rarefun)
    out <- lapply(seq_len(nr), function(x) {
      result <- NULL
      capture.output(result <- rarefun(x))
      return(result)
    })
  }
  
  df <- do.call(rbind, out)
  labels <- data.frame(x = tot, y = S, Sample = rownames(x), stringsAsFactors = FALSE)
  
  # Get sample data
  sample_df <- NULL
  if (!is.null(phyloseq::sample_data(physeq_object, FALSE))) {
    sdf <- methods::as(phyloseq::sample_data(physeq_object), "data.frame")
    sdf$Sample <- rownames(sdf)
    sample_df <- sdf
    data <- merge(df, sdf, by = "Sample", all.x = TRUE)
    labels <- merge(labels, sdf, by = "Sample")
  } else {
    data <- df
  }
  
  sample_order <- rownames(x)

  add_sample_mapping <- function(target_df, spec, column_name) {
    if (is.null(spec)) {
      return(target_df)
    }

    if (length(spec) == 1 && is.character(spec) && !is.null(sample_df) && spec %in% colnames(sample_df)) {
      return(target_df)
    }

    if (length(spec) == 1) {
      target_df[[column_name]] <- spec
      return(target_df)
    }

    if (length(spec) != length(sample_order)) {
      stop(sprintf("`%s` must be a metadata column name, a single value, or a vector with one value per sample.", column_name))
    }

    map_df <- data.frame(
      Sample = sample_order,
      value = spec,
      stringsAsFactors = FALSE
    )
    colnames(map_df)[2] <- column_name
    merge(target_df, map_df, by = "Sample", all.x = TRUE)
  }

  data <- add_sample_mapping(data, color, ".plot_color")
  labels <- add_sample_mapping(labels, label, ".plot_label")

  color_mapping <- NULL
  color_constant <- NULL
  if (!is.null(color)) {
    if (length(color) == 1 && is.character(color) && !is.null(sample_df) && color %in% colnames(sample_df)) {
      color_mapping <- color
    } else if (".plot_color" %in% colnames(data)) {
      if (length(unique(stats::na.omit(data$.plot_color))) <= 1) {
        color_constant <- unique(stats::na.omit(data$.plot_color))[1]
      } else {
        color_mapping <- ".plot_color"
      }
    }
  }

  label_mapping <- NULL
  if (!is.null(label)) {
    if (length(label) == 1 && is.character(label) && !is.null(sample_df) && label %in% colnames(sample_df)) {
      label_mapping <- label
    } else if (".plot_label" %in% colnames(labels)) {
      label_mapping <- ".plot_label"
    }
  }

  base_mapping <- ggplot2::aes_string(x = "Size", y = ".S", group = "Sample")
  if (!is.null(color_mapping)) {
    base_mapping <- ggplot2::aes_string(x = "Size", y = ".S", group = "Sample", color = color_mapping)
  }

  p <- ggplot2::ggplot(data = data, base_mapping)
  p <- p + ggplot2::labs(x = "Sequence Sample Size", y = "Species Richness")
  p <- p + ggplot2::theme_classic()

  if (!is.null(xlimit)) {
    p <- p + ggplot2::coord_cartesian(xlim = c(NA, xlimit))
  }

  if (!is.null(label_mapping)) {
    p <- p + ggplot2::geom_text(data = labels,
                                ggplot2::aes_string(x = "x",
                                                    y = "y",
                                                    label = label_mapping,
                                                    color = color_mapping),
                                size = 4, hjust = 0)
  }

  if (is.null(color_constant)) {
    p <- p + ggplot2::geom_line()
  } else {
    p <- p + ggplot2::geom_line(color = color_constant)
  }

  if (se) { ## add standard error if available
    ribbon_mapping <- if (is.null(color_mapping)) {
      ggplot2::aes_string(ymin = ".S - .se", ymax = ".S + .se", group = "Sample")
    } else {
      ggplot2::aes_string(ymin = ".S - .se", ymax = ".S + .se", group = "Sample", fill = color_mapping)
    }

    if (is.null(color_constant)) {
      p <- p + ggplot2::geom_ribbon(ribbon_mapping, alpha = 0.2, inherit.aes = FALSE)
    } else {
      p <- p + ggplot2::geom_ribbon(ribbon_mapping, alpha = 0.2, fill = color_constant, inherit.aes = FALSE)
    }
  }
  if (plot) {
    suppressWarnings(plot(p))
  }
  invisible(p)
}
