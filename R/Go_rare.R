#' Generate Rarefaction Curves for Microbiome Data
#'
#' This function creates rarefaction curves for microbiome data using a Phyloseq object.
#' It provides an assessment of species richness across varying sequencing depths.
#'
#' @param physeq_object A Phyloseq object containing OTU/ASV counts and associated sample data.
#' @param step Step size for rarefaction curve plotting.
#' @param label Optional labels for the curves.
#' @param color Color to be used for the curves.
#' @param xlimit Maximum limit for the x-axis (sequence sample size).
#' @param plot Logical; if TRUE, the function plots the rarefaction curve.
#' @param parallel Logical; if TRUE, uses parallel processing for faster computation.
#' @param se Logical; if TRUE, includes standard error in the rarefaction curves.
#'
#' @return A ggplot object representing the rarefaction curve.
#'         This object is returned invisibly and plotted if 'plot' is TRUE.
#'
#' @examples
#' # Assuming 'ps' is a Phyloseq object
#' Go_rare_markdown(ps, step = 10, label = "SampleLabel", color = "red", xlimit = 1000, plot = TRUE)
#'
#' @export
#' @importFrom phyloseq phyloseq otu_table taxa_are_rows sample_data
#' @importFrom vegan rarefy
#' @importFrom parallel mclapply
#' @importFrom ggplot2 ggplot aes_string labs geom_text geom_line geom_ribbon theme_classic


Go_rare <- function(physeq_object, step = 10, label = NULL, color = NULL, xlimit, plot = TRUE, parallel = FALSE, se = TRUE) {
  
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
  
  # Get sample data
  if (!is.null(phyloseq::sample_data(physeq_object, FALSE))) {
    sdf <- methods::as(phyloseq::sample_data(physeq_object), "data.frame")
    sdf$Sample <- rownames(sdf)
    data <- merge(df, sdf, by = "Sample")
    labels <- data.frame(x = tot, y = S, Sample = rownames(x))
    labels <- merge(labels, sdf, by = "Sample")
  }
  
  # Add, any custom-supplied plot-mapped variables
  if ( length(color) > 1 ) {
    data$color <- color
    names(data)[names(data) == "color"] <- deparse(substitute(color))
    color <- deparse(substitute(color))
  }
  
  if ( length(label) > 1 ) {
    labels$label <- label
    names(labels)[names(labels) == "label"] <- deparse(substitute(label))
    label <- deparse(substitute(label))
  }
  
  p <- ggplot2::ggplot(data = data,
                       ggplot2::aes_string(x = "Size",
                                           y = ".S",
                                           group = "Sample",
                                           color = color))
  
  p <- p + ggplot2::labs(x = "Sequence Sample Size", y = "Species Richness")
  p <- p + xlim(NA, xlimit) + theme_classic()
  
  if (!is.null(label)) {
    p <- p + ggplot2::geom_text(data = labels,
                                ggplot2::aes_string(x = "x",
                                                    y = "y",
                                                    label = label,
                                                    color = color),
                                size = 4, hjust = 0)
  }
  
  p <- p + ggplot2::geom_line()
  if (se) { ## add standard error if available
    p <- p +
      ggplot2::geom_ribbon(ggplot2::aes_string(ymin = ".S - .se",
                                               ymax = ".S + .se",
                                               color = NULL,
                                               fill = color),
                           alpha = 0.2)
  }
  if (plot) {
    suppressWarnings(plot(p))
  }
  invisible(p)
}
