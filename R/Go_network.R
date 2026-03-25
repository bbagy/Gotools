#' Go_network
#'
#' Constructs and visualizes a feature association network from one to three
#' abundance tables, with subgroup filtering from sample metadata.
#'
#' @param tab1_path Path to the CSV file containing the first table of data.
#' @param tab2_path Optional path to the CSV file containing the second table of data.
#' @param tab3_path Optional path to the CSV file containing the third table of data.
#' @param Sampledata Path to the CSV file containing sample metadata, or a data.frame.
#' @param mainGroup The name of the main grouping variable in the sample metadata.
#' @param subgroup The name of the subgroup for filtering within the main group.
#' @param tab1_name Descriptive name for the first data table used in the network visualization.
#' @param tab2_name Optional descriptive name for the second data table.
#' @param tab3_name Optional descriptive name for the third data table.
#' @param project Project name used for output directories and file names.
#' @param cutoff Correlation coefficient threshold for including edges in the network.
#' @param pval_threshold Optional threshold for p-values to filter edges based on statistical significance.
#' @param qval_threshold Optional threshold for FDR-adjusted q-values to filter edges based on statistical significance.
#' @param method Network backend. Currently supports `"clr_spearman"` and `"kendall"`.
#'   `"sparcc"` and `"spring"` are reserved interface values and will stop with an informative message.
#' @param prevalence Minimum fraction of samples with non-zero abundance required to keep a feature.
#' @param min_nonzero Minimum number of non-zero samples required to keep a feature.
#'   If `NULL`, it is derived from `prevalence`.
#' @param pseudo_count Pseudocount added before CLR transformation.
#' @param min_variance Minimum variance required to keep a feature after preprocessing.
#' @param node_font Font size for node labels in the network plot.
#' @param node_names Boolean to indicate whether to display node names or IDs.
#' @param name Optional name for additional naming in output files.
#' @param width Width of the output plot in inches.
#' @param height Height of the output plot in inches.
#'
#' @return Invisibly returns a list containing the filtered matrix, association table,
#'   edge table, node map, and network object. Also saves network plots and tables.
#'
#' @details
#' This refactored version keeps the original "single function to finished output"
#' workflow, but uses a microbiome-oriented preprocessing path:
#' \itemize{
#'   \item subgroup filtering using sample metadata
#'   \item duplicate-column cleanup after merging tables
#'   \item prevalence / non-zero filtering
#'   \item optional CLR transformation for compositional data
#'   \item unique undirected edge testing and BH correction on the final edge set
#' }
#' For microbiome data, the default method is `clr_spearman`, which is generally more
#' appropriate than running rank correlation directly on raw compositional abundances.
#'
#' @examples
#' \dontrun{
#' Go_network(
#'   tab1_path = "path/to/data1.csv",
#'   Sampledata = "path/to/sampledata.csv",
#'   mainGroup = "Treatment",
#'   subgroup = "Control",
#'   tab1_name = "16S",
#'   project = "MyProject",
#'   method = "clr_spearman",
#'   prevalence = 0.1,
#'   qval_threshold = 0.1,
#'   node_font = 0.8
#' )
#' }
#' @export

Go_network <- function(
    tab1_path,
    tab2_path = NULL,
    tab3_path = NULL,
    Sampledata,
    mainGroup,
    subgroup,
    tab1_name,
    tab2_name = NULL,
    tab3_name = NULL,
    project = "Go_network",
    cutoff = 0.3,
    pval_threshold = NULL,
    qval_threshold = NULL,
    method = c("clr_spearman", "kendall", "sparcc", "spring"),
    prevalence = 0.1,
    min_nonzero = NULL,
    pseudo_count = 0.5,
    min_variance = 0,
    node_font = 0.8,
    node_names = TRUE,
    name = NULL,
    width = 10,
    height = 10
) {

  method <- match.arg(method)

  required_pkgs <- c("dplyr", "igraph", "purrr", "RColorBrewer")
  missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_pkgs) > 0) {
    stop("Missing required packages: ", paste(missing_pkgs, collapse = ", "))
  }
  if (identical(method, "clr_spearman") && !requireNamespace("compositions", quietly = TRUE)) {
    stop("`method = 'clr_spearman'` requires the 'compositions' package.")
  }
  if (method %in% c("sparcc", "spring")) {
    stop("`method = '", method, "'` is reserved in this build, but the backend is not yet implemented. ",
         "Use `method = 'clr_spearman'` or `method = 'kendall'` for now.")
  }

  process_table <- function(tab) {
    tab <- as.data.frame(tab, check.names = FALSE, stringsAsFactors = FALSE)
    numeric_cols <- vapply(tab, is.numeric, logical(1))
    if (any(numeric_cols)) {
      tab$RowSum <- rowSums(tab[, numeric_cols, drop = FALSE], na.rm = TRUE)
    } else {
      tab$RowSum <- 0
    }

    if ("Species" %in% colnames(tab)) {
      tab$names <- paste(tab$Species, tab$RowSum, sep = "_")
      tab$names <- make.unique(tab$names)
    } else {
      tab$names <- rownames(tab)
    }

    rownames(tab) <- gsub("\\[|\\]", "", tab$names)

    for (rank in c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "RowSum", "names")) {
      if (rank %in% colnames(tab)) {
        tab[, rank] <- NULL
      }
    }

    as.data.frame(t(tab), stringsAsFactors = FALSE, check.names = FALSE)
  }

  read_and_process <- function(file_path) {
    if (is.null(file_path)) {
      return(NULL)
    }
    if (!file.exists(file_path)) {
      warning(sprintf("File not found: %s", file_path))
      return(NULL)
    }

    tab <- suppressWarnings(
      utils::read.csv(file_path, row.names = 1, check.names = FALSE, fill = TRUE, strip.white = TRUE)
    )

    if (nrow(tab) == 0 || ncol(tab) == 0) {
      warning(sprintf("Warning: %s is empty or malformed. Returning NULL.", file_path))
      return(NULL)
    }

    process_table(tab)
  }

  sanitize_tag <- function(x) {
    gsub("[^A-Za-z0-9._-]+", "-", x)
  }

  merge_feature_tables <- function(tab_list) {
    nonnull <- Filter(Negate(is.null), tab_list)
    if (length(nonnull) == 0) {
      return(NULL)
    }
    merged <- nonnull[[1]]
    if (length(nonnull) == 1) {
      merged$Row.names <- rownames(merged)
      return(merged)
    }
    for (i in 2:length(nonnull)) {
      next_tab <- nonnull[[i]]
      next_tab$Row.names <- rownames(next_tab)
      merged <- merge(merged, next_tab, by = "Row.names", all = TRUE)
    }
    merged
  }

  coalesce_duplicate_columns <- function(dat) {
    dat <- as.data.frame(dat, check.names = FALSE, stringsAsFactors = FALSE)
    base_names <- gsub("\\.(x|y)$", "", colnames(dat))
    unique_names <- unique(base_names)
    out <- lapply(unique_names, function(nm) {
      idx <- which(base_names == nm)
      block <- dat[, idx, drop = FALSE]
      if (ncol(block) == 1) {
        return(block[[1]])
      }
      block_num <- suppressWarnings(as.data.frame(lapply(block, as.numeric), check.names = FALSE))
      conflict_rows <- apply(block_num, 1, function(x) {
        x <- x[!is.na(x)]
        length(unique(x)) > 1
      })
      if (any(conflict_rows)) {
        warning(sprintf("Duplicate feature '%s' had conflicting values across merged tables; first non-missing value was kept.", nm))
      }
      dplyr::coalesce(!!!block_num)
    })
    out <- as.data.frame(out, check.names = FALSE, stringsAsFactors = FALSE)
    colnames(out) <- unique_names
    out
  }

  prepare_matrix <- function(dat, prevalence, min_nonzero, min_variance) {
    mat <- as.matrix(dat)
    storage.mode(mat) <- "numeric"
    mat[is.na(mat)] <- 0

    if (!is.numeric(prevalence) || length(prevalence) != 1 || prevalence < 0 || prevalence > 1) {
      stop("`prevalence` must be a single number between 0 and 1.")
    }
    if (is.null(min_nonzero)) {
      min_nonzero <- max(1, ceiling(nrow(mat) * prevalence))
    }

    prevalence_prop <- colMeans(mat > 0, na.rm = TRUE)
    keep_prev <- prevalence_prop >= prevalence
    keep_nonzero <- colSums(mat > 0, na.rm = TRUE) >= min_nonzero
    keep_var <- apply(mat, 2, stats::var, na.rm = TRUE) > min_variance
    keep <- keep_prev & keep_nonzero & keep_var

    mat <- mat[, keep, drop = FALSE]
    list(
      matrix = mat,
      prevalence = prevalence_prop,
      min_nonzero = min_nonzero
    )
  }

  clr_transform <- function(mat, pseudo_count) {
    mat_pos <- mat + pseudo_count
    row_totals <- rowSums(mat_pos, na.rm = TRUE)
    if (any(!is.finite(row_totals)) || any(row_totals <= 0)) {
      stop("CLR transformation failed because one or more samples had non-positive total abundance after pseudocount addition.")
    }
    comp <- sweep(mat_pos, 1, row_totals, "/")
    compositions::clr(comp)
  }

  compute_association_table <- function(mat, method) {
    features <- colnames(mat)
    if (length(features) < 2) {
      return(data.frame(
        Source = character(0),
        Target = character(0),
        Correlation = numeric(0),
        p_value = numeric(0),
        q_value = numeric(0),
        stringsAsFactors = FALSE
      ))
    }

    pairs <- utils::combn(features, 2, simplify = FALSE)
    assoc_list <- lapply(pairs, function(pair) {
      test_method <- if (identical(method, "kendall")) "kendall" else "spearman"
      test_result <- suppressWarnings(
        stats::cor.test(mat[, pair[1]], mat[, pair[2]], method = test_method, exact = FALSE)
      )
      data.frame(
        Source = pair[1],
        Target = pair[2],
        Correlation = unname(test_result$estimate),
        p_value = test_result$p.value,
        stringsAsFactors = FALSE
      )
    })

    cor_results <- dplyr::bind_rows(assoc_list)
    cor_results$q_value <- stats::p.adjust(cor_results$p_value, method = "fdr")
    cor_results
  }

  build_edges <- function(cor_results, cutoff, pval_threshold, qval_threshold) {
    edges <- cor_results %>%
      dplyr::filter(is.finite(Correlation), abs(Correlation) >= cutoff)

    sig <- "none"
    sigval <- sprintf("|assoc| >= %.2f", cutoff)
    signame <- "FDR"

    if (!is.null(qval_threshold)) {
      edges <- edges %>% dplyr::filter(q_value < qval_threshold)
      sig <- "FDR"
      sigval <- sprintf("FDR < %s", qval_threshold)
    } else if (!is.null(pval_threshold)) {
      edges <- edges %>% dplyr::filter(p_value < pval_threshold)
      sig <- "p"
      sigval <- sprintf("p < %s", pval_threshold)
    } else {
      sig <- "p_FDR"
      sigval <- "FDR highlighted"
      signame <- "FDR"
    }

    list(edges = edges, sig = sig, sigval = sigval, signame = signame)
  }

  build_node_colors <- function(network, tab1, tab2, tab3, tab1_name, tab2_name, tab3_name) {
    if (!is.null(tab1) && is.null(tab2) && is.null(tab3)) {
      all_bacteria_tab1 <- colnames(tab1)
      genus_names <- sub("[ _.].*", "", all_bacteria_tab1)
      unique_genera <- unique(genus_names)
      num_genera <- length(unique_genera)

      if (num_genera <= 1) {
        genus_colors <- stats::setNames("deepskyblue", unique_genera)
      } else if (num_genera > 9) {
        genus_colors <- stats::setNames(grDevices::colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))(num_genera),
                                        unique_genera)
      } else {
        genus_colors <- stats::setNames(RColorBrewer::brewer.pal(num_genera, "Set1"), unique_genera)
      }

      igraph::V(network)$Genus <- genus_names[match(igraph::V(network)$name, all_bacteria_tab1)]
      node_colors <- vapply(igraph::V(network)$Genus, function(genus) {
        genus_colors[[genus]] %||% "gray"
      }, character(1))

      return(list(
        node_colors = node_colors,
        legend_labels = unique_genera,
        legend_colors = genus_colors[unique_genera],
        legend_title = "Genus Groups"
      ))
    }

    color_mapping <- stats::setNames(
      c("deepskyblue", "gold", "forestgreen", "gray"),
      c(tab1_name, tab2_name %||% "Table2", tab3_name %||% "Table3", "Unknown")
    )
    node_colors <- vapply(igraph::V(network)$Type, function(type) {
      color_mapping[[type]] %||% "gray"
    }, character(1))
    used_types <- unique(igraph::V(network)$Type)
    legend_labels <- used_types[used_types %in% names(color_mapping)]
    if (length(legend_labels) == 0) {
      legend_labels <- "Unknown"
    }
    legend_colors <- vapply(legend_labels, function(type) color_mapping[[type]] %||% "gray", character(1))

    list(
      node_colors = node_colors,
      legend_labels = legend_labels,
      legend_colors = legend_colors,
      legend_title = "Node Types"
    )
  }

  if (is.character(Sampledata)) {
    sampledata <- utils::read.csv(Sampledata, row.names = 1, check.names = FALSE)
  } else {
    sampledata <- as.data.frame(Sampledata, check.names = FALSE, stringsAsFactors = FALSE)
  }

  if (!mainGroup %in% colnames(sampledata)) {
    stop("`mainGroup` not found in `Sampledata`.")
  }

  tab1 <- read_and_process(tab1_path)
  if (is.null(tab1)) {
    stop("`tab1_path` is required and must point to a readable table.")
  }
  tab2 <- read_and_process(tab2_path)
  tab3 <- read_and_process(tab3_path)

  common_samples <- rownames(tab1)
  if (!is.null(tab2)) common_samples <- intersect(common_samples, rownames(tab2))
  if (!is.null(tab3)) common_samples <- intersect(common_samples, rownames(tab3))
  common_samples <- intersect(common_samples, rownames(sampledata))

  if (length(common_samples) < 3) {
    stop("Fewer than 3 common samples remained after aligning the input tables and sample metadata.")
  }

  tab1 <- tab1[common_samples, , drop = FALSE]
  if (!is.null(tab2)) tab2 <- tab2[common_samples, , drop = FALSE]
  if (!is.null(tab3)) tab3 <- tab3[common_samples, , drop = FALSE]
  sampledata <- sampledata[common_samples, , drop = FALSE]
  sampledata$Row.names <- rownames(sampledata)

  merged_table <- merge_feature_tables(list(tab1, tab2, tab3))
  final_merged_table <- merge(
    merged_table,
    sampledata[, c("Row.names", mainGroup), drop = FALSE],
    by = "Row.names",
    all = FALSE
  )
  rownames(final_merged_table) <- final_merged_table$Row.names

  all_features <- unique(c(colnames(tab1), colnames(tab2), colnames(tab3)))
  all_features <- all_features[!is.na(all_features)]
  bacteria_map <- data.frame(
    Bacteria = all_features,
    ID = as.character(seq_along(all_features)),
    stringsAsFactors = FALSE
  )

  message(sprintf("[Go_network] total features before filtering: %d", length(all_features)))

  data <- final_merged_table %>%
    dplyr::filter(.data[[mainGroup]] == subgroup) %>%
    dplyr::select(-dplyr::all_of(mainGroup), -Row.names)

  if (nrow(data) < 3) {
    stop("Fewer than 3 samples remained in the selected subgroup.")
  }

  data <- as.data.frame(lapply(data, function(x) suppressWarnings(as.numeric(trimws(as.character(x))))),
                        check.names = FALSE, stringsAsFactors = FALSE)
  data <- coalesce_duplicate_columns(data)

  prep <- prepare_matrix(data, prevalence = prevalence, min_nonzero = min_nonzero, min_variance = min_variance)
  data_mat <- prep$matrix
  min_nonzero_used <- prep$min_nonzero

  if (ncol(data_mat) < 2) {
    stop("Fewer than 2 features remained after prevalence / variance filtering.")
  }

  assoc_input <- switch(
    method,
    clr_spearman = clr_transform(data_mat, pseudo_count = pseudo_count),
    kendall = data_mat
  )

  cor_results <- compute_association_table(assoc_input, method = method)
  edge_info <- build_edges(cor_results, cutoff = cutoff, pval_threshold = pval_threshold, qval_threshold = qval_threshold)
  edges_unique <- edge_info$edges
  sig <- edge_info$sig
  sigval <- edge_info$sigval
  signame <- edge_info$signame

  feature_presence <- unique(c(edges_unique$Source, edges_unique$Target))
  if (length(feature_presence) == 0) {
    warning("No edges passed the current cutoff / significance criteria. Returning tables without plotting a network.")
  }

  nodes <- data.frame(name = feature_presence, Type = NA_character_, stringsAsFactors = FALSE)
  if (!is.null(tab1)) nodes$Type[nodes$name %in% colnames(tab1)] <- tab1_name
  if (!is.null(tab2)) nodes$Type[nodes$name %in% colnames(tab2)] <- tab2_name
  if (!is.null(tab3) && !is.null(tab3_name)) nodes$Type[nodes$name %in% colnames(tab3)] <- tab3_name
  nodes$Type[is.na(nodes$Type)] <- "Unknown"

  network <- igraph::graph_from_data_frame(
    d = edges_unique[, c("Source", "Target", "Correlation", "p_value", "q_value"), drop = FALSE],
    vertices = nodes,
    directed = FALSE
  )
  if (igraph::gsize(network) > 0) {
    loop_indices <- which(igraph::which_loop(network))
    if (length(loop_indices) > 0) {
      network <- igraph::delete_edges(network, loop_indices)
    }
  }

  out <- file.path(sprintf("%s_%s", project, format(Sys.Date(), "%y%m%d")))
  out_pdf <- file.path(out, "pdf")
  out_table <- file.path(out, "table")
  out_network <- file.path(out_table, "network")
  for (d in c(out, out_pdf, out_table, out_network)) {
    if (!dir.exists(d)) dir.create(d, recursive = TRUE)
  }

  filter_summary <- data.frame(
    subgroup = subgroup,
    method = method,
    n_samples = nrow(data_mat),
    n_features_before = ncol(data),
    n_features_after = ncol(data_mat),
    prevalence = prevalence,
    min_nonzero = min_nonzero_used,
    pseudo_count = pseudo_count,
    cutoff = cutoff,
    stringsAsFactors = FALSE
  )

  utils::write.csv(filter_summary,
                   file.path(out_network, sprintf("network.%s.%s.filter_summary.%s.csv",
                                                  mainGroup, sanitize_tag(subgroup), format(Sys.Date(), "%y%m%d"))),
                   row.names = FALSE)
  utils::write.csv(as.data.frame(data_mat, check.names = FALSE),
                   file.path(out_network, sprintf("network.%s.%s.filtered_matrix.%s.csv",
                                                  mainGroup, sanitize_tag(subgroup), format(Sys.Date(), "%y%m%d"))))
  utils::write.csv(cor_results,
                   file.path(out_network, sprintf("network.%s.%s.associations.%s.csv",
                                                  mainGroup, sanitize_tag(subgroup), format(Sys.Date(), "%y%m%d"))),
                   row.names = FALSE)
  utils::write.csv(edges_unique,
                   file.path(out_network, sprintf("network.%s.%s.edges.%s.csv",
                                                  mainGroup, sanitize_tag(subgroup), format(Sys.Date(), "%y%m%d"))),
                   row.names = FALSE)

  if (igraph::gsize(network) > 0) {
    color_info <- build_node_colors(network, tab1, tab2, tab3, tab1_name, tab2_name, tab3_name)
    node_colors <- color_info$node_colors
    legend_labels <- color_info$legend_labels
    legend_colors <- color_info$legend_colors
    legend_title <- color_info$legend_title

    bacteria_map_filtered <- bacteria_map %>%
      dplyr::filter(Bacteria %in% igraph::V(network)$name)

    node_labels <- if (node_names) {
      bacteria_map_filtered$ID[match(igraph::V(network)$name, bacteria_map_filtered$Bacteria)]
    } else {
      igraph::V(network)$name
    }

    edge_styles <- if (sig == "FDR") {
      ifelse(igraph::E(network)$q_value < (qval_threshold %||% 0.05), 1, 2)
    } else if (sig == "p") {
      ifelse(igraph::E(network)$p_value < (pval_threshold %||% 0.05), 1, 2)
    } else {
      ifelse(igraph::E(network)$q_value < 0.05, 1, 2)
    }

    if (!is.null(dev.list())) dev.off()

    pdf(
      file.path(
        out_pdf,
        sprintf(
          "network.%s.%s.(%s).(%s).%s.%s.%s%s%s%s%s.pdf",
          mainGroup,
          sanitize_tag(subgroup),
          method,
          cutoff,
          sig,
          ifelse(node_names, "IDs", "Names"),
          ifelse(is.null(tab1_name), "", paste0(sanitize_tag(tab1_name), ".")),
          ifelse(is.null(tab2_name), "", paste0(sanitize_tag(tab2_name), ".")),
          ifelse(is.null(tab3_name), "", paste0(sanitize_tag(tab3_name), ".")),
          ifelse(is.null(name), "", paste0(sanitize_tag(name), ".")),
          format(Sys.Date(), "%y%m%d")
        )
      ),
      width = width,
      height = height
    )

    set.seed(123)
    plot(
      network,
      vertex.color = node_colors,
      vertex.label = node_labels,
      vertex.label.cex = node_font,
      vertex.size = 5 + 10 * (igraph::degree(network) / max(igraph::degree(network))),
      edge.width = abs(igraph::E(network)$Correlation) * 5,
      edge.color = ifelse(igraph::E(network)$Correlation > 0, "red", "blue"),
      edge.lty = edge_styles,
      edge.curved = 0.15,
      main = sprintf("%s-%s Network (%s; %s)", mainGroup, subgroup, method, sigval)
    )

    legend("bottomright",
           legend = legend_labels,
           col = legend_colors,
           pch = 19,
           pt.cex = 1.5,
           bty = "n",
           title = legend_title)

    legend("bottomleft",
           legend = c("Positive association", "Negative association"),
           col = c("red", "blue"),
           lty = 1,
           lwd = 3,
           bty = "n",
           title = "Edge sign")

    if (node_names) {
      legend("topleft",
             legend = paste(bacteria_map_filtered$ID, bacteria_map_filtered$Bacteria, sep = " = "),
             cex = 0.6,
             bty = "n")
    }

    if (sig == "p_FDR") {
      legend("topright",
             legend = c(sprintf("%s < 0.05", signame), sprintf("%s >= 0.05", signame)),
             lty = c(1, 2),
             lwd = 3,
             bty = "n",
             title = sprintf("%s Significance", signame))
    }

    dev.off()
  }

  if (igraph::gsize(network) > 0) {
    centrality_data <- data.frame(
      Node = igraph::V(network)$name,
      Degree = igraph::degree(network),
      Betweenness = igraph::betweenness(network),
      Closeness = igraph::closeness(network),
      Eigenvector = igraph::eigen_centrality(network)$vector,
      stringsAsFactors = FALSE
    )
  } else {
    centrality_data <- data.frame(
      Node = character(0),
      Degree = numeric(0),
      Betweenness = numeric(0),
      Closeness = numeric(0),
      Eigenvector = numeric(0),
      stringsAsFactors = FALSE
    )
  }

  sorted_centrality <- centrality_data[order(-centrality_data$Degree), , drop = FALSE]
  utils::write.csv(sorted_centrality,
                   file.path(out_network, sprintf("%s.%s.sorted_centrality.%s.csv",
                                                  mainGroup, sig, format(Sys.Date(), "%y%m%d"))),
                   row.names = FALSE)

  invisible(list(
    filtered_matrix = data_mat,
    association_table = cor_results,
    edge_table = edges_unique,
    bacteria_map = bacteria_map,
    network = network,
    centrality = sorted_centrality
  ))
}
