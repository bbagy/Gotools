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
#' @param subgroup A single level of `mainGroup` to plot, or `NULL` (default) to
#'   plot all levels side by side in one PDF with a shared legend.
#' @param tab1_name Descriptive name for the first data table used in the network visualization.
#' @param tab2_name Optional descriptive name for the second data table.
#' @param tab3_name Optional descriptive name for the third data table.
#' @param project Project name used for output directories and file names.
#' @param cutoff Correlation coefficient threshold for including edges in the network.
#' @param pval_threshold Optional threshold for p-values to filter edges based on statistical significance.
#' @param qval_threshold Optional threshold for FDR-adjusted q-values to filter edges based on statistical significance.
#' @param method Network backend. Supported methods are `"clr_spearman"` and
#'   `"clr_kendall"`. Both apply CLR transformation before association testing,
#'   which is generally more appropriate for compositional microbiome data than
#'   correlating raw relative abundances.
#' @param prevalence Minimum fraction of samples with non-zero abundance required to keep a feature.
#' @param min_nonzero Minimum number of non-zero samples required to keep a feature.
#'   If `NULL`, it is derived from `prevalence`.
#' @param pseudo_count Pseudocount added before CLR transformation.
#' @param min_variance Minimum variance required to keep a feature after preprocessing.
#' @param node_font Font size for node labels in the network plot.
#' @param node_names Boolean indicating whether to display node names (`TRUE`)
#'   or numeric IDs (`FALSE`).
#' @param name Optional name for additional naming in output files.
#' @param width Width of the output plot in inches.
#' @param height Height of the output plot in inches.
#'
#' @return Invisibly returns a list containing the filtered matrix, association
#'   input used for testing, method summary, association table, edge table, node
#'   map, and network object. Also saves network plots and tables.
#'
#' @details
#' This refactored version keeps the original "single function to finished output"
#' workflow, but uses a microbiome-oriented preprocessing path:
#' \itemize{
#'   \item subgroup filtering using sample metadata
#'   \item duplicate-column cleanup after merging tables
#'   \item prevalence / non-zero filtering
#'   \item CLR transformation for compositional data before association testing
#'   \item unique undirected edge testing and BH correction on the final edge set
#' }
#' For microbiome data, the default method is `clr_spearman`, which is usually
#' the most straightforward and reviewer-friendly choice for compositional 16S
#' abundance data.
#'
#' @examples
#' \dontrun{
#' # Plot all subgroups side by side (subgroup = NULL is the default)
#' Go_network(
#'   tab1_path = "path/to/data1.csv",
#'   Sampledata = "path/to/sampledata.csv",
#'   mainGroup = "VRE_byID",
#'   tab1_name = "16S",
#'   project = "MyProject",
#'   method = "clr_spearman",
#'   prevalence = 0.1,
#'   qval_threshold = 0.1
#' )
#'
#' # Plot a single subgroup
#' Go_network(
#'   tab1_path = "path/to/data1.csv",
#'   Sampledata = "path/to/sampledata.csv",
#'   mainGroup = "VRE_byID",
#'   subgroup = "Yes",
#'   tab1_name = "16S",
#'   project = "MyProject",
#'   method = "clr_spearman",
#'   prevalence = 0.1,
#'   qval_threshold = 0.1
#' )
#' }
#' @export

Go_network <- function(
    tab1_path,
    tab2_path = NULL,
    tab3_path = NULL,
    Sampledata,
    mainGroup,
    subgroup = NULL,
    tab1_name,
    tab2_name = NULL,
    tab3_name = NULL,
    project = "Go_network",
    cutoff = 0.3,
    pval_threshold = NULL,
    qval_threshold = NULL,
    method = c("clr_spearman", "clr_kendall"),
    prevalence = 0.1,
    min_nonzero = NULL,
    pseudo_count = 0.5,
    min_variance = 0,
    node_font = 0.8,
    node_names = FALSE,
    target_bacteria = NULL,
    name = NULL,
    seed = 123,
    width = 10,
    height = 10
) {

  method <- match.arg(method)

  required_pkgs <- c("dplyr", "igraph", "RColorBrewer")
  missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_pkgs) > 0) {
    stop("Missing required packages: ", paste(missing_pkgs, collapse = ", "))
  }
  if (!requireNamespace("compositions", quietly = TRUE)) {
    stop("CLR-based network inference requires the 'compositions' package.")
  }
  resolve_method <- function(method) {
    mapping <- switch(
      method,
      clr_spearman = list(analysis_method = "clr_spearman", transform = "clr", correlation = "spearman"),
      clr_kendall = list(analysis_method = "clr_kendall", transform = "clr", correlation = "kendall")
    )

    if (is.null(mapping)) {
      stop("Unsupported method: ", method)
    }

    mapping
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

  sample_key <- "..gotools_sample_id.."

  merge_feature_tables <- function(tab_list) {
    nonnull <- Filter(Negate(is.null), tab_list)
    if (length(nonnull) == 0) {
      return(NULL)
    }
    merged <- nonnull[[1]]
    merged[[sample_key]] <- rownames(merged)
    if (length(nonnull) == 1) {
      return(merged)
    }
    for (i in 2:length(nonnull)) {
      next_tab <- nonnull[[i]]
      next_tab[[sample_key]] <- rownames(next_tab)
      merged <- merge(merged, next_tab, by = sample_key, all = TRUE)
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

  compute_association_table <- function(mat, correlation_method) {
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
      test_result <- suppressWarnings(
        stats::cor.test(mat[, pair[1]], mat[, pair[2]], method = correlation_method, exact = FALSE)
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

  method_info <- resolve_method(method)
  analysis_method <- method_info$analysis_method
  transform_method <- method_info$transform
  correlation_method <- method_info$correlation

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
        genus_colors <- stats::setNames("#66C2A5", unique_genera)
      } else if (num_genera > 8) {
        genus_colors <- stats::setNames(grDevices::colorRampPalette(RColorBrewer::brewer.pal(8, "Set2"))(num_genera),
                                        unique_genera)
      } else {
        genus_colors <- stats::setNames(RColorBrewer::brewer.pal(max(3, num_genera), "Set2")[seq_len(num_genera)],
                                        unique_genera)
      }

      igraph::V(network)$Genus <- genus_names[match(igraph::V(network)$name, all_bacteria_tab1)]
      node_colors <- vapply(igraph::V(network)$Genus, function(genus) {
        genus_colors[[genus]] %||% "#BEBEBE"
      }, character(1))

      return(list(
        node_colors = node_colors,
        legend_labels = unique_genera,
        legend_colors = genus_colors[unique_genera],
        legend_title = "Genus Groups"
      ))
    }

    type_names <- c(tab1_name, tab2_name %||% "Table2", tab3_name %||% "Table3", "Unknown")
    n_types <- length(type_names)
    type_palette <- c("#7FCDBB", "#FDB97D", "#C994C7", "#BEBEBE")
    color_mapping <- stats::setNames(type_palette[seq_len(n_types)], type_names)

    node_colors <- vapply(igraph::V(network)$Type, function(type) {
      color_mapping[[type]] %||% "#BEBEBE"
    }, character(1))
    used_types <- unique(igraph::V(network)$Type)
    legend_labels <- used_types[used_types %in% names(color_mapping)]
    if (length(legend_labels) == 0) {
      legend_labels <- "Unknown"
    }
    legend_colors <- vapply(legend_labels, function(type) color_mapping[[type]] %||% "#BEBEBE", character(1))

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
  sampledata[[sample_key]] <- rownames(sampledata)

  merged_table <- merge_feature_tables(list(tab1, tab2, tab3))
  final_merged_table <- merge(
    merged_table,
    sampledata[, c(sample_key, mainGroup), drop = FALSE],
    by = sample_key,
    all = FALSE
  )
  rownames(final_merged_table) <- final_merged_table[[sample_key]]

  all_features <- unique(c(colnames(tab1), colnames(tab2), colnames(tab3)))
  all_features <- all_features[!is.na(all_features)]
  bacteria_map <- data.frame(
    Bacteria = all_features,
    ID = as.character(seq_along(all_features)),
    stringsAsFactors = FALSE
  )

  message(sprintf("[Go_network] total features before filtering: %d", length(all_features)))

  # ── Determine subgroup to plot (always single) ────────────────────────────
  subgroups_to_plot <- if (!is.null(subgroup)) subgroup else "all"

  # ── Create output directories ──────────────────────────────────────────────
  out         <- file.path(sprintf("%s_%s", project, format(Sys.Date(), "%y%m%d")))
  out_pdf     <- file.path(out, "pdf")
  out_table   <- file.path(out, "table")
  out_network <- file.path(out_table, "network")
  for (d in c(out, out_pdf, out_table, out_network)) {
    if (!dir.exists(d)) dir.create(d, recursive = TRUE)
  }

  # ── Single-subgroup computation ───────────────────────────────────────────
  sg_results <- lapply(subgroups_to_plot, function(sg) {
    if (sg == "all") {
      data <- final_merged_table %>%
        dplyr::select(-dplyr::all_of(mainGroup), -dplyr::all_of(sample_key))
    } else {
      data <- final_merged_table %>%
        dplyr::filter(.data[[mainGroup]] == sg) %>%
        dplyr::select(-dplyr::all_of(mainGroup), -dplyr::all_of(sample_key))
    }

    if (nrow(data) < 3) {
      warning(sprintf("[Go_network] Fewer than 3 samples in subgroup '%s'. Skipping.", sg))
      return(NULL)
    }

    data <- as.data.frame(
      lapply(data, function(x) suppressWarnings(as.numeric(trimws(as.character(x))))),
      check.names = FALSE, stringsAsFactors = FALSE
    )
    data <- coalesce_duplicate_columns(data)

    prep             <- prepare_matrix(data, prevalence = prevalence, min_nonzero = min_nonzero, min_variance = min_variance)
    data_mat         <- prep$matrix
    min_nonzero_used <- prep$min_nonzero

    if (ncol(data_mat) < 2) {
      warning(sprintf("[Go_network] Fewer than 2 features after filtering in subgroup '%s'. Skipping.", sg))
      return(NULL)
    }

    assoc_input <- switch(
      transform_method,
      clr  = clr_transform(data_mat, pseudo_count = pseudo_count),
      none = data_mat
    )

    cor_results  <- compute_association_table(assoc_input, correlation_method = correlation_method)
    edge_info    <- build_edges(cor_results, cutoff = cutoff,
                                pval_threshold = pval_threshold, qval_threshold = qval_threshold)
    edges_unique <- edge_info$edges
    sig          <- edge_info$sig
    sigval       <- edge_info$sigval
    signame      <- edge_info$signame

    feature_presence <- unique(c(edges_unique$Source, edges_unique$Target))
    if (length(feature_presence) == 0) {
      warning(sprintf("[Go_network] No edges passed thresholds for subgroup '%s'.", sg))
    }

    nodes <- data.frame(name = feature_presence, Type = NA_character_, stringsAsFactors = FALSE)
    if (!is.null(tab1)) nodes$Type[nodes$name %in% colnames(tab1)] <- tab1_name
    if (!is.null(tab2)) nodes$Type[nodes$name %in% colnames(tab2)] <- tab2_name
    if (!is.null(tab3) && !is.null(tab3_name)) nodes$Type[nodes$name %in% colnames(tab3)] <- tab3_name
    nodes$Type[is.na(nodes$Type)] <- "Unknown"

    network <- igraph::graph_from_data_frame(
      d        = edges_unique[, c("Source", "Target", "Correlation", "p_value", "q_value"), drop = FALSE],
      vertices = nodes,
      directed = FALSE
    )
    if (igraph::gsize(network) > 0) {
      loop_idx <- which(igraph::which_loop(network))
      if (length(loop_idx) > 0) network <- igraph::delete_edges(network, loop_idx)
    }

    filter_summary <- data.frame(
      subgroup           = sg,
      method_requested   = method,
      method_used        = analysis_method,
      transformation     = transform_method,
      correlation_method = correlation_method,
      n_samples          = nrow(data_mat),
      n_features_before  = ncol(data),
      n_features_after   = ncol(data_mat),
      prevalence         = prevalence,
      min_nonzero        = min_nonzero_used,
      pseudo_count       = pseudo_count,
      cutoff             = cutoff,
      stringsAsFactors   = FALSE
    )

    utils::write.csv(filter_summary,
                     file.path(out_network, sprintf("network.%s.%s.filter_summary.%s.csv",
                                                    mainGroup, sanitize_tag(sg), format(Sys.Date(), "%y%m%d"))),
                     row.names = FALSE)
    utils::write.csv(as.data.frame(data_mat, check.names = FALSE),
                     file.path(out_network, sprintf("network.%s.%s.filtered_matrix.%s.csv",
                                                    mainGroup, sanitize_tag(sg), format(Sys.Date(), "%y%m%d"))))
    utils::write.csv(as.data.frame(assoc_input, check.names = FALSE),
                     file.path(out_network, sprintf("network.%s.%s.association_input.%s.csv",
                                                    mainGroup, sanitize_tag(sg), format(Sys.Date(), "%y%m%d"))))
    utils::write.csv(cor_results,
                     file.path(out_network, sprintf("network.%s.%s.associations.%s.csv",
                                                    mainGroup, sanitize_tag(sg), format(Sys.Date(), "%y%m%d"))),
                     row.names = FALSE)
    utils::write.csv(edges_unique,
                     file.path(out_network, sprintf("network.%s.%s.edges.%s.csv",
                                                    mainGroup, sanitize_tag(sg), format(Sys.Date(), "%y%m%d"))),
                     row.names = FALSE)

    list(
      sg             = sg,
      data_mat       = data_mat,
      assoc_input    = assoc_input,
      cor_results    = cor_results,
      edges_unique   = edges_unique,
      sig            = sig,
      sigval         = sigval,
      signame        = signame,
      network        = network,
      filter_summary = filter_summary
    )
  })

  sg_results <- Filter(Negate(is.null), sg_results)

  if (length(sg_results) == 0) {
    warning("[Go_network] No subgroups produced a valid network. No PDF generated.")
    return(invisible(NULL))
  }

  # ── Shared bacteria map (union of all networks' nodes) ─────────────────────
  all_net_nodes <- unique(unlist(lapply(sg_results, function(r) {
    if (igraph::gsize(r$network) > 0) igraph::V(r$network)$name else character(0)
  })))
  bacteria_map_filtered <- bacteria_map %>%
    dplyr::filter(Bacteria %in% all_net_nodes)

  # ── Build per-network visual parameters ───────────────────────────────────
  sg_visuals <- lapply(sg_results, function(res) {
    net <- res$network
    if (igraph::gsize(net) == 0) {
      return(list(color_info = NULL, node_colors = NULL, node_sizes = NULL,
                  frame_colors = NULL, label_fonts = NULL, node_labels = NULL,
                  edge_styles = NULL, norm_lay = NULL,
                  has_target = FALSE, target_idx = NULL))
    }

    color_info  <- build_node_colors(net, tab1, tab2, tab3, tab1_name, tab2_name, tab3_name)
    node_colors <- color_info$node_colors

    node_sizes   <- 5 + 10 * (igraph::degree(net) / max(igraph::degree(net)))
    frame_colors <- rep("white", igraph::vcount(net))
    label_fonts  <- rep(3L, igraph::vcount(net))
    has_target   <- FALSE
    target_idx   <- NULL

    if (!is.null(target_bacteria) && length(target_bacteria) > 0) {
      tidx <- igraph::V(net)$name %in% target_bacteria
      if (any(tidx)) {
        has_target            <- TRUE
        target_idx            <- tidx
        frame_colors[tidx]    <- "#E67E22"
        node_sizes[tidx]      <- node_sizes[tidx] * 1.5
        label_fonts[tidx]     <- 4L
      }
    }

    node_labels <- if (node_names) {
      igraph::V(net)$name
    } else {
      bacteria_map_filtered$ID[match(igraph::V(net)$name, bacteria_map_filtered$Bacteria)]
    }

    edge_styles <- if (res$sig == "FDR") {
      ifelse(igraph::E(net)$q_value < (qval_threshold %||% 0.05), 1, 2)
    } else if (res$sig == "p") {
      ifelse(igraph::E(net)$p_value < (pval_threshold %||% 0.05), 1, 2)
    } else {
      ifelse(igraph::E(net)$q_value < 0.05, 1, 2)
    }

    set.seed(seed)
    raw_lay  <- igraph::layout_with_fr(net)
    norm_lay <- {
      x  <- raw_lay[, 1]; y  <- raw_lay[, 2]
      rx <- range(x);     ry <- range(y)
      if (diff(rx) > 0) x <- (x - rx[1]) / diff(rx) * 2 - 1 else x[] <- 0
      if (diff(ry) > 0) y <- (y - ry[1]) / diff(ry) * 2 - 1 else y[] <- 0
      cbind(x, y)
    }

    list(
      color_info   = color_info,
      node_colors  = node_colors,
      node_sizes   = node_sizes,
      frame_colors = frame_colors,
      label_fonts  = label_fonts,
      node_labels  = node_labels,
      edge_styles  = edge_styles,
      norm_lay     = norm_lay,
      has_target   = has_target,
      target_idx   = target_idx
    )
  })

  # ── PDF output ─────────────────────────────────────────────────────────────
  any_edges <- any(vapply(sg_results, function(r) igraph::gsize(r$network) > 0, logical(1)))

  if (any_edges) {
    n_sg            <- length(sg_results)
    show_bact_panel <- !node_names && nrow(bacteria_map_filtered) > 0

    legend_cm <- 2.5 * 2.54   # 2.5 inches → fixed legend column (lcm)
    pdf_width <- width

    if (!is.null(dev.list())) dev.off()

    sg_tag   <- if (!is.null(subgroup)) sanitize_tag(subgroup) else "all"
    first_sg <- sg_results[[which(vapply(sg_results, function(r) igraph::gsize(r$network) > 0, logical(1)))[1]]]

    pdf(
      file.path(
        out_pdf,
        sprintf(
          "network.%s.%s.(%s).(%s).%s.%s.%s%s%s%s%s.pdf",
          mainGroup, sg_tag,
          analysis_method, cutoff,
          first_sg$sig,
          ifelse(node_names, "Names", "IDs"),
          ifelse(is.null(tab1_name), "", paste0(sanitize_tag(tab1_name), ".")),
          ifelse(is.null(tab2_name), "", paste0(sanitize_tag(tab2_name), ".")),
          ifelse(is.null(tab3_name), "", paste0(sanitize_tag(tab3_name), ".")),
          ifelse(is.null(name),      "", paste0(sanitize_tag(name), ".")),
          format(Sys.Date(), "%y%m%d")
        )
      ),
      width = pdf_width, height = height
    )

    if (show_bact_panel) {
      # Top row: network(s) + fixed legend column; bottom row: bacteria full-width
      legend_id   <- n_sg + 1L
      bacteria_id <- n_sg + 2L
      mat <- rbind(
        c(seq_len(n_sg), legend_id),
        c(rep(bacteria_id, n_sg), bacteria_id)
      )
      widths  <- c(rep(1, n_sg), lcm(legend_cm))
      heights <- c(height * 0.70, height * 0.30)
    } else {
      legend_id <- n_sg + 1L
      mat <- rbind(
        c(seq_len(n_sg), legend_id),
        c(rep(0, n_sg), 0)
      )
      widths  <- c(rep(1, n_sg), lcm(legend_cm))
      heights <- c(height * 0.85, height * 0.15)
    }
    graphics::layout(mat, widths = widths, heights = heights)

    # ── Panels: per-subgroup networks ────────────────────────────────────
    for (k in seq_along(sg_results)) {
      res <- sg_results[[k]]
      vis <- sg_visuals[[k]]
      net <- res$network

      graphics::par(mar = c(0.3, 0.3, 1.5, 0.3))

      if (igraph::gsize(net) == 0) {
        graphics::plot.new()
        graphics::title(main = sprintf("%s — %s\n(no edges)", mainGroup, res$sg), cex.main = 0.9)
        next
      }

      x_rng <- range(vis$norm_lay[, 1], na.rm = TRUE)
      y_rng <- range(vis$norm_lay[, 2], na.rm = TRUE)
      x_pad <- max(0.06, diff(x_rng) * 0.08)
      y_pad <- max(0.06, diff(y_rng) * 0.08)

      plot(
        net,
        layout             = vis$norm_lay,
        rescale            = FALSE,
        xlim               = c(x_rng[1] - x_pad, x_rng[2] + x_pad),
        ylim               = c(y_rng[1] - y_pad, y_rng[2] + y_pad),
        vertex.color       = vis$node_colors,
        vertex.label       = vis$node_labels,
        vertex.label.cex   = node_font,
        vertex.label.font  = vis$label_fonts,
        vertex.frame.color = vis$frame_colors,
        vertex.frame.width = ifelse(vis$frame_colors == "white", 1, 2.5),
        vertex.size        = vis$node_sizes,
        edge.width         = abs(igraph::E(net)$Correlation) * 3,
        edge.color         = ifelse(igraph::E(net)$Correlation > 0, "#C0392B", "#1A5276"),
        edge.lty           = vis$edge_styles,
        edge.curved        = 0.15,
        main               = sprintf("%s%s\n(%s; %s)",
                                      mainGroup,
                                      if (res$sg == "all") "" else sprintf(" — %s", res$sg),
                                      analysis_method, res$sigval)
      )

      if (vis$has_target) {
        for (i in which(vis$target_idx)) {
          nx    <- vis$norm_lay[i, 1]
          ny    <- vis$norm_lay[i, 2]
          y_off <- if (ny < 0) 0.14 else -0.14
          graphics::text(
            x      = nx, y = ny + y_off,
            labels = igraph::V(net)$name[i],
            cex    = node_font * 0.85,
            font   = 4,
            col    = "#E67E22",
            adj    = c(0.5, if (ny < 0) 0 else 1)
          )
        }
      }
    }

    # ── Panel: shared legend (right top, panel 2) ────────────────────────
    first_valid_k <- which(vapply(sg_results, function(r) igraph::gsize(r$network) > 0, logical(1)))[1]
    vis_leg       <- sg_visuals[[first_valid_k]]
    res_leg       <- sg_results[[first_valid_k]]

    graphics::par(mar = c(0.05, 0.05, 0.05, 0.05))
    graphics::plot.new()

    # Node types — top
    graphics::legend(
      x = 0.05, y = 0.95,
      legend = vis_leg$color_info$legend_labels,
      col    = vis_leg$color_info$legend_colors,
      pch    = 19,
      pt.cex = 1.35,
      bty    = "n",
      cex    = 0.92,
      title  = vis_leg$color_info$legend_title,
      xjust  = 0, yjust = 1,
      horiz  = FALSE
    )

    # Edge sign — middle
    graphics::legend(
      x = 0.05, y = 0.60,
      legend = c("Positive association", "Negative association"),
      col    = c("#C0392B", "#1A5276"),
      lty    = 1,
      lwd    = 3,
      bty    = "n",
      cex    = 0.92,
      title  = "Edge sign",
      xjust  = 0, yjust = 1,
      horiz  = FALSE
    )

    # FDR significance — bottom
    if (res_leg$sig == "p_FDR") {
      graphics::legend(
        x = 0.05, y = 0.30,
        legend = c(sprintf("%s < 0.05",  res_leg$signame),
                   sprintf("%s >= 0.05", res_leg$signame)),
        lty    = c(1, 2),
        lwd    = 3,
        bty    = "n",
        cex    = 0.92,
        title  = sprintf("%s Significance", res_leg$signame),
        xjust  = 0, yjust = 1,
        horiz  = FALSE
      )
    }

    # ── Panel: shared bacteria list (bottom full-width, panel 3) ─────────
    if (show_bact_panel) {
      graphics::par(mar = c(0.05, 0.2, 0.05, 0.2))
      graphics::plot.new()

      bacteria_legend_colors <- vapply(bacteria_map_filtered$Bacteria, function(b) {
        for (k in seq_along(sg_results)) {
          v_names <- igraph::V(sg_results[[k]]$network)$name
          idx     <- match(b, v_names)
          if (!is.na(idx) && !is.null(sg_visuals[[k]]$node_colors)) {
            return(sg_visuals[[k]]$node_colors[idx])
          }
        }
        "#888888"
      }, character(1))

      graphics::legend(
        x = 0.02, y = 0.98,
        legend    = paste(bacteria_map_filtered$ID, bacteria_map_filtered$Bacteria, sep = " = "),
        cex       = 0.68,
        bty       = "n",
        text.font = 3,
        text.col  = bacteria_legend_colors,
        title     = "Node ID",
        ncol      = max(2, min(4, ceiling(nrow(bacteria_map_filtered) / 20))),
        xjust     = 0, yjust = 1
      )
    }

    graphics::layout(1)
    dev.off()
  }

  # ── Centrality (all subgroups combined) ────────────────────────────────────
  all_centrality <- dplyr::bind_rows(lapply(sg_results, function(res) {
    net <- res$network
    if (igraph::gsize(net) == 0) return(NULL)
    cd <- data.frame(
      Subgroup    = res$sg,
      Node        = igraph::V(net)$name,
      Degree      = igraph::degree(net),
      Betweenness = igraph::betweenness(net),
      Closeness   = igraph::closeness(net),
      Eigenvector = igraph::eigen_centrality(net)$vector,
      stringsAsFactors = FALSE
    )
    cd[order(-cd$Degree), , drop = FALSE]
  }))

  first_sig <- sg_results[[1]]$sig
  utils::write.csv(all_centrality,
                   file.path(out_network, sprintf("%s.%s.sorted_centrality.%s.csv",
                                                  mainGroup, first_sig, format(Sys.Date(), "%y%m%d"))),
                   row.names = FALSE)

  invisible(sg_results)
}
