#' Sample clustering with ordination and PERMANOVA
#'
#' @param psIN  Phyloseq object **or** numeric matrix/data.frame (samples = rows).
#' @param project  Project name string.
#' @param method  Clustering method: \code{"DMM"}, \code{"PAM"}, \code{"kmeans"},
#'   or \code{"hierarchical"}. Ignored when \code{k} is a number.
#' @param transform  Count transformation: \code{"auto"} (compositional for phyloseq,
#'   log1p for matrix), \code{"compositional"}, \code{"clr"}, \code{"log"}, \code{"none"}.
#' @param k  \code{NULL} = auto-select optimal k. A positive integer = force this many
#'   clusters (method ignored; PAM used for assignment).
#' @param max_k  Upper bound for k search when \code{k = NULL}.
#' @param distance_metric  Distance for ordination and PERMANOVA (e.g. \code{"bray"}).
#' @param plot  Ordination type: \code{"PCoA"} or \code{"PCA"}.
#' @param mycols  Color vector for clusters.
#' @param ellipse  \code{TRUE} = 95\% normal ellipses; \code{FALSE} = none;
#'   column name = ellipse grouped by that variable.
#' @param marginal  Add marginal panels via ggExtra (\code{TRUE}/\code{FALSE}).
#' @param marginal_type  \code{"density"}, \code{"boxplot"}, or \code{"histogram"}.
#' @param marginal_size  Relative size of marginal panels.
#' @param marginal_alpha  Transparency of marginal fills.
#' @param statistics  Run global + pairwise PERMANOVA and annotate plot.
#' @param addnumber  Append \code{(n=N)} to cluster legend labels.
#' @param name  If supplied, used as the cluster column name in sample_data
#'   (default: \code{"Cluster"}) and as a suffix in output file names.
#' @param height  PDF height in inches.
#' @param width   PDF width in inches.
#'
#' @return Invisibly returns \code{psIN} with a cluster column added to sample_data
#'   (phyloseq input), or a data.frame of cluster assignments (matrix input).
#' @export
Go_cluster <- function(
    psIN,
    project,
    method          = "DMM",
    transform       = "auto",
    k               = NULL,
    max_k           = 6,
    distance_metric = "bray",
    plot            = "PCoA",
    mycols          = NULL,
    ellipse         = TRUE,
    marginal        = FALSE,
    marginal_type   = c("density", "boxplot", "histogram"),
    marginal_size   = 5,
    marginal_alpha  = 0.35,
    statistics      = TRUE,
    addnumber       = FALSE,
    name            = NULL,
    height          = 4,
    width           = 4
) {

  marginal_type <- match.arg(marginal_type)

  # cluster column name: use name if provided
  cluster_col <- if (!is.null(name)) name else "Cluster"

  # ── method descriptions ───────────────────────────────────────────
  .method_desc <- c(
    DMM          = "Dirichlet Multinomial Mixture: probabilistic model for count data; optimal k by Laplace criterion. Best for 16S / metagenomics.",
    PAM          = "Partitioning Around Medoids: robust to outliers; optimal k by mean silhouette width. Works for any data type.",
    kmeans       = "k-means: fast Euclidean clustering; optimal k by mean silhouette width. Good for normalized RNA-seq / continuous data.",
    hierarchical = "Hierarchical clustering (Ward.D2): dendrogram-based; k by silhouette after tree cutting. Useful for exploratory analysis."
  )

  # ── output directories ────────────────────────────────────────────
  date_sfx  <- format(Sys.Date(), "%y%m%d")
  name_sfx  <- if (!is.null(name)) paste0(".", name) else ""
  out_base  <- sprintf("%s_%s", project, date_sfx)
  out_pdf   <- file.path(out_base, "pdf")
  out_clust <- file.path(out_base, "table", "cluster")   # ← new
  for (d in c(out_pdf, out_clust))
    if (!dir.exists(d)) dir.create(d, recursive = TRUE)

  # ── 1. input type ─────────────────────────────────────────────────
  is_ps <- inherits(psIN, "phyloseq")

  if (is_ps) {
    otu_raw    <- as(otu_table(psIN), "matrix")
    if (taxa_are_rows(psIN)) otu_raw <- t(otu_raw)  # samples = rows
    meta_df    <- data.frame(sample_data(psIN))
    input_type <- "phyloseq"
  } else {
    otu_raw <- as.matrix(psIN)
    if (nrow(otu_raw) < ncol(otu_raw)) otu_raw <- t(otu_raw)
    meta_df    <- NULL
    input_type <- "matrix"
  }

  # ── 2. transformation ─────────────────────────────────────────────
  if (transform == "auto")
    transform <- if (is_ps) "compositional" else "log"

  otu_trans <- switch(transform,
    compositional = sweep(otu_raw, 1, rowSums(otu_raw), "/"),
    clr = {
      pos <- otu_raw + 1e-6
      log(pos) - rowMeans(log(pos))
    },
    log  = log1p(otu_raw),
    none = otu_raw,
    stop(sprintf("[Go_cluster] Unknown transform: '%s'", transform))
  )

  # ── 3. clustering ─────────────────────────────────────────────────
  forced_k <- !is.null(k) && is.numeric(k)

  if (forced_k) {
    best_k <- as.integer(k)
    message(sprintf(
      "\n[Go_cluster] k = %d forced → method ignored; PAM used for assignment.",
      best_k
    ))
    clust_labels <- cluster::pam(otu_trans, k = best_k)$clustering
    method_used  <- sprintf("PAM (k=%d forced)", best_k)

  } else {
    method <- match.arg(method, c("DMM", "PAM", "kmeans", "hierarchical"))
    message(sprintf("\n[Go_cluster] Method: %s\n  → %s",
                    method, .method_desc[method]))

    if (method == "DMM") {
      if (!requireNamespace("DirichletMultinomial", quietly = TRUE))
        stop("[Go_cluster] Install DirichletMultinomial via BiocManager.")
      # otu_trans (compositional 0~1) × 1000 → integer overflow 방지
      count_mat <- matrix(
        as.integer(round(otu_trans * 1000)),
        nrow = nrow(otu_trans), ncol = ncol(otu_trans),
        dimnames = dimnames(otu_trans)
      )
      fit_list  <- lapply(seq_len(max_k), function(ki)
        DirichletMultinomial::dmn(count = count_mat, k = ki))
      lap    <- sapply(fit_list, DirichletMultinomial::laplace)
      best_k <- which.min(lap)
      clust_labels <- apply(
        DirichletMultinomial::mixture(fit_list[[best_k]]), 1, which.max)
      message(sprintf(
        "[Go_cluster] DMM: tested k=1~%d → best k=%d (Laplace=%.2f)",
        max_k, best_k, min(lap)
      ))

    } else if (method == "PAM") {
      sil <- sapply(2:max_k, function(ki) {
        res <- cluster::pam(otu_trans, k = ki)
        mean(cluster::silhouette(res)[, 3])
      })
      best_k       <- which.max(sil) + 1
      clust_labels <- cluster::pam(otu_trans, k = best_k)$clustering
      message(sprintf(
        "[Go_cluster] PAM: tested k=2~%d → best k=%d (mean silhouette=%.3f)",
        max_k, best_k, max(sil)
      ))

    } else if (method == "kmeans") {
      d_sil <- dist(otu_trans)
      sil <- sapply(2:max_k, function(ki) {
        set.seed(42)
        km <- kmeans(otu_trans, centers = ki, nstart = 25)
        mean(cluster::silhouette(km$cluster, d_sil)[, 3])
      })
      best_k <- which.max(sil) + 1
      set.seed(42)
      clust_labels <- kmeans(otu_trans, centers = best_k, nstart = 25)$cluster
      message(sprintf(
        "[Go_cluster] k-means: tested k=2~%d → best k=%d (mean silhouette=%.3f)",
        max_k, best_k, max(sil)
      ))

    } else {  # hierarchical
      d_hc <- dist(otu_trans)
      hc   <- hclust(d_hc, method = "ward.D2")
      sil  <- sapply(2:max_k, function(ki)
        mean(cluster::silhouette(cutree(hc, k = ki), d_hc)[, 3]))
      best_k       <- which.max(sil) + 1
      clust_labels <- cutree(hc, k = best_k)
      message(sprintf(
        "[Go_cluster] Hierarchical (ward.D2): tested k=2~%d → best k=%d (mean silhouette=%.3f)",
        max_k, best_k, max(sil)
      ))
    }
    method_used <- method
  }

  # ── 4. attach cluster labels → sample_data ────────────────────────
  names(clust_labels) <- rownames(otu_trans)

  if (is_ps) {
    meta_df[[cluster_col]] <- factor(clust_labels[rownames(meta_df)])
    sample_data(psIN)      <- sample_data(meta_df)
  } else {
    meta_df <- stats::setNames(
      data.frame(factor(clust_labels)), cluster_col)
  }

  # ── 5. distance matrix ────────────────────────────────────────────
  dist_mat <- tryCatch(
    vegan::vegdist(otu_trans, method = distance_metric),
    error = function(e) {
      message(sprintf(
        "[Go_cluster] vegdist('%s') failed; falling back to euclidean.",
        distance_metric
      ))
      dist(otu_trans)
    }
  )

  # ── 6. ordination + variance explained ───────────────────────────
  .pcoa_pct <- function(eig_vec) {
    # robust: use only positive eigenvalues as denominator
    den <- sum(eig_vec[eig_vec > 0], na.rm = TRUE)
    if (!is.finite(den) || den <= 0)
      den <- sum(abs(eig_vec), na.rm = TRUE)
    pmax(eig_vec[1:2], 0) / den * 100
  }

  if (plot == "PCoA") {
    pcoa_res <- ape::pcoa(dist_mat)
    coord_df <- as.data.frame(pcoa_res$vectors[, 1:2, drop = FALSE])
    colnames(coord_df) <- c("Axis_1", "Axis_2")

    eig_vals <- pcoa_res$values$Eigenvalues
    pct      <- .pcoa_pct(eig_vals)
    xlab <- sprintf("PCoA 1 (%.1f%%)", pct[1])
    ylab <- sprintf("PCoA 2 (%.1f%%)", pct[2])

    # variance table (all axes)
    all_pct <- pmax(eig_vals, 0) /
               sum(eig_vals[eig_vals > 0], na.rm = TRUE) * 100
    variance_df <- data.frame(
      Axis               = paste0("PCoA", seq_along(all_pct)),
      Eigenvalue         = round(eig_vals, 4),
      Variance_Explained = round(all_pct, 2),
      Cumulative         = round(cumsum(all_pct), 2)
    )

  } else {  # PCA
    pca_res  <- prcomp(otu_trans, scale. = TRUE)
    coord_df <- as.data.frame(pca_res$x[, 1:2, drop = FALSE])
    colnames(coord_df) <- c("Axis_1", "Axis_2")
    imp  <- summary(pca_res)$importance
    pct  <- imp[2, 1:2] * 100
    xlab <- sprintf("PC1 (%.1f%%)", imp[2, 1] * 100)
    ylab <- sprintf("PC2 (%.1f%%)", imp[2, 2] * 100)

    variance_df <- data.frame(
      Axis               = paste0("PC", seq_len(ncol(pca_res$x))),
      Variance_Explained = round(imp[2, ] * 100, 2),
      Cumulative         = round(imp[3, ] * 100, 2)
    )
  }

  # save variance table
  fn_var <- file.path(out_clust, sprintf(
    "variance_explained.%s.%s%s.%s.csv",
    method_used, project, name_sfx, date_sfx
  ))
  utils::write.csv(variance_df, fn_var, row.names = FALSE)
  message(sprintf(
    "[Go_cluster] %s 1: %.1f%%  |  %s 2: %.1f%%  (variance explained)",
    plot, pct[1], plot, pct[2]
  ))

  # ── 7. build plot dataframe ───────────────────────────────────────
  plot_df              <- cbind(coord_df, meta_df)
  plot_df[[cluster_col]] <- factor(plot_df[[cluster_col]])

  # mapping CSV는 addnumber 적용 전 원본 레이블로 저장
  fn_map <- file.path(out_clust, sprintf(
    "mapping.%s.%s%s.%s.csv",
    method_used, project, name_sfx, date_sfx
  ))
  utils::write.csv(plot_df, fn_map, row.names = TRUE)

  # addnumber는 plot legend에만 적용
  if (isTRUE(addnumber)) {
    for (lv in levels(plot_df[[cluster_col]])) {
      n <- sum(plot_df[[cluster_col]] == lv, na.rm = TRUE)
      levels(plot_df[[cluster_col]])[
        levels(plot_df[[cluster_col]]) == lv] <-
        sprintf("%s (n=%d)", lv, n)
    }
  }

  # ── 8. ggplot ─────────────────────────────────────────────────────
  p <- ggplot2::ggplot(
      plot_df,
      ggplot2::aes_string(x = "Axis_1", y = "Axis_2", color = cluster_col)
    ) +
    ggplot2::geom_point(size = 2, alpha = 0.8) +
    ggplot2::labs(
      x     = xlab,
      y     = ylab,
      title = sprintf("Sample clustering  |  method=%s  |  k=%d",
                      method_used, best_k)
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      strip.background     = ggplot2::element_blank(),
      legend.position      = "bottom",
      legend.title         = ggplot2::element_blank(),
      legend.justification = "left",
      legend.box           = "vertical",
      legend.box.margin    = ggplot2::margin(0, 0, 0, -1, "cm"),
      legend.spacing.y     = ggplot2::unit(0.02, "cm"),
      legend.key.height    = ggplot2::unit(0.25, "cm"),
      legend.key.width     = ggplot2::unit(0.35, "cm"),
      panel.grid           = ggplot2::element_blank(),
      legend.key           = ggplot2::element_blank(),
      panel.background     = ggplot2::element_rect(
        fill = "white", colour = "black", linewidth = 0.7),
      aspect.ratio         = 1,
      plot.title           = ggplot2::element_text(size = 8, face = "bold")
    ) +
    ggplot2::geom_vline(xintercept = 0, linewidth = 0.1) +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.1) +
    ggplot2::scale_x_continuous(
      expand = ggplot2::expansion(mult = c(0.03, 0.03))) +
    ggplot2::scale_y_continuous(
      expand = ggplot2::expansion(mult = c(0.06, 0.03)))

  if (!is.null(mycols))
    p <- p + ggplot2::scale_color_manual(values = mycols)

  if (isTRUE(ellipse)) {
    p <- p + ggplot2::stat_ellipse(type = "norm", linetype = 2)
  } else if (!isFALSE(ellipse) && !is.null(ellipse)) {
    p <- p + ggplot2::stat_ellipse(
      ggplot2::aes_string(group = ellipse, color = ellipse),
      type = "norm", linetype = 2
    )
  }

  # ── 9. PERMANOVA (global + pairwise) ─────────────────────────────
  if (isTRUE(statistics)) {
    perm_df <- data.frame(
      row.names = rownames(otu_trans),
      Cluster   = factor(clust_labels[rownames(otu_trans)])
    )
    names(perm_df) <- cluster_col

    n_samp   <- nrow(perm_df)
    n_groups <- length(unique(perm_df[[cluster_col]]))

    if (n_samp >= 3 && n_groups >= 2) {

      # --- global PERMANOVA ---
      set.seed(1)
      form_g <- stats::as.formula(sprintf("dist_mat ~ %s", cluster_col))
      ad_g   <- vegan::adonis2(form_g, data = perm_df,
                               permutations = 999, by = "terms")
      R2_g   <- round(ad_g[1, "R2"], 3)
      p_g    <- ad_g[1, "Pr(>F)"]

      fn_global <- file.path(out_clust, sprintf(
        "PERMANOVA_global.%s.%s%s.%s.csv",
        method_used, project, name_sfx, date_sfx
      ))
      utils::write.csv(as.data.frame(ad_g), fn_global, row.names = TRUE)

      message(sprintf(
        "[Go_cluster] Global PERMANOVA: R2=%.3f  p=%.3f  (%s)",
        R2_g, p_g, distance_metric
      ))

      # annotate plot
      ann <- data.frame(
        x = -Inf, y = -Inf,
        label = sprintf("%s\nR2=%.3f\nPERMANOVA p=%.3f",
                        distance_metric, R2_g, p_g)
      )
      p <- p + ggplot2::geom_text(
        data = ann,
        ggplot2::aes(x = x, y = y, label = label),
        size = 3, hjust = -0.005, vjust = -0.3,
        lineheight = 0.95, inherit.aes = FALSE
      )

      # --- pairwise PERMANOVA ---
      dist_mat_full <- as.matrix(dist_mat)
      lvls  <- levels(perm_df[[cluster_col]])
      pairs <- utils::combn(lvls, 2, simplify = FALSE)

      pair_rows <- lapply(pairs, function(pr) {
        idx <- perm_df[[cluster_col]] %in% pr
        sub_dist   <- stats::as.dist(dist_mat_full[idx, idx])
        sub_labels <- perm_df[idx, , drop = FALSE]
        sub_labels[[cluster_col]] <- droplevels(sub_labels[[cluster_col]])

        if (nrow(sub_labels) < 3 ||
            length(unique(sub_labels[[cluster_col]])) < 2)
          return(NULL)

        set.seed(1)
        form_p <- stats::as.formula(sprintf("sub_dist ~ %s", cluster_col))
        ad_p   <- vegan::adonis2(form_p, data = sub_labels,
                                 permutations = 999, by = "terms")
        data.frame(
          Group1  = pr[1],
          Group2  = pr[2],
          n1      = sum(perm_df[[cluster_col]] == pr[1]),
          n2      = sum(perm_df[[cluster_col]] == pr[2]),
          R2      = round(ad_p[1, "R2"],   3),
          F_stat  = round(ad_p[1, "F"],    3),
          p_value = ad_p[1, "Pr(>F)"],
          stringsAsFactors = FALSE
        )
      })

      pair_df <- do.call(rbind, Filter(Negate(is.null), pair_rows))
      if (!is.null(pair_df) && nrow(pair_df) > 0) {
        pair_df$p_adj <- stats::p.adjust(pair_df$p_value, method = "BH")
        pair_df$sig   <- dplyr::case_when(
          pair_df$p_adj < 0.001 ~ "***",
          pair_df$p_adj < 0.01  ~ "**",
          pair_df$p_adj < 0.05  ~ "*",
          TRUE                  ~ "ns"
        )

        fn_pair <- file.path(out_clust, sprintf(
          "PERMANOVA_pairwise.%s.%s%s.%s.csv",
          method_used, project, name_sfx, date_sfx
        ))
        utils::write.csv(pair_df, fn_pair, row.names = FALSE)
        message(sprintf(
          "[Go_cluster] Pairwise PERMANOVA: %d pairs  →  %s",
          nrow(pair_df), fn_pair
        ))
        print(pair_df)
      }
    }
  }

  # ── 10. marginal ─────────────────────────────────────────────────
  if (isTRUE(marginal)) {
    if (!requireNamespace("ggExtra", quietly = TRUE)) {
      message("[Go_cluster] ggExtra not found; marginal skipped.")
    } else {
      p <- ggExtra::ggMarginal(
        p,
        type        = marginal_type,
        size        = marginal_size,
        alpha       = marginal_alpha,
        groupColour = TRUE,
        groupFill   = TRUE
      )
    }
  }

  # ── 11. save PDF ─────────────────────────────────────────────────
  fn_pdf <- file.path(out_pdf, sprintf(
    "cluster.%s.%s%s.%s.pdf",
    method_used, project, name_sfx, date_sfx
  ))
  if (inherits(p, "ggExtraPlot")) {
    grDevices::pdf(fn_pdf, height = height, width = width)
    print(p)
    grDevices::dev.off()
  } else {
    ggplot2::ggsave(fn_pdf, plot = p, device = grDevices::cairo_pdf,
                   height = height, width = width)
  }

  # ── 12. summary message ──────────────────────────────────────────
  message(sprintf(paste0(
    "\n============================================\n",
    "  Go_cluster  —  done\n",
    "--------------------------------------------\n",
    "  Input        : %s\n",
    "  Method       : %s\n",
    "  Transform    : %s\n",
    "  Distance     : %s\n",
    "  Clusters     : k = %d\n",
    "  Column name  : %s\n",
    "  PDF          : %s\n",
    "  table/cluster: %s\n",
    "============================================\n"
  ), input_type, method_used, transform, distance_metric,
     best_k, cluster_col, fn_pdf, out_clust))

  if (is_ps) invisible(psIN) else
    invisible(stats::setNames(
      data.frame(factor(clust_labels)), cluster_col))
}
