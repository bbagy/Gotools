#' Sample clustering with ordination and optional PERMANOVA
#'
#' @param psIN        Phyloseq object or numeric matrix/data.frame (samples = rows).
#' @param project     Project name string.
#' @param method      \code{"DMM"}, \code{"PAM"}, \code{"kmeans"}, \code{"hierarchical"}.
#'   Ignored when \code{k} is supplied.
#' @param transform   \code{"auto"}, \code{"compositional"}, \code{"clr"},
#'   \code{"log"}, \code{"none"}.
#' @param k           \code{NULL} = auto-select; integer = force k (PAM used).
#' @param max_k       Upper bound for k search.
#' @param distance_metric  Distance for ordination/PERMANOVA.
#' @param plot        \code{"PCoA"} or \code{"PCA"}.
#' @param mycols      Color vector for clusters.
#' @param ellipse     \code{TRUE} = 95\% ellipse; \code{FALSE} = none;
#'   string = group column name for ellipse.
#' @param marginal    Add marginal density panels via ggExtra.
#' @param marginal_type  \code{"density"}, \code{"boxplot"}, \code{"histogram"}.
#' @param marginal_size  Relative size of marginal panels.
#' @param marginal_alpha Transparency of marginal fills.
#' @param statistics  Run PERMANOVA and annotate plot.
#' @param addnumber   Append \code{(n=N)} to cluster legend labels.
#' @param sample_labels \code{FALSE} = none; \code{TRUE} = rownames;
#'   string = metadata column to use as labels.
#' @param name        Cluster column name suffix in output files.
#' @param height      PDF height in inches.
#' @param width       PDF width in inches.
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
    sample_labels   = FALSE,
    name            = NULL,
    height          = 4,
    width           = 4
) {
  marginal_type <- match.arg(marginal_type)
  cluster_col   <- if (!is.null(name)) name else "Cluster"
  date_sfx      <- format(Sys.Date(), "%y%m%d")
  name_sfx      <- if (!is.null(name)) paste0(".", name) else ""
  out_base      <- sprintf("%s_%s", project, date_sfx)
  out_pdf       <- file.path(out_base, "pdf")
  out_clust     <- file.path(out_base, "table", "cluster")
  for (d in c(out_pdf, out_clust)) if (!dir.exists(d)) dir.create(d, recursive = TRUE)

  # ── 1. input ──────────────────────────────────────────────────────
  is_ps <- inherits(psIN, "phyloseq")
  if (is_ps) {
    otu_raw <- as(otu_table(psIN), "matrix")
    if (taxa_are_rows(psIN)) otu_raw <- t(otu_raw)
    meta_df <- data.frame(sample_data(psIN))
  } else {
    otu_raw <- as.matrix(psIN)
    if (nrow(otu_raw) < ncol(otu_raw)) otu_raw <- t(otu_raw)
    meta_df <- NULL
  }

  # ── 2. transform ──────────────────────────────────────────────────
  if (transform == "auto") transform <- if (is_ps) "compositional" else "log"
  otu_trans <- switch(transform,
    compositional = sweep(otu_raw, 1, rowSums(otu_raw), "/"),
    clr  = { p <- otu_raw + 1e-6; log(p) - rowMeans(log(p)) },
    log  = log1p(otu_raw),
    none = otu_raw,
    stop(sprintf("[Go_cluster] Unknown transform: '%s'", transform))
  )

  # ── 3. clustering ─────────────────────────────────────────────────
  if (!is.null(k) && is.numeric(k)) {
    best_k <- as.integer(k)
    clust_labels <- cluster::pam(otu_trans, k = best_k)$clustering
    method_used  <- sprintf("PAM (k=%d forced)", best_k)
    message(sprintf("[Go_cluster] k=%d forced → PAM", best_k))
  } else {
    method <- match.arg(method, c("DMM", "PAM", "kmeans", "hierarchical"))
    message(sprintf("[Go_cluster] Method: %s", method))

    if (method == "DMM") {
      if (!requireNamespace("DirichletMultinomial", quietly = TRUE))
        stop("[Go_cluster] Install DirichletMultinomial via BiocManager.")
      count_mat <- matrix(as.integer(round(otu_trans * 1000)),
                          nrow = nrow(otu_trans), ncol = ncol(otu_trans),
                          dimnames = dimnames(otu_trans))
      fit_list  <- lapply(seq_len(max_k), function(ki)
        DirichletMultinomial::dmn(count = count_mat, k = ki))
      best_k    <- which.min(sapply(fit_list, DirichletMultinomial::laplace))
      clust_labels <- apply(DirichletMultinomial::mixture(fit_list[[best_k]]), 1, which.max)

    } else if (method == "PAM") {
      sil    <- sapply(2:max_k, function(ki)
        mean(cluster::silhouette(cluster::pam(otu_trans, k = ki))[, 3]))
      best_k <- which.max(sil) + 1
      clust_labels <- cluster::pam(otu_trans, k = best_k)$clustering

    } else if (method == "kmeans") {
      d_sil  <- dist(otu_trans)
      sil    <- sapply(2:max_k, function(ki) {
        set.seed(42); km <- kmeans(otu_trans, centers = ki, nstart = 25)
        mean(cluster::silhouette(km$cluster, d_sil)[, 3])
      })
      best_k <- which.max(sil) + 1
      set.seed(42); clust_labels <- kmeans(otu_trans, centers = best_k, nstart = 25)$cluster

    } else {  # hierarchical
      d_hc   <- dist(otu_trans)
      hc     <- hclust(d_hc, method = "ward.D2")
      sil    <- sapply(2:max_k, function(ki)
        mean(cluster::silhouette(cutree(hc, k = ki), d_hc)[, 3]))
      best_k <- which.max(sil) + 1
      clust_labels <- cutree(hc, k = best_k)
    }
    method_used <- method
    message(sprintf("[Go_cluster] best k = %d", best_k))
  }

  # ── 4. attach cluster labels ──────────────────────────────────────
  names(clust_labels) <- rownames(otu_trans)
  if (is_ps) {
    meta_df[[cluster_col]] <- factor(clust_labels[rownames(meta_df)])
    sample_data(psIN)      <- sample_data(meta_df)
  } else {
    meta_df <- stats::setNames(data.frame(factor(clust_labels)), cluster_col)
  }

  # ── 5. distance ───────────────────────────────────────────────────
  dist_mat <- tryCatch(
    vegan::vegdist(otu_trans, method = distance_metric),
    error = function(e) { message("[Go_cluster] vegdist failed; using euclidean."); dist(otu_trans) }
  )

  # ── 6. ordination ─────────────────────────────────────────────────
  if (plot == "PCoA") {
    pcoa_res <- ape::pcoa(dist_mat)
    coord_df <- as.data.frame(pcoa_res$vectors[, 1:2, drop = FALSE])
    colnames(coord_df) <- c("Axis_1", "Axis_2")
    eig  <- pcoa_res$values$Eigenvalues
    den  <- sum(eig[eig > 0], na.rm = TRUE)
    pct  <- pmax(eig[1:2], 0) / den * 100
    xlab <- sprintf("PCoA 1 (%.1f%%)", pct[1])
    ylab <- sprintf("PCoA 2 (%.1f%%)", pct[2])
  } else {  # PCA
    nz <- apply(otu_trans, 2, var) > 0
    if (any(!nz)) {
      message(sprintf("[Go_cluster] PCA: removing %d zero-variance features.", sum(!nz)))
      otu_trans <- otu_trans[, nz, drop = FALSE]
    }
    pca_res  <- prcomp(otu_trans, scale. = TRUE)
    coord_df <- as.data.frame(pca_res$x[, 1:2, drop = FALSE])
    colnames(coord_df) <- c("Axis_1", "Axis_2")
    imp  <- summary(pca_res)$importance
    xlab <- sprintf("PC1 (%.1f%%)", imp[2, 1] * 100)
    ylab <- sprintf("PC2 (%.1f%%)", imp[2, 2] * 100)
  }

  # ── 7. plot dataframe ─────────────────────────────────────────────
  plot_df              <- cbind(coord_df, meta_df)
  plot_df[[cluster_col]] <- factor(plot_df[[cluster_col]])

  utils::write.csv(plot_df, file.path(out_clust, sprintf(
    "mapping.%s.%s%s.%s.csv", method_used, project, name_sfx, date_sfx
  )), row.names = TRUE)

  if (isTRUE(addnumber)) {
    for (lv in levels(plot_df[[cluster_col]])) {
      n <- sum(plot_df[[cluster_col]] == lv, na.rm = TRUE)
      levels(plot_df[[cluster_col]])[levels(plot_df[[cluster_col]]) == lv] <-
        sprintf("%s (n=%d)", lv, n)
    }
  }

  if (!isFALSE(sample_labels))
    plot_df$.label_ <- if (isTRUE(sample_labels)) rownames(plot_df) else plot_df[[sample_labels]]

  # ── 8. plot ───────────────────────────────────────────────────────
  p <- ggplot2::ggplot(plot_df,
      ggplot2::aes(x = Axis_1, y = Axis_2,
                   color = .data[[cluster_col]])) +
    ggplot2::geom_point(size = 2, alpha = 0.8) +
    ggplot2::geom_vline(xintercept = 0, linewidth = 0.1) +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.1) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = 0.03)) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.06, 0.03))) +
    ggplot2::labs(x = xlab, y = ylab,
      title = sprintf("Clustering  |  %s  |  k=%d", method_used, best_k)) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      panel.grid        = ggplot2::element_blank(),
      legend.position   = "bottom",
      legend.title      = ggplot2::element_blank(),
      legend.key        = ggplot2::element_blank(),
      panel.background  = ggplot2::element_rect(fill = "white", colour = "black", linewidth = 0.7),
      aspect.ratio      = 1,
      plot.title        = ggplot2::element_text(size = 8, face = "bold")
    )

  if (!is.null(mycols))
    p <- p + ggplot2::scale_color_manual(values = mycols)

  if (isTRUE(ellipse)) {
    p <- p + ggplot2::stat_ellipse(type = "norm", linetype = 2)
  } else if (is.character(ellipse)) {
    p <- p + ggplot2::stat_ellipse(
      ggplot2::aes(group = .data[[ellipse]], color = .data[[ellipse]]),
      type = "norm", linetype = 2)
  }

  if (!isFALSE(sample_labels))
    p <- p + ggrepel::geom_text_repel(
      ggplot2::aes(label = .label_),
      size = 2.5, max.overlaps = 30, segment.size = 0.3, show.legend = FALSE)

  # ── 9. PERMANOVA ─────────────────────────────────────────────────
  if (isTRUE(statistics)) {
    perm_df <- stats::setNames(
      data.frame(factor(clust_labels[rownames(otu_trans)])), cluster_col)

    if (nrow(perm_df) >= 3 && length(unique(perm_df[[cluster_col]])) >= 2) {
      set.seed(1)
      ad <- vegan::adonis2(
        stats::as.formula(sprintf("dist_mat ~ %s", cluster_col)),
        data = perm_df, permutations = 999, by = "terms")
      R2 <- round(ad[1, "R2"], 3); pval <- ad[1, "Pr(>F)"]

      utils::write.csv(as.data.frame(ad), file.path(out_clust, sprintf(
        "PERMANOVA.%s.%s%s.%s.csv", method_used, project, name_sfx, date_sfx
      )), row.names = TRUE)

      p <- p + ggplot2::annotate("text", x = -Inf, y = -Inf,
        label = sprintf("%s\nR2=%.3f\np=%.3f", distance_metric, R2, pval),
        size = 3, hjust = -0.05, vjust = -0.3, lineheight = 0.95)

      # pairwise
      dm   <- as.matrix(dist_mat)
      lvls <- levels(perm_df[[cluster_col]])
      pw   <- do.call(rbind, Filter(Negate(is.null), lapply(
        utils::combn(lvls, 2, simplify = FALSE), function(pr) {
          idx <- perm_df[[cluster_col]] %in% pr
          sd  <- stats::as.dist(dm[idx, idx])
          sl  <- droplevels(perm_df[idx, , drop = FALSE])
          if (nrow(sl) < 3 || length(unique(sl[[cluster_col]])) < 2) return(NULL)
          set.seed(1)
          r <- vegan::adonis2(stats::as.formula(sprintf("sd ~ %s", cluster_col)),
                              data = sl, permutations = 999, by = "terms")
          data.frame(Group1=pr[1], Group2=pr[2],
                     n1=sum(perm_df[[cluster_col]]==pr[1]),
                     n2=sum(perm_df[[cluster_col]]==pr[2]),
                     R2=round(r[1,"R2"],3), F_stat=round(r[1,"F"],3),
                     p_value=r[1,"Pr(>F)"])
        })))
      if (!is.null(pw) && nrow(pw) > 0) {
        pw$p_adj <- stats::p.adjust(pw$p_value, method = "BH")
        pw$sig   <- dplyr::case_when(pw$p_adj<0.001~"***",pw$p_adj<0.01~"**",pw$p_adj<0.05~"*",TRUE~"ns")
        utils::write.csv(pw, file.path(out_clust, sprintf(
          "PERMANOVA_pairwise.%s.%s%s.%s.csv", method_used, project, name_sfx, date_sfx
        )), row.names = FALSE)
        print(pw)
      }
    }
  }

  # ── 10. marginal ──────────────────────────────────────────────────
  if (isTRUE(marginal) && requireNamespace("ggExtra", quietly = TRUE))
    p <- ggExtra::ggMarginal(p, type = marginal_type, size = marginal_size,
                              alpha = marginal_alpha,
                              groupColour = TRUE, groupFill = TRUE)

  # ── 11. save ──────────────────────────────────────────────────────
  fn_pdf <- file.path(out_pdf, sprintf("cluster.%s.%s%s.%s.pdf",
                                        method_used, project, name_sfx, date_sfx))
  if (inherits(p, "ggExtraPlot")) {
    grDevices::pdf(fn_pdf, height = height, width = width); print(p); grDevices::dev.off()
  } else {
    ggplot2::ggsave(fn_pdf, plot = p, device = grDevices::cairo_pdf,
                    height = height, width = width)
  }
  message(sprintf("[Go_cluster] done → %s  (k=%d)", fn_pdf, best_k))

  if (is_ps) invisible(psIN) else
    invisible(stats::setNames(data.frame(factor(clust_labels)), cluster_col))
}
