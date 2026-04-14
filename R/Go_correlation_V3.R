
#' Correlation Analysis Between Microbiome Taxa and External Data
#'
#' Calculates pairwise Spearman (or Pearson/Kendall) correlations between
#' aggregated microbiome taxa and an external data table (e.g., metabolomics,
#' clinical variables). Taxa are filtered by prevalence and top N abundance
#' before correlation. FDR correction is applied across all tests at once.
#'
#' @param psIN Phyloseq object containing microbiome data.
#' @param project Project name for output file naming.
#' @param rank Taxonomic rank to aggregate taxa (e.g., "Genus", "Species").
#' @param Table Data frame of external variables (rows = samples, cols = variables).
#' @param method Correlation method: "spearman" (default), "pearson", or "kendall".
#' @param top Number of top taxa by mean abundance to include. Default 50.
#' @param min_prev Minimum prevalence (fraction of samples) for taxa to be included. Default 0.1.
#' @param pvalue Significance threshold for display. Default 0.05.
#' @param padj Logical. Apply BH-FDR correction across all tests. Default TRUE.
#' @param group Optional column name in sample_data for stratified analysis. If provided, correlation is calculated separately per group level.
#' @param orders Optional vector to order group levels.
#' @param xangle Angle for x-axis labels. Default 90.
#' @param ncols Number of columns for facet wrap (group panels). Default NULL (auto).
#' @param name Optional suffix for output file name.
#' @param height Height of the output PDF.
#' @param width Width of the output PDF.
#'
#' @return Saves a PDF heatmap and a CSV of correlation results. Returns the result data frame invisibly.
#'
#' @details
#' Taxa are aggregated to the specified rank, filtered by prevalence (min_prev),
#' then the top N by mean relative abundance are selected. Log(1+x) transformation
#' is applied if data are not already on a relative scale. FDR (BH) correction is
#' applied across all taxa x variable pairs (and groups if stratified).
#'
#' @examples
#' Go_correlation(psIN = ps, project = "MyProject", rank = "Genus",
#'                Table = metabolite_df, top = 30, min_prev = 0.1,
#'                method = "spearman", padj = TRUE, height = 10, width = 8)
#'
#' @param patchwork Logical. If \code{TRUE}, skip saving and return the plot object(s) for use with \code{Gg_patchwork()} or the \pkg{patchwork} package. Default \code{FALSE}.
#' @export

Go_correlation <- function(psIN, project,
                           rank,
                           Table,
                           method   = "spearman",
                           top      = 50,
                           min_prev = 0.1,
                           pvalue   = 0.05,
                           padj     = TRUE,
                           group    = NULL,
                           orders   = NULL,
                           xangle   = 90,
                           ncols    = NULL,
                           name     = NULL,
                           height, width,
                           patchwork = FALSE) {

  if (!is.null(dev.list())) dev.off()

  # output dirs
  out      <- file.path(sprintf("%s_%s", project, format(Sys.Date(), "%y%m%d")))
  out_path <- file.path(sprintf("%s_%s/pdf",   project, format(Sys.Date(), "%y%m%d")))
  out_tab  <- file.path(sprintf("%s_%s/table", project, format(Sys.Date(), "%y%m%d")))
  for (d in c(out, out_path, out_tab)) if (!dir.exists(d)) dir.create(d, recursive = TRUE)

  # ── 1. aggregate to rank ──────────────────────────────────────────────────
  ps_glom <- tryCatch(
    tax_glom(psIN, taxrank = rank, NArm = TRUE),
    error = function(e) {
      stop(sprintf("tax_glom failed for rank '%s': %s", rank, e$message))
    }
  )

  otu_mat <- as.data.frame(t(otu_table(ps_glom)))   # rows = samples, cols = taxa
  tax_mat <- as.data.frame(tax_table(ps_glom))

  # label columns by rank name
  colnames(otu_mat) <- tax_mat[[rank]]

  # ── 2. prevalence filter ─────────────────────────────────────────────────
  prev <- colMeans(otu_mat > 0)
  otu_mat <- otu_mat[, prev >= min_prev, drop = FALSE]
  if (ncol(otu_mat) == 0) stop("No taxa passed the prevalence filter. Lower min_prev.")
  message(sprintf("[Go_correlation] %d taxa passed prevalence filter (>= %.0f%%)",
                  ncol(otu_mat), min_prev * 100))

  # ── 3. top N by mean abundance ───────────────────────────────────────────
  mean_abund <- colMeans(otu_mat)
  top_n <- min(top, ncol(otu_mat))
  otu_mat <- otu_mat[, order(mean_abund, decreasing = TRUE)[seq_len(top_n)], drop = FALSE]
  message(sprintf("[Go_correlation] Using top %d taxa by mean abundance.", top_n))

  # ── 4. log transform if count data ───────────────────────────────────────
  if (max(otu_mat, na.rm = TRUE) > 1) {
    otu_mat <- log1p(otu_mat)
  }

  # ── 5. align samples between psIN and Table ───────────────────────────────
  Table <- as.data.frame(Table)
  Table[] <- lapply(Table, function(x) suppressWarnings(as.numeric(as.character(x))))
  common_samps <- intersect(rownames(otu_mat), rownames(Table))
  if (length(common_samps) < 5) stop("Fewer than 5 common samples between psIN and Table.")
  otu_mat <- otu_mat[common_samps, , drop = FALSE]
  Table   <- Table[common_samps, , drop = FALSE]

  # remove external variables with too many NAs (>50%)
  ok_vars <- colMeans(is.na(Table)) < 0.5
  Table   <- Table[, ok_vars, drop = FALSE]
  if (ncol(Table) == 0) stop("No external variables with sufficient data.")

  # ── 6. define groups ─────────────────────────────────────────────────────
  map <- as.data.frame(sample_data(psIN))
  map <- map[common_samps, , drop = FALSE]

  if (!is.null(group) && group %in% colnames(map)) {
    grp_vec <- as.character(map[[group]])
    grp_levels <- if (!is.null(orders)) intersect(orders, unique(grp_vec)) else sort(unique(grp_vec))
  } else {
    grp_vec    <- rep("All samples", nrow(otu_mat))
    grp_levels <- "All samples"
  }

  # ── 7. correlation ────────────────────────────────────────────────────────
  results <- data.frame()
  for (grp in grp_levels) {
    idx <- which(grp_vec == grp)
    if (length(idx) < 5) {
      message(sprintf("[Go_correlation] Skipping group '%s': < 5 samples.", grp))
      next
    }
    om  <- otu_mat[idx, , drop = FALSE]
    tab <- Table[idx,   , drop = FALSE]

    for (taxa_j in colnames(om)) {
      for (var_i in colnames(tab)) {
        x  <- om[[taxa_j]]
        y  <- tab[[var_i]]
        ok <- complete.cases(x, y)
        if (sum(ok) < 5) next
        ct <- tryCatch(
          cor.test(x[ok], y[ok], method = method),
          error = function(e) NULL
        )
        if (is.null(ct)) next
        results <- rbind(results, data.frame(
          group    = grp,
          taxa     = taxa_j,
          variable = var_i,
          r        = as.numeric(ct$estimate),
          p        = ct$p.value,
          stringsAsFactors = FALSE
        ))
      }
    }
  }

  if (nrow(results) == 0) stop("No correlation results produced.")

  # ── 8. FDR correction (all tests at once) ─────────────────────────────────
  results$padj <- p.adjust(results$p, method = "BH")
  p_use <- if (isTRUE(padj)) results$padj else results$p
  results$sig <- cut(p_use,
                     breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
                     labels = c("***", "**", "*", ""))
  results$sig[p_use >= pvalue] <- ""

  message(sprintf("[Go_correlation] %d / %d pairs significant (p%s < %.2f).",
                  sum(p_use < pvalue, na.rm = TRUE), nrow(results),
                  if (isTRUE(padj)) "adj" else "", pvalue))

  # ── 9. save table ─────────────────────────────────────────────────────────
  csv_name <- sprintf("%s/correlation.%s.%s.%s%s.csv",
                      out_tab, rank, method,
                      ifelse(is.null(name), "", paste0(name, ".")),
                      format(Sys.Date(), "%y%m%d"))
  write.csv(results, csv_name, row.names = FALSE, quote = FALSE)

  # ── 10. heatmap ───────────────────────────────────────────────────────────
  # order taxa by mean |r| for cleaner visualization
  taxa_order <- results %>%
    dplyr::group_by(taxa) %>%
    dplyr::summarise(mean_r = mean(abs(r), na.rm = TRUE)) %>%
    dplyr::arrange(dplyr::desc(mean_r)) %>%
    dplyr::pull(taxa)
  results$taxa     <- factor(results$taxa,     levels = rev(taxa_order))
  results$variable <- factor(results$variable, levels = colnames(Table))
  results$group    <- factor(results$group,    levels = grp_levels)

  n_facets <- length(grp_levels)
  nc <- if (!is.null(ncols)) ncols else min(n_facets, 4)

  p <- ggplot(results, aes(x = variable, y = taxa, fill = r)) +
    geom_tile(color = "white", linewidth = 0.3) +
    geom_text(aes(label = sig), color = "black", size = 3, vjust = 0.75) +
    scale_fill_gradient2(low  = "#0C5BB0FF",
                         mid  = "white",
                         high = "#EE0011FF",
                         midpoint = 0,
                         limits = c(-1, 1),
                         name = method) +
    labs(x = NULL, y = NULL,
         title = sprintf("Correlation: %s vs external variables", rank),
         caption = sprintf("Filter: top %d taxa, prevalence >= %.0f%% | %s | FDR: %s",
                           top_n, min_prev * 100, method,
                           if (isTRUE(padj)) "BH" else "none")) +
    theme_bw() +
    theme(strip.background  = element_blank(),
          strip.text        = element_text(size = 8, face = "bold"),
          axis.text.x       = element_text(angle = xangle, hjust = 1, vjust = 0.5, size = 7),
          axis.text.y       = element_text(size = 7, face = "italic"),
          plot.title        = element_text(size = 9, face = "bold"),
          plot.caption      = element_text(size = 6, color = "grey50"),
          legend.key.height = unit(0.4, "cm"),
          panel.grid        = element_blank())

  if (n_facets > 1) {
    p <- p + facet_wrap(~ group, ncol = nc, scales = "free_x")
  }

  if (isTRUE(patchwork)) return(invisible(p))
  pdf_name <- sprintf("%s/correlation.%s.%s.%s%s%s.pdf",
                      out_path, rank, method,
                      ifelse(!is.null(group), paste0(group, "."), ""),
                      ifelse(is.null(name), "", paste0(name, ".")),
                      format(Sys.Date(), "%y%m%d"))
  pdf(pdf_name, height = height, width = width)
  print(p)
  dev.off()

  message(sprintf("[Go_correlation] Done. PDF: %s", pdf_name))
  invisible(results)
}
