#' Stacked bar chart of microbial composition with phylum-aware colouring
#'
#' Produces stacked bar charts of taxonomic relative or absolute abundance,
#' one bar per sample, coloured by a phylum-aware palette that is built
#' automatically from the phyla present in the input data (up to 10 distinct
#' hues, with subtle within-phylum shading for multiple taxa).
#'
#' Compared with `Go_barchart()`, this function emphasises taxonomic context:
#' taxa from the same phylum receive tightly grouped shades of the same base
#' colour so that phylogenetic groupings are immediately visible. The legend/panel layout
#' system is shared with `Go_barchart()`.
#'
#' @param psIN Phyloseq object.
#' @param cate.vars Character vector of metadata variables used to facet the
#'   plot (one page per variable).
#' @param project Project name used for output directory and file naming.
#' @param taxanames Character vector of taxonomic ranks to visualise
#'   (e.g. `c("Phylum", "Genus")`).  Each rank produces a separate PDF.
#' @param orders Optional character vector specifying factor level order for
#'   group/facet variables and the x-axis label.
#' @param simple Logical.  If `TRUE`, produces unfaceted plots.  Default
#'   `FALSE`.
#' @param relative Logical.  `TRUE` (default) plots relative abundance;
#'   `FALSE` plots raw counts.  The cutoff is always evaluated on relative
#'   abundance regardless of this setting.
#' @param y_axis Optional y-axis label.  When `NULL`, defaults to
#'   `"Relative abundance"` or `"Absolute abundance"`.
#' @param x_label Metadata column to use as the x-axis.  Defaults to
#'   `"SampleIDfactor"` (samples ordered by dominant-taxon abundance).
#' @param show_x_text Logical; if `TRUE`, shows x-axis sample labels. Defaults
#'   to `FALSE` because these plots are typically used for global pattern
#'   inspection rather than per-sample identification.
#' @param facet Character vector of additional faceting variables combined
#'   with each element of `cate.vars`.
#' @param legend Legend position passed to `ggplot2::theme()`.  Default
#'   `"bottom"`.  (Layout is always computed automatically.)
#' @param cutoff Minimum mean relative abundance threshold; taxa below this
#'   are collapsed into `"[1_#Other]"`.  Default `0.005`.
#' @param name Optional string appended to output file names.
#' @param ncol Optional integer fixing the number of facet columns.  When
#'   `NULL` (default) the layout is determined automatically.
#' @param stack_from Character scalar controlling stack direction.
#'   `"bottom"` (default) places the most abundant taxon at the bottom of each
#'   bar so the dominant taxa form a stable baseline across all samples.
#'   `"top"` reverses this so the most abundant taxon is drawn at the top.
#' @param order_other Character scalar controlling whether samples with larger
#'   `"[1_#Other]"` fractions are nudged toward one side within the existing
#'   phylum backbone. `"none"` (default) keeps the original order,
#'   `"right"` moves higher-Other samples to the right, `"left"` moves them to
#'   the left, and `"auto"` uses `A`-style ordering for `Phylum/Class/Order`
#'   and `B`-style ordering for `Family/Genus` and finer ranks.
#' @param mark Optional metadata column name.  Samples where this column
#'   equals `"Yes"` receive a `" *"` suffix on the x-axis label.
#' @param height Nominal chart height in inches (legend height is added
#'   automatically).
#' @param width PDF width in inches.
#'
#' @return Called for its side effects (PDF + CSV files written to disk).
#'   Returns `invisible(NULL)`.
#'
#' @examples
#' \dontrun{
#' Go_colbarchart(
#'   psIN      = ps,
#'   cate.vars = c("TreatmentGroup", "TimePoint"),
#'   project   = "MyProject",
#'   taxanames = c("Phylum", "Genus"),
#'   orders    = c("Control", "Treatment"),
#'   relative  = TRUE,
#'   cutoff    = 0.005,
#'   height    = 6,
#'   width     = 10
#' )
#' }
#'
#' @export
Go_colbarchart <- function(psIN,
                            cate.vars,
                            project,
                            taxanames,
                            orders   = NULL,
                            simple   = FALSE,
                            relative = TRUE,
                            y_axis   = NULL,
                            x_label  = NULL,
                            show_x_text = FALSE,
                            facet    = NULL,
                            legend   = "bottom",
                            cutoff   = 0.005,
                            name     = NULL,
                            ncol             = NULL,
                            stack_from       = c("bottom", "top"),
                            order_other      = c("none", "auto", "right", "left"),
                            mark             = NULL,
                            height,
                            width) {

  if (!is.null(dev.list())) dev.off()

  # ── output directories ───────────────────────────────────────────────────────
  out_base <- sprintf("%s_%s", project, format(Sys.Date(), "%y%m%d"))
  out_pdf  <- file.path(out_base, "pdf")
  out_taxa <- file.path(out_base, "table", "taxa")
  for (d in c(out_base, out_pdf, file.path(out_base, "table"), out_taxa)) {
    if (!dir.exists(d)) dir.create(d, recursive = TRUE)
  }

  # ── defaults ─────────────────────────────────────────────────────────────────
  if (is.null(x_label)) x_label <- "SampleIDfactor"
  if (is.function(name)) name <- NULL
  stack_from <- match.arg(stack_from)
  order_other <- match.arg(order_other)

  tt <- try(orders,  TRUE); if (inherits(tt, "try-error")) orders  <- NULL

  # ── normalise rank names (case-insensitive) ──────────────────────────────────
  actual_ranks <- colnames(tax_table(psIN))
  taxanames <- vapply(taxanames, function(tn) {
    m <- which(tolower(actual_ranks) == tolower(tn))
    if (length(m)) actual_ranks[m[1]] else tn
  }, character(1), USE.NAMES = FALSE)

  mapping_sel <- data.frame(sample_data(psIN))

  # ════════════════════════════════════════════════════════════════════════════
  # Internal helpers
  # ════════════════════════════════════════════════════════════════════════════

  # ── phylum-aware colour palette (replaces Go_color) ─────────────────────────
  # Detects phyla present in the data, assigns up to 10 distinct base hues,
  # and applies only subtle within-phylum shading so related taxa stay grouped.
  get_phylum_hue_map <- function(phyla) {
    base_hues <- c(0.40,  # green
                   0.80,  # purple
                   0.60,  # blue
                   0.20,  # orange
                   0.00,  # red
                   0.10,  # brown
                   0.50,  # cyan
                   0.90,  # pink
                   0.70,  # violet
                   0.15)  # yellow

    phyla <- unique(phyla)
    phyla <- phyla[!is.na(phyla) & nzchar(phyla) & phyla != "[Other]"]
    n_ph  <- length(phyla)
    if (n_ph == 0) return(stats::setNames(numeric(0), character(0)))

    if (n_ph > length(base_hues)) {
      base_hues <- seq(0, 1, length.out = n_ph + 1)[seq_len(n_ph)]
    }
    stats::setNames(base_hues[seq_len(n_ph)], phyla)
  }

  assign_phylum_base_colors <- function(phyla) {
    hue_map <- get_phylum_hue_map(phyla)
    if (length(hue_map) == 0) return(stats::setNames(character(0), character(0)))
    stats::setNames(
      vapply(names(hue_map), function(ph) {
        grDevices::hsv(hue_map[[ph]], 0.95, 1.0)
      }, character(1)),
      names(hue_map)
    )
  }

  assign_phylum_colors <- function(phylum_vec, taxa_names, rank_name = NULL) {
    phyla   <- unique(phylum_vec)
    hue_map <- get_phylum_hue_map(phyla)
    n_ph    <- length(hue_map)

    if (n_ph == 0) {
      cols <- grDevices::colorRampPalette(c("#AAAAAA", "#333333"))(length(taxa_names))
      colors <- stats::setNames(cols, taxa_names)
      colors["[1_#Other]"] <- grDevices::hsv(0, 0, 0.75)
      return(colors)
    }

    colors <- stats::setNames(rep(NA_character_, length(taxa_names)), taxa_names)
    rank_settings <- switch(
      rank_name,
      "Phylum"  = list(s_min = 0.98, v_min = 1.00),
      "Class"   = list(s_min = 0.82, v_min = 0.98),
      "Order"   = list(s_min = 0.66, v_min = 0.96),
      "Family"  = list(s_min = 0.50, v_min = 0.94),
      "Genus"   = list(s_min = 0.34, v_min = 0.92),
      "Species" = list(s_min = 0.26, v_min = 0.90),
      "ASV"     = list(s_min = 0.20, v_min = 0.88),
      list(s_min = 0.60, v_min = 0.95)
    )
    for (ph in names(hue_map)) {
      idx    <- which(phylum_vec == ph)
      n      <- length(idx)
      h_vals <- rep(hue_map[[ph]], n)
      s_vals <- seq(rank_settings$s_min, 1.0, length.out = max(n, 1L))
      v_vals <- seq(rank_settings$v_min, 1.0, length.out = max(n, 1L))
      for (j in seq_along(idx)) {
        colors[taxa_names[idx[j]]] <- grDevices::hsv(h_vals[j], s_vals[j], v_vals[j])
      }
    }
    colors["[1_#Other]"] <- grDevices::hsv(0, 0, 0.75)
    colors[!is.na(colors)]  # drop any NAs (unmapped)
  }

  # ── legend / panel layout estimator (mirrors Go_barchart) ───────────────────
  estimate_layout <- function(df_plot, tax_var, mvar,
                              facet_vars = NULL, x_var = "SampleIDfactor",
                              user_ncol = NULL) {
    if (!is.null(facet_vars) && mvar %in% facet_vars) return(NULL)
    fvars      <- if (is.null(facet_vars)) character(0) else
                  setdiff(as.character(facet_vars), "SampleType")
    panel_vars <- unique(c(fvars, if (isTRUE(simple)) character(0) else mvar))

    panel_ids <- if (length(panel_vars) == 0) {
      factor(rep("all", nrow(df_plot)))
    } else {
      interaction(df_plot[, panel_vars, drop = FALSE],
                  drop = TRUE, lex.order = TRUE)
    }

    facet_ncol <- if (!is.null(user_ncol)) max(1L, as.integer(user_ncol)) else NA_integer_

    labels    <- unique(as.character(df_plot[[tax_var]]))
    labels    <- labels[!is.na(labels) & nzchar(labels)]
    cc        <- length(labels)
    max_nc    <- if (cc > 0) max(nchar(labels)) else 0L
    leg_lim   <- if (tax_var %in% c("Genus", "Species")) 3L else 4L
    leg_ncol  <- min(leg_lim, max(1L, cc))
    leg_rows  <- max(1L, ceiling(cc / leg_ncol))
    leg_key   <- if (max_nc > 40) 0.14 else if (max_nc > 25) 0.18 else 0.20

    list(facet_ncol       = facet_ncol,
         legend_ncol      = leg_ncol,
         legend_text_size = 10,
         legend_key_size  = leg_key,
         legend_height    = 0.18 * leg_rows + 0.35)
  }

  # ── legend extraction & panel assembly (mirrors Go_barchart) ────────────────
  extract_legend_grob <- function(p) cowplot::get_legend(p)

  build_combined <- function(plot_obj, legend_grob, chart_h, legend_h) {
    pm <- plot_obj + ggplot2::theme(legend.position = "none")
    if (is.null(legend_grob)) return(pm)
    cowplot::plot_grid(pm, legend_grob,
                       ncol = 1,
                       rel_heights = c(chart_h, legend_h),
                       align = "v", axis = "lr")
  }

  # ── PDF path builder ─────────────────────────────────────────────────────────
  build_pdf_path <- function(rank_name) {
    rank_tag  <- if (length(taxanames) > 1) paste0(rank_name, ".") else ""
    abund_tag <- if (isTRUE(relative)) "relative" else "absolute"
    sprintf("%s/colbarchart.%s.%s.%s%s(%s).%s%s.%s.pdf",
            out_pdf, abund_tag, project,
            ifelse(is.null(facet), "", paste0(paste(facet, collapse = "."), ".")),
            ifelse(is.null(name),  "", paste0(name, ".")),
            cutoff, rank_tag, format(Sys.Date(), "%y%m%d"), stack_from)
  }

  # ════════════════════════════════════════════════════════════════════════════
  # Main loop — one PDF per taxonomic rank
  # ════════════════════════════════════════════════════════════════════════════
  for (i in seq_along(taxanames)) {
    rank <- taxanames[i]
    order_other_rank <- if (identical(order_other, "auto")) {
      if (rank %in% c("Family", "Genus", "Species", "ASV")) "right" else "none"
    } else {
      order_other
    }

    # ── OTU table (handle taxa-as-rows or taxa-as-cols) ──────────────────────
    otu_raw <- as.data.frame(otu_table(psIN))
    test_tt <- try(
      otu_raw[[rank]] <- getTaxonomy(
        otus = rownames(otu_raw), tax_tab = tax_table(psIN),
        taxRanks = colnames(tax_table(psIN)), level = rank),
      TRUE
    )
    if (inherits(test_tt, "try-error")) {
      otu_raw        <- as.data.frame(t(otu_table(psIN)))
      otu_raw[[rank]] <- getTaxonomy(
        otus = rownames(otu_raw), tax_tab = tax_table(psIN),
        taxRanks = colnames(tax_table(psIN)), level = rank)
    }

    otu_raw$PhylumCol <- getTaxonomy(
      otus = rownames(otu_raw), tax_tab = tax_table(psIN),
      taxRanks = colnames(tax_table(psIN)), level = "Phylum")

    if (ncol(otu_raw) <= 2) next

    # ── aggregate ────────────────────────────────────────────────────────────
    value_cols <- setdiff(colnames(otu_raw), c(rank, "PhylumCol"))
    group_key <- interaction(otu_raw[[rank]], otu_raw$PhylumCol, drop = TRUE, lex.order = TRUE)
    agg_mat <- rowsum(as.matrix(otu_raw[, value_cols, drop = FALSE]), group = group_key, reorder = FALSE)
    group_map <- unique(data.frame(
      group_key = group_key,
      rank = otu_raw[[rank]],
      PhylumCol = otu_raw$PhylumCol,
      stringsAsFactors = FALSE
    ))
    group_map <- group_map[match(rownames(agg_mat), group_map$group_key), ]
    agg <- data.frame(
      rank = group_map$rank,
      PhylumCol = group_map$PhylumCol,
      agg_mat,
      row.names = NULL,
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
    genera    <- agg$rank
    PhylumCol <- agg$PhylumCol
    agg$rank      <- NULL
    agg$PhylumCol <- NULL

    # ── cutoff (always on relative) ──────────────────────────────────────────
    agg_rel              <- normalizeByCols(agg)
    inds_grey            <- which(rowMeans(agg_rel) < cutoff)
    genera[inds_grey]    <- "[1_#Other]"
    PhylumCol[inds_grey] <- "[Other]"

    # ── choose values for plotting ────────────────────────────────────────────
    agg_plot              <- if (isTRUE(relative)) agg_rel else agg
    agg_plot[[rank]]      <- genera
    agg_plot$PhylumCol    <- PhylumCol

    # ── save taxa table (relative, non-Other) ────────────────────────────────
    agg_save        <- agg_rel
    agg_save[[rank]] <- genera
    agg_save         <- subset(agg_save, agg_save[[rank]] != "[1_#Other]")
    write.csv(agg_save, quote = FALSE, row.names = FALSE,
              file = sprintf("%s/%s.taxa_relative_abundance.(%s).%s.%s%s.%s.csv",
                             out_taxa, project, cutoff, rank,
                             ifelse(is.null(name), "", paste0(name, ".")),
                             format(Sys.Date(), "%y%m%d"), stack_from))

    # ── sample ordering (top_taxa) ────────────────────────────────────────────
    # Always derive a reference sample order from Phylum-level composition.
    # Finer ranks (Class/Order/Family/...) reuse that phylum order so their
    # panels keep the same broad structure instead of re-shuffling by sub-rank.
    non_other <- genera != "[1_#Other]"
    agg_no    <- agg_rel[non_other, , drop = FALSE]
    gen_no    <- genera[non_other]
    phy_no    <- PhylumCol[non_other]

    phy_group_key <- interaction(otu_raw$PhylumCol, drop = TRUE, lex.order = TRUE)
    phy_mat_raw <- rowsum(as.matrix(otu_raw[, value_cols, drop = FALSE]), group = phy_group_key, reorder = FALSE)
    phy_rel_ref <- normalizeByCols(phy_mat_raw)
    phylum_totals <- if (nrow(phy_rel_ref) > 0) rowSums(phy_rel_ref) else numeric(0)
    phylum_order  <- names(sort(phylum_totals, decreasing = TRUE))


    taxon_totals <- if (nrow(agg_no) > 0) tapply(rowSums(agg_no), gen_no, sum) else numeric(0)
    taxon_info <- unique(data.frame(
      taxon  = gen_no,
      phylum = phy_no,
      stringsAsFactors = FALSE
    ))
    taxon_info$phylum <- factor(taxon_info$phylum, levels = phylum_order)
    taxon_info$taxon_total <- unname(taxon_totals[taxon_info$taxon])
    taxon_info <- taxon_info[order(taxon_info$phylum, -taxon_info$taxon_total, taxon_info$taxon), ]
    taxon_order <- taxon_info$taxon

    ord_df <- data.frame(
      SampleID = colnames(agg_rel),
      dom_phylum = factor(rep(NA_character_, ncol(agg_rel)), levels = phylum_order),
      dom_phy_abund = rep(0, ncol(agg_rel)),
      stringsAsFactors = FALSE
    )
    sample_order_ref <- if (nrow(phy_rel_ref) == 0) {
      colnames(agg_rel)
    } else {
      phy_mat <- phy_rel_ref[phylum_order[phylum_order %in% rownames(phy_rel_ref)], , drop = FALSE]
      dom_phy_idx   <- apply(phy_mat, 2, which.max)
      dom_phylum    <- rownames(phy_mat)[dom_phy_idx]
      dom_phy_abund <- apply(phy_mat, 2, max)
      ord_df <- data.frame(
        SampleID  = colnames(agg_rel),
        dom_phylum = factor(dom_phylum, levels = phylum_order),
        dom_phy_abund = dom_phy_abund,
        stringsAsFactors = FALSE
      )
      phy_abund_df <- as.data.frame(t(phy_mat))
      colnames(phy_abund_df) <- paste0("phy_", make.names(colnames(phy_abund_df), unique = TRUE))
      ord_df <- cbind(ord_df, phy_abund_df[ord_df$SampleID, , drop = FALSE])
      order_args <- c(
        list(ord_df$dom_phylum, -ord_df$dom_phy_abund),
        lapply(rev(colnames(phy_abund_df)), function(nm) -ord_df[[nm]])
      )
      ord_df <- ord_df[do.call(order, order_args), ]
      ord_df$SampleID
    }

    sample_order <- sample_order_ref
    if (!identical(order_other_rank, "none")) {
      other_prop <- if (any(genera == "[1_#Other]")) {
        colSums(agg_rel[genera == "[1_#Other]", , drop = FALSE])
      } else {
        stats::setNames(rep(0, ncol(agg_rel)), colnames(agg_rel))
      }
      ord_other_df <- data.frame(
        SampleID = sample_order_ref,
        ref_rank = seq_along(sample_order_ref),
        dom_phylum = ord_df$dom_phylum[match(sample_order_ref, ord_df$SampleID)],
        other_prop = unname(other_prop[sample_order_ref]),
        stringsAsFactors = FALSE
      )
      block_size <- 12L
      ord_other_df$within_phylum_rank <- ave(
        ord_other_df$ref_rank,
        ord_other_df$dom_phylum,
        FUN = seq_along
      )
      ord_other_df$ref_block <- ave(
        ord_other_df$within_phylum_rank,
        ord_other_df$dom_phylum,
        FUN = function(x) (seq_along(x) - 1L) %/% block_size
      )
      ord_other_df <- ord_other_df[
        order(
          ord_other_df$dom_phylum,
          ord_other_df$ref_block,
          if (identical(order_other_rank, "right")) ord_other_df$other_prop else -ord_other_df$other_prop,
          ord_other_df$ref_rank
        ),
      ]
      sample_order <- ord_other_df$SampleID
    }

    # ── melt to long format ───────────────────────────────────────────────────
    df <- reshape2::melt(agg_plot, variable.name = "SampleID")
    long_key <- interaction(df[[rank]], df$PhylumCol, df$SampleID, drop = TRUE, lex.order = TRUE)
    df2_mat <- rowsum(matrix(df$value, ncol = 1), group = long_key, reorder = FALSE)
    long_map <- unique(data.frame(
      long_key = long_key,
      rank = df[[rank]],
      PhylumCol = df$PhylumCol,
      SampleID = df$SampleID,
      stringsAsFactors = FALSE
    ))
    long_map <- long_map[match(rownames(df2_mat), long_map$long_key), ]
    df2 <- data.frame(
      rank = long_map$rank,
      PhylumCol = long_map$PhylumCol,
      SampleID = long_map$SampleID,
      value = as.numeric(df2_mat[, 1]),
      row.names = NULL,
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
    colnames(df2)[colnames(df2) == "rank"] <- rank
    df2$SampleID       <- as.character(df2$SampleID)
    df2$SampleIDfactor <- factor(df2$SampleID, levels = sample_order)

    # ── add metadata columns ─────────────────────────────────────────────────
    for (mvar in cate.vars) {
      df2[[mvar]] <- mapping_sel[df2$SampleID, mvar]
      df2[[mvar]] <- if (!is.null(orders)) factor(df2[[mvar]], levels = orders) else
                     factor(df2[[mvar]])
    }
    if (!is.null(facet)) {
      for (fa in facet) {
        df2[[fa]] <- mapping_sel[df2$SampleID, fa]
        df2[[fa]] <- if (!is.null(orders)) factor(df2[[fa]], levels = orders) else
                     factor(df2[[fa]])
      }
    }

    # ── x-axis label + optional mark ─────────────────────────────────────────
    if (!x_label %in% c("SampleID", "SampleIDfactor")) {
      df2[[x_label]] <- mapping_sel[df2$SampleID, x_label]
      if (!is.null(mark)) {
        df2[[mark]] <- mapping_sel[df2$SampleID, mark]
        df2[[x_label]] <- ifelse(
          df2[[mark]] == "Yes",
          paste0(df2[[x_label]], " *"),
          as.character(df2[[x_label]])
        )
        new_orders <- unlist(lapply(orders, function(o) c(o, paste0(o, " *"))))
        df2[[x_label]] <- factor(df2[[x_label]], levels = new_orders)
      } else {
        df2[[x_label]] <- factor(df2[[x_label]], levels = orders)
      }
    }

    # ── clean taxa name strings ───────────────────────────────────────────────
    df2[[rank]] <- gsub("\\s*\\(ex [^)]+\\)", "", as.character(df2[[rank]]))

    # ── colour palette ────────────────────────────────────────────────────────
    cdf_sel <- unique(df2[df2[[rank]] != "[1_#Other]", c("PhylumCol", rank)])
    cdf_sel <- cdf_sel[order(cdf_sel$PhylumCol), ]

    coloring <- assign_phylum_colors(
      phylum_vec = cdf_sel$PhylumCol,
      taxa_names = cdf_sel[[rank]],
      rank_name = rank
    )

    # Keep stack order grouped by Phylum so same-colour families/classes stay
    # contiguous inside each bar. position_stack(reverse=) still controls
    # whether the dominant blocks sit near the top or the bottom.
    non_other_mask  <- df2[[rank]] != "[1_#Other]"
    taxon_abund_sum <- tapply(df2$value[non_other_mask],
                              df2[[rank]][non_other_mask], sum)
    taxon_stack_map <- unique(df2[non_other_mask, c(rank, "PhylumCol")])
    taxon_stack_map$PhylumCol <- factor(taxon_stack_map$PhylumCol, levels = phylum_order)
    taxon_stack_map$taxon_abund <- unname(taxon_abund_sum[as.character(taxon_stack_map[[rank]])])
    taxon_stack_map <- taxon_stack_map[
      order(taxon_stack_map$PhylumCol, -taxon_stack_map$taxon_abund, taxon_stack_map[[rank]]),
    ]
    level_order <- c(as.character(taxon_stack_map[[rank]]), "[1_#Other]")
    df2$.fill_var <- factor(df2[[rank]], levels = level_order)

    # ── layout plan ──────────────────────────────────────────────────────────
    layout_plan <- list()
    for (mvar in cate.vars) {
      li <- estimate_layout(df2, rank, mvar, facet, x_label, ncol)
      if (!is.null(li)) layout_plan[[mvar]] <- li
    }
    if (length(layout_plan) == 0) next

    pdf_h <- max(vapply(layout_plan, function(x) height + x$legend_height, numeric(1))) +
             if (!identical(rank, "Phylum")) 0.25 else 0
    grDevices::pdf(build_pdf_path(rank), height = pdf_h, width = width)

    # ── plot per cate.var ─────────────────────────────────────────────────────
    for (mvar in cate.vars) {
      cl <- layout_plan[[mvar]]
      if (is.null(cl)) next
      rank_legend_ncol <- if (!identical(rank, "Phylum")) 3L else cl$legend_ncol

      p <- ggplot2::ggplot(
        df2,
        ggplot2::aes(x    = .data[[x_label]],
                     y    = value,
                     fill = .fill_var)
      ) +
        ggplot2::geom_bar(
          stat = "identity",
          position = ggplot2::position_stack(
            reverse = identical(stack_from, "top")
          )
        ) +
        ggplot2::theme_classic() +
        ggplot2::scale_fill_manual(values = coloring, name = rank, drop = FALSE) +
        ggplot2::guides(fill = ggplot2::guide_legend(
          ncol = rank_legend_ncol, byrow = TRUE,
          title.position = "top", title.hjust = 0)) +
        ggplot2::theme(
          legend.position   = "bottom",
          legend.text       = ggplot2::element_text(face   = "italic",
                                                    size   = cl$legend_text_size),
          legend.title      = ggplot2::element_text(
                                size = cl$legend_text_size,
                                margin = ggplot2::margin(0, 0, 0.08, 0, "cm")
                              ),
          legend.key.height = grid::unit(cl$legend_key_size, "cm"),
          legend.key.width  = grid::unit(cl$legend_key_size, "cm"),
          legend.spacing.y  = grid::unit(0.02, "cm"),
          legend.margin     = ggplot2::margin(0, 0, 0, 0, "cm"),
          axis.title.x      = ggplot2::element_blank(),
          axis.text.x       = if (isTRUE(show_x_text)) {
                                ggplot2::element_text(angle = 90, vjust = 0.5,
                                                      hjust = 1, size = 8)
                              } else {
                                ggplot2::element_blank()
                              },
          axis.ticks.x      = if (isTRUE(show_x_text)) {
                                ggplot2::element_line()
                              } else {
                                ggplot2::element_blank()
                              }
        ) +
        ggplot2::ggtitle(sprintf(
          "%s barplots of %s%s (cutoff < %s)%s",
          rank, mvar,
          ifelse(is.null(name), "", paste0("-", name)),
          cutoff,
          ifelse(is.null(mark), "", paste0("  mark: ", mark))
        ))
      if (!identical(rank, "Phylum")) {
        p <- p + ggplot2::theme(
          legend.justification = "left",
          legend.box.just = "left",
          legend.box.spacing = grid::unit(0, "cm"),
          legend.box.margin = ggplot2::margin(0, 0, 0, -0.45, "cm"),
          legend.spacing.x = grid::unit(0, "cm")
        )
      }

      # y-axis label
      if (!is.null(y_axis)) {
        p <- p + ggplot2::labs(y = y_axis)
      } else if (isTRUE(relative)) {
        p <- p + ggplot2::labs(y = "Relative abundance") +
                 ggplot2::ylim(c(-0.1, 1.01))
      } else {
        p <- p + ggplot2::labs(y = "Absolute abundance")
      }

      # faceting
      if (!isTRUE(simple)) {
        if (!is.null(facet)) {
          fvars     <- setdiff(as.character(facet), "SampleType")
          f_formula <- as.formula(
            sprintf("~ %s", paste(c(fvars, mvar), collapse = " + ")))
          if (is.na(cl$facet_ncol)) {
            p <- p + ggplot2::facet_grid(f_formula, scales = "free_x",
                                         space = "free_x")
          } else {
            p <- p + ggplot2::facet_wrap(f_formula, scales = "free_x",
                                         ncol = cl$facet_ncol)
          }
        } else {
          f_formula <- as.formula(sprintf("~ %s", mvar))
          if (is.na(cl$facet_ncol)) {
            p <- p + ggplot2::facet_grid(f_formula, scales = "free_x",
                                         space = "free_x")
          } else {
            p <- p + ggplot2::facet_wrap(f_formula, scales = "free_x",
                                         ncol = cl$facet_ncol)
          }
        }
      }

      legend_grob <- extract_legend_grob(p)
      if (!identical(rank, "Phylum")) {
        phylum_levels <- phylum_order[phylum_order %in% unique(as.character(df2$PhylumCol))]
        phylum_colors <- assign_phylum_base_colors(phylum_levels)
        if (length(phylum_colors) > 0) {
          phylum_legend_df <- data.frame(
            PhylumCol = factor(names(phylum_colors), levels = names(phylum_colors)),
            x = 1,
            y = 1,
            stringsAsFactors = FALSE
          )
          p_phy_legend <- ggplot2::ggplot(
            phylum_legend_df,
            ggplot2::aes(x = x, y = y, fill = PhylumCol)
          ) +
            ggplot2::geom_col() +
            ggplot2::scale_fill_manual(values = phylum_colors, name = "Phylum", drop = FALSE) +
            ggplot2::guides(fill = ggplot2::guide_legend(
              ncol = 2L, byrow = TRUE,
              title.position = "top", title.hjust = 0)) +
            ggplot2::theme_void() +
            ggplot2::theme(
              legend.position   = "bottom",
              legend.justification = "left",
              legend.box.just = "left",
              legend.box.margin = ggplot2::margin(0, 0, 0, 0, "cm"),
              legend.title      = ggplot2::element_text(
                                    size = cl$legend_text_size,
                                    margin = ggplot2::margin(0, 0, 0.08, 0, "cm")
                                  ),
              legend.text       = ggplot2::element_text(face = "italic",
                                                        size = cl$legend_text_size),
              legend.key.height = grid::unit(cl$legend_key_size, "cm"),
              legend.key.width  = grid::unit(cl$legend_key_size, "cm"),
              legend.spacing.x  = grid::unit(0.02, "cm"),
              legend.spacing.y  = grid::unit(0.02, "cm"),
              legend.margin     = ggplot2::margin(0, 0, 0, 0, "cm")
            )
          phylum_legend_grob <- extract_legend_grob(p_phy_legend)
          phy_legend_rows <- max(1L, ceiling(length(phylum_colors) / 2L))
          phy_raise <- max(0.16, 0.52 - 0.08 * phy_legend_rows)
          phylum_legend_panel <- cowplot::ggdraw() +
            cowplot::draw_grob(
              phylum_legend_grob,
              x = 0, y = phy_raise, width = 1, height = 1 - phy_raise,
              hjust = 0, vjust = 0
            )
          legend_grob <- cowplot::plot_grid(
            phylum_legend_panel, legend_grob,
            nrow = 1,
            rel_widths = c(0.66, 1.34),
            align = "v",
            axis = "tb"
          )
        }
      }
      combined    <- build_combined(p, legend_grob, height,
                                    cl$legend_height + if (!identical(rank, "Phylum")) 0.10 else 0)
      print(combined)
    }
    grDevices::dev.off()
  }
  invisible(NULL)
}
