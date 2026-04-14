#' Generate NMDS Biplot for Microbiome Data
#'
#' @param psIN Phyloseq object containing microbiome data.
#' @param project Name of the project for labeling.
#' @param cate.vars Optional categorical sample-data variables to plot. If
#'   \code{NULL}, the function auto-detects categorical columns from
#'   \code{sample_data(psIN)}.
#' @param orders Optional ordering of factor levels applied to \code{cate.vars}
#'   and \code{facet} via \code{intersect()}.
#' @param distance_metrics Vector of distance metrics to be used for NMDS.
#' @param data_type Deprecated compatibility argument. Ignored.
#' @param biplot Optional taxa selector for arrow variables. Accepts numeric
#'   indices or taxa names at the aggregated \code{TaxaLevel}. If \code{NULL},
#'   all aggregated taxa are eligible.
#' @param shapes Point shape value. Scalar shapes are drawn directly.
#' @param TaxaLevel Taxonomic rank used for abundance aggregation and taxa arrows.
#' @param ID Optional sample ID column for point labels.
#' @param facet Optional faceting variable from \code{sample_data(psIN)}.
#' @param name Optional name for the output plot.
#' @param height Height of the plot.
#' @param width Width of the plot.
#' @param con.vars Optional numeric sample-data columns to include as
#'   supplementary envfit arrows.
#' @param taxa_allow Logical. If \code{TRUE}, include taxa arrows from the
#'   aggregated \code{TaxaLevel}; if \code{FALSE}, hide taxa arrows.
#'
#' @details
#' This function creates NMDS biplots from a \code{phyloseq} object using
#' \code{sample_data(psIN)} directly. Sample distances are determined only by
#' the aggregated microbiome abundance matrix at \code{TaxaLevel}. Optional
#' continuous variables in \code{con.vars} are added as supplementary envfit
#' arrows and do not affect the NMDS coordinates themselves.
#'
#' @return A PDF file containing the NMDS biplot.
#'
#' @examples
#' \dontrun{
#' Go_biplot(
#'   psIN = ps_object,
#'   project = "MyMicrobiomeStudy",
#'   cate.vars = c("Treatment"),
#'   orders = c("Control", "Case"),
#'   distance_metrics = c("bray"),
#'   TaxaLevel = "Genus",
#'   facet = "Site",
#'   con.vars = c("Acetate", "Propionate", "Butyrate"),
#'   taxa_allow = TRUE,
#'   height = 6,
#'   width = 8
#' )
#' }
#'
#' @param patchwork Logical. If \code{TRUE}, skip saving and return the plot object(s) for use with \code{Gg_patchwork()} or the \pkg{patchwork} package. Default \code{FALSE}.
#' @export

Go_biplot <- function(psIN,
                      project,
                      cate.vars = NULL,
                      orders = NULL,
                      distance_metrics = c("bray"),
                      data_type = NULL,
                      biplot = NULL,
                      shapes = 16,
                      TaxaLevel = "Genus",
                      ID = NULL,
                      facet = NULL,
                      name = NULL,
                      height = 6,
                      width = 8,
                      con.vars = NULL,
                      taxa_allow = TRUE,
                      patchwork = FALSE) {

  if (!is.null(dev.list())) dev.off()

  palette_base <- c(
    "#1170aa", "#fc7d0b", "#76B7B2", "#E15759", "#59A14F",
    "#EDC948", "#B07AA1", "#FF9DA7", "#9C755F", "#BAB0AC"
  )

  out <- sprintf("%s_%s", project, format(Sys.Date(), "%y%m%d"))
  out_path <- file.path(out, "pdf")
  if (!dir.exists(out)) dir.create(out, recursive = TRUE)
  if (!dir.exists(out_path)) dir.create(out_path, recursive = TRUE)

  if (!inherits(psIN, "phyloseq")) {
    stop("`psIN` must be a phyloseq object.")
  }

  mapping.sel <- as.data.frame(sample_data(psIN), stringsAsFactors = FALSE)
  if (nrow(mapping.sel) == 0) {
    stop("`sample_data(psIN)` is empty.")
  }

  if (!is.null(con.vars) && length(con.vars) > 0) {
    missing_con <- setdiff(con.vars, names(mapping.sel))
    if (length(missing_con) > 0) {
      stop("`con.vars` not found in sample_data(psIN): ", paste(missing_con, collapse = ", "))
    }
    non_numeric <- con.vars[!vapply(mapping.sel[, con.vars, drop = FALSE], is.numeric, logical(1))]
    if (length(non_numeric) > 0) {
      stop("`con.vars` must be numeric columns in sample_data(psIN): ", paste(non_numeric, collapse = ", "))
    }
  }

  detect_categorical_vars <- function(df, exclude = character(0)) {
    keep <- setdiff(names(df), exclude)
    keep[vapply(df[, keep, drop = FALSE], function(x) {
      if (is.factor(x) || is.character(x) || is.logical(x)) {
        return(TRUE)
      }
      if (is.integer(x)) {
        return(length(unique(stats::na.omit(x))) <= 10)
      }
      FALSE
    }, logical(1))]
  }

  normalize_levels <- function(x, orders = NULL) {
    x_chr <- as.character(x)
    x_chr[x_chr == ""] <- NA_character_
    seen <- unique(stats::na.omit(x_chr))
    if (length(seen) == 0) {
      return(factor(x_chr))
    }
    if (is.null(orders) || length(orders) == 0) {
      return(factor(x_chr, levels = seen))
    }
    lvls <- c(intersect(as.character(orders), seen), setdiff(seen, as.character(orders)))
    factor(x_chr, levels = lvls)
  }

  ranks <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  otu_mat <- as(otu_table(psIN), "matrix")
  if (!phyloseq::taxa_are_rows(psIN)) {
    otu_mat <- t(otu_mat)
  }
  otu_df <- as.data.frame(otu_mat, stringsAsFactors = FALSE)
  otu_df$TaxonLabel <- getTaxonomy(
    otus = rownames(otu_df),
    tax_tab = tax_table(psIN),
    taxRanks = ranks,
    level = TaxaLevel
  )
  otu_df$TaxonLabel[is.na(otu_df$TaxonLabel) | otu_df$TaxonLabel == ""] <- rownames(otu_df)[is.na(otu_df$TaxonLabel) | otu_df$TaxonLabel == ""]
  agg <- stats::aggregate(. ~ TaxonLabel, data = otu_df, FUN = sum, na.action = na.pass)
  rownames(agg) <- agg$TaxonLabel
  agg <- agg[, setdiff(names(agg), "TaxonLabel"), drop = FALSE]
  agg_t <- as.data.frame(t(as.matrix(agg)), stringsAsFactors = FALSE)

  resolve_taxa_columns <- function(agg_t, selector) {
    if (is.null(selector) || identical(selector, TRUE)) {
      return(colnames(agg_t))
    }
    if (identical(selector, FALSE)) {
      return(character(0))
    }
    if (is.numeric(selector)) {
      idx <- unique(as.integer(selector))
      idx <- idx[idx >= 1 & idx <= ncol(agg_t)]
      return(colnames(agg_t)[idx])
    }
    selector <- as.character(selector)
    selector[selector %in% colnames(agg_t)]
  }

  build_envfit_table <- function(sample_ids, agg_t, taxa_cols, con_df = NULL) {
    env_parts <- list()
    if (length(taxa_cols) > 0) {
      taxa_df <- agg_t[sample_ids, taxa_cols, drop = FALSE]
      env_parts$taxa <- taxa_df
    }
    if (!is.null(con_df) && ncol(con_df) > 0) {
      env_parts$con <- con_df[sample_ids, , drop = FALSE]
    }
    if (length(env_parts) == 0) {
      return(NULL)
    }
    out <- do.call(cbind, env_parts)
    out <- out[, vapply(out, function(x) stats::var(as.numeric(x), na.rm = TRUE) > 0, logical(1)), drop = FALSE]
    if (ncol(out) == 0) {
      return(NULL)
    }
    out
  }

  build_arrow_df <- function(ordi, envfit_table, taxa_cols, con.vars) {
    if (is.null(envfit_table) || ncol(envfit_table) == 0) {
      return(NULL)
    }
    bio.fit <- vegan::envfit(ordi, envfit_table, perm = 999)
    vec <- vegan::scores(bio.fit, display = "vectors")
    if (is.null(vec) || nrow(vec) == 0) {
      return(NULL)
    }
    vec <- as.data.frame(vec * vegan:::ordiArrowMul(vec))
    vec$label <- rownames(vec)
    vec$type <- ifelse(vec$label %in% con.vars, "con", "taxa")
    rownames(vec) <- NULL
    vec
  }

  exclude_auto <- unique(c(con.vars, ID, facet))
  if (is.null(cate.vars) || length(cate.vars) == 0) {
    cate.vars <- detect_categorical_vars(mapping.sel, exclude = exclude_auto)
  }
  if (length(cate.vars) == 0) {
    stop("No categorical variables found for `cate.vars`. Provide `cate.vars` explicitly.")
  }

  missing_cate <- setdiff(cate.vars, names(mapping.sel))
  if (length(missing_cate) > 0) {
    stop("`cate.vars` not found in sample_data(psIN): ", paste(missing_cate, collapse = ", "))
  }
  if (!is.null(facet) && length(facet) == 1 && !facet %in% names(mapping.sel)) {
    stop("`facet` not found in sample_data(psIN): ", facet)
  }
  if (!is.null(ID) && length(ID) == 1 && !ID %in% names(mapping.sel)) {
    stop("`ID` not found in sample_data(psIN): ", ID)
  }

  pdf_file <- if (!is.null(name) && nzchar(name)) {
    sprintf("12_biplot.%s.%s.%s.pdf", project, name, format(Sys.Date(), "%y%m%d"))
  } else {
    sprintf("12_biplot.%s.%s.pdf", project, format(Sys.Date(), "%y%m%d"))
  }
  plotlist_pw <- list()
  if (!isTRUE(patchwork)) grDevices::pdf(file.path(out_path, pdf_file), height = height, width = width)

  taxa_cols <- if (isTRUE(taxa_allow)) resolve_taxa_columns(agg_t, biplot) else character(0)

  for (mvar in cate.vars) {
    if (!is.null(facet) && length(facet) == 1 && identical(mvar, facet)) {
      next
    }
    for (distance_metric in distance_metrics) {
      mapping_use <- mapping.sel
      mapping_use[[mvar]] <- normalize_levels(mapping_use[[mvar]], orders = orders)
      if (!is.null(facet) && length(facet) == 1) {
        mapping_use[[facet]] <- normalize_levels(mapping_use[[facet]], orders = orders)
      }

      sample_keep <- rownames(mapping_use)[!is.na(mapping_use[[mvar]])]
      if (length(sample_keep) < 3) {
        next
      }

      psIN.na <- phyloseq::prune_samples(sample_keep, psIN)
      mapping_na <- mapping_use[sample_keep, , drop = FALSE]

      message(sprintf("##-- %s / %s (without NA: %d/%d) --##",
                      mvar, distance_metric, nrow(mapping_na), nrow(mapping.sel)))

      ordi <- phyloseq::ordinate(psIN.na, method = "NMDS", distance = distance_metric)
      df_sites <- as.data.frame(vegan::scores(ordi, display = "sites"))
      df_sites$SampleID <- rownames(df_sites)
      df_sites <- merge(df_sites, cbind(SampleID = rownames(mapping_na), mapping_na),
                        by = "SampleID", sort = FALSE)

      env_con <- if (!is.null(con.vars) && length(con.vars) > 0) {
        mapping_na[, con.vars, drop = FALSE]
      } else {
        NULL
      }
      envfit_table <- build_envfit_table(sample_ids = rownames(mapping_na), agg_t = agg_t,
                                         taxa_cols = taxa_cols, con_df = env_con)
      df_biofit <- build_arrow_df(ordi = ordi, envfit_table = envfit_table,
                                  taxa_cols = taxa_cols,
                                  con.vars = if (is.null(con.vars)) character(0) else con.vars)

      mvar_levels <- levels(mapping_use[[mvar]])
      mvar_levels <- mvar_levels[mvar_levels %in% unique(as.character(stats::na.omit(mapping_na[[mvar]])))]
      plot_cols <- rep(palette_base, length.out = max(1, length(mvar_levels)))
      names(plot_cols) <- mvar_levels

      p <- ggplot2::ggplot(df_sites, ggplot2::aes(x = NMDS1, y = NMDS2, colour = .data[[mvar]])) +
        ggplot2::geom_point(shape = shapes, size = 2) +
        ggplot2::stat_ellipse(level = 0.95, type = "t", linetype = 3, linewidth = 0.4, show.legend = FALSE) +
        ggplot2::theme_bw() +
        ggplot2::ggtitle(sprintf("%s - %s - NMDS - %s", project, mvar, distance_metric)) +
        ggplot2::scale_color_manual(values = plot_cols, drop = FALSE)

      if (!is.null(df_biofit) && nrow(df_biofit) > 0) {
        p <- p +
          ggplot2::geom_segment(
            data = df_biofit,
            ggplot2::aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2, linetype = type),
            inherit.aes = FALSE,
            arrow = grid::arrow(length = grid::unit(0.16, "cm")),
            color = "#808080",
            alpha = 0.75
          ) +
          ggplot2::geom_text(
            data = transform(df_biofit, NMDS1 = NMDS1 * 1.08, NMDS2 = NMDS2 * 1.08),
            ggplot2::aes(x = NMDS1, y = NMDS2, label = label),
            inherit.aes = FALSE,
            color = "#606060",
            size = 3
          ) +
          ggplot2::scale_linetype_manual(values = c(taxa = "solid", con = "dashed"), guide = "none")
      }

      if (!is.null(ID) && length(ID) == 1) {
        p <- p + ggrepel::geom_text_repel(ggplot2::aes(label = .data[[ID]]), size = 2, show.legend = FALSE)
      }

      if (!is.null(facet) && length(facet) == 1) {
        ncol_facet <- length(unique(stats::na.omit(mapping_na[[facet]])))
        p <- p + ggplot2::facet_wrap(stats::as.formula(sprintf("~ %s", facet)),
                                     scales = "free_x", ncol = max(1, ncol_facet))
      }

      plotlist_pw[[length(plotlist_pw) + 1]] <- p
      if (!isTRUE(patchwork)) print(p)
    }
  }

  if (isTRUE(patchwork)) return(invisible(plotlist_pw))
  grDevices::dev.off()
  invisible(file.path(out_path, pdf_file))
}
