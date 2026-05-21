#' Mixture Analysis: BKMR + QGCOMP
#'
#' Models how a mixture of variables predicts continuous outcomes using
#' Bayesian Kernel Machine Regression (BKMR) and Quantile G-computation (QGCOMP).
#'
#' mix_vars sources (auto-detected in order):
#'   1. sample_data columns  → used directly, scaled
#'   2. tax_table$ASV_ID     → CLR-transformed from OTU table (16S)
#'   3. taxa_names(psIN)     → CLR-transformed from OTU table (KO, PWY, pathway)
#'
#' Both BKMR and QGCOMP are always run for every outcome.
#'
#' @param psIN Phyloseq object. OTU table provides microbiome/pathway feature
#'   values; sample_data provides outcomes and covariates.
#' @param project Project name string.
#' @param mix_vars Character vector of mixture variable names. Can be any
#'   combination of: ASV_IDs (tax_table$ASV_ID), KO/PWY names (taxa_names),
#'   or sample_data column names (e.g. metals, clinical variables).
#' @param outcomes Character vector of outcome variable names in sample_data.
#'   Required.
#' @param covariates Character vector of covariate names in sample_data, or NULL.
#' @param bkmr_iter Integer. MCMC iterations for BKMR (default 5000).
#' @param seed Integer. Random seed.
#' @param mycols Length-2 color vector: [1] negative, [2] positive.
#' @param panel_width Fixed panel width in inches (y-label area auto-expands;
#'   see GOTOOLS_PLOT_CONVENTIONS.md).
#' @param font Base font size in points.
#' @param name Optional report token from ensure_report_token().
#'
#' @return Invisible named list: results[[outcome]]$bkmr and $qgcomp data frames.
#'
#' @export

Go_mixture <- function(psIN,
                       project,
                       mix_vars,
                       outcomes,
                       covariates  = NULL,
                       bkmr_iter   = 5000,
                       seed        = 42,
                       mycols      = c("#74ADD1", "#E31A1C"),
                       panel_width = 2.5,
                       font        = 10,
                       name        = NULL) {

  # ---- layout helpers (GOTOOLS_PLOT_CONVENTIONS.md) ----
  calc_label_width_in <- function(labels, font_size_pt) {
    labels   <- as.character(labels[!is.na(labels)])
    if (length(labels) == 0) return(1)
    longest  <- labels[which.max(nchar(labels))]
    width_in <- tryCatch(
      grid::convertWidth(
        grid::grobWidth(grid::textGrob(longest, gp = grid::gpar(fontsize = font_size_pt))),
        "in", valueOnly = TRUE),
      error = function(e) NA_real_)
    if (!is.finite(width_in)) width_in <- max(1, nchar(longest) * 0.08)
    width_in
  }

  fix_panel_width_grob <- function(plot_obj, panel_width_in) {
    gt         <- ggplot2::ggplotGrob(plot_obj)
    panel_rows <- gt$layout[gt$layout$name == "panel", , drop = FALSE]
    if (nrow(panel_rows) > 0) {
      panel_cols <- unique(unlist(Map(seq.int, panel_rows$l, panel_rows$r)))
      gt$widths[panel_cols] <- grid::unit(panel_width_in / length(panel_cols), "in")
    }
    gt
  }

  auto_pdf <- function(p, labels, filename, panel_w, font_pt, n_features) {
    plot_height <- max(3.0, min(8.5, 1.2 + 0.15 * n_features))
    grob        <- fix_panel_width_grob(p, panel_w)
    grob_width  <- tryCatch(
      grid::convertWidth(sum(grob$widths), "in", valueOnly = TRUE),
      error = function(e) NA_real_)
    if (!is.finite(grob_width))
      grob_width <- panel_w + max(1.0, min(5.0, calc_label_width_in(labels, font_pt) + 0.3))
    pdf(filename, width = grob_width, height = plot_height)
    grid::grid.draw(grob)
    dev.off()
    message(sprintf("[Go_mixture] %s (%.2f x %.2f in)", basename(filename), grob_width, plot_height))
  }

  # ---- output dirs ----
  dir_root   <- sprintf("%s_%s", project, format(Sys.Date(), "%y%m%d"))
  pdf_dir    <- file.path(dir_root, "pdf")
  tab_dir    <- file.path(dir_root, "table")
  dir.create(pdf_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

  name_token <- if (!is.null(name) && nzchar(name)) paste0(".", name) else ""
  date_str   <- format(Sys.Date(), "%y%m%d")

  # ---- extract metadata and taxonomy ----
  meta <- data.frame(phyloseq::sample_data(psIN), check.names = FALSE)
  tt   <- as.data.frame(phyloseq::tax_table(psIN), stringsAsFactors = FALSE)
  tnames <- phyloseq::taxa_names(psIN)

  missing_out <- setdiff(outcomes, colnames(meta))
  if (length(missing_out) > 0)
    stop("[Go_mixture] outcomes not found in sample_data: ", paste(missing_out, collapse = ", "))

  # ---- auto-detect source for each mix_var ----
  # source: "sampledata" | "asv_id" | "taxa_names"
  detect_source <- function(v) {
    if (v %in% colnames(meta))                                   return("sampledata")
    if ("ASV_ID" %in% colnames(tt) && v %in% tt$ASV_ID)         return("asv_id")
    if (v %in% tnames)                                           return("taxa_names")
    stop("[Go_mixture] '", v, "' not found in sample_data, tax_table$ASV_ID, or taxa_names.")
  }
  sources <- setNames(vapply(mix_vars, detect_source, character(1)), mix_vars)
  message(sprintf("[Go_mixture] mix_vars sources detected:"))
  for (s in unique(sources))
    message(sprintf("  %s: %s", s, paste(names(sources)[sources == s], collapse = ", ")))

  # ---- CLR matrix for OTU-derived variables ----
  otu_vars <- names(sources)[sources %in% c("asv_id", "taxa_names")]
  clr_cols <- NULL   # named character: mix_var -> column name in CLR matrix

  if (length(otu_vars) > 0) {
    ps_clr  <- microbiome::transform(psIN, "clr")
    otu_clr <- as(phyloseq::otu_table(ps_clr), "matrix")
    if (phyloseq::taxa_are_rows(ps_clr)) otu_clr <- t(otu_clr)
    clr_mat <- as.data.frame(otu_clr, check.names = FALSE)

    clr_cols <- character(length(otu_vars))
    names(clr_cols) <- otu_vars
    for (v in otu_vars) {
      if (sources[v] == "asv_id") {
        seq_id <- rownames(tt)[match(v, tt$ASV_ID)]
        if (is.na(seq_id) || !seq_id %in% colnames(clr_mat))
          stop("[Go_mixture] ASV_ID '", v, "' could not be mapped to OTU table.")
        clr_cols[v] <- seq_id
      } else {
        # taxa_names: var IS the column name
        if (!v %in% colnames(clr_mat))
          stop("[Go_mixture] '", v, "' not found in CLR matrix columns.")
        clr_cols[v] <- v
      }
    }
  }

  # ---- feature labels for plots ----
  make_label <- function(v) {
    src <- sources[v]
    if (src == "sampledata") return(v)
    if (src == "asv_id") {
      sp_col <- if ("Species" %in% colnames(tt)) "Species" else NULL
      if (!is.null(sp_col)) {
        seq_id <- rownames(tt)[match(v, tt$ASV_ID)]
        sp <- trimws(as.character(tt[seq_id, sp_col]))
        if (nzchar(sp) && !sp %in% c("NA", "NA NA"))
          return(sprintf("%s (%s)", v, sp))
      }
      return(v)
    }
    # taxa_names: truncate long KO/PWY names
    if (nchar(v) > 45) return(paste0(substr(v, 1, 42), "..."))
    return(v)
  }
  feature_labels <- setNames(vapply(mix_vars, make_label, character(1)), mix_vars)

  # ---- align samples and build Z matrix ----
  all_samps <- rownames(meta)
  if (length(otu_vars) > 0) all_samps <- intersect(all_samps, rownames(clr_mat))

  # build Z: columns in mix_vars order
  Z_list <- lapply(mix_vars, function(v) {
    if (sources[v] == "sampledata") {
      vals <- as.numeric(meta[all_samps, v])
      as.numeric(scale(vals))
    } else {
      as.numeric(clr_mat[all_samps, clr_cols[v]])
    }
  })
  Z <- do.call(cbind, Z_list)
  rownames(Z) <- all_samps
  colnames(Z) <- mix_vars

  meta_sub <- meta[all_samps, , drop = FALSE]

  # covariate matrix
  if (!is.null(covariates)) {
    cov_mat <- model.matrix(
      as.formula(paste("~", paste(covariates, collapse = " + "))),
      data = meta_sub)[, -1, drop = FALSE]
  } else {
    cov_mat <- matrix(0, nrow = nrow(meta_sub), ncol = 1,
                      dimnames = list(rownames(meta_sub), "intercept_only"))
  }

  message(sprintf("[Go_mixture] %d samples | %d mix_vars | %d outcomes",
                  length(all_samps), length(mix_vars), length(outcomes)))

  results <- list()

  for (ny in outcomes) {
    y_ny         <- as.numeric(meta_sub[[ny]])
    complete_idx <- !is.na(y_ny) & apply(!is.na(Z), 1, all)
    message(sprintf("  Outcome [%s]: %d complete cases", ny, sum(complete_idx)))
    if (sum(complete_idx) < 10) {
      message("  -> Skipping: fewer than 10 complete cases.")
      next
    }

    y_c  <- y_ny[complete_idx]
    Z_c  <- Z[complete_idx, , drop = FALSE]
    X_c  <- cov_mat[complete_idx, , drop = FALSE]

    ############################################################
    #=========  BKMR  =========#
    ############################################################
    set.seed(seed)
    fit_b <- tryCatch(
      bkmr::kmbayes(y = y_c, Z = Z_c, X = X_c,
                    iter = bkmr_iter, verbose = FALSE, varsel = TRUE),
      error = function(e) { message("  BKMR error: ", e$message); NULL }
    )

    if (!is.null(fit_b)) {
      pip_df <- data.frame(
        Variable = mix_vars,
        Label    = as.character(feature_labels[mix_vars]),
        Source   = as.character(sources[mix_vars]),
        PIP      = bkmr::ExtractPIPs(fit_b)$PIP,
        row.names = NULL
      )
      results[[ny]]$bkmr <- pip_df

      p_b <- ggplot2::ggplot(pip_df,
               ggplot2::aes(x = PIP, y = reorder(Label, PIP), fill = Source)) +
        ggplot2::geom_col(color = "white") +
        ggplot2::scale_fill_manual(
          values = setNames(
            scales::hue_pal()(length(unique(pip_df$Source))),
            sort(unique(pip_df$Source))
          )
        ) +
        ggplot2::geom_vline(xintercept = 0.5, linetype = "dashed", color = "navy") +
        ggplot2::theme_bw(base_size = font) +
        ggplot2::labs(title    = sprintf("BKMR: Mixture → %s", ny),
                      subtitle = "PIP > 0.5: important contributor",
                      x        = "Posterior Inclusion Probability (PIP)",
                      y        = NULL, fill = "Source")

      auto_pdf(p_b, pip_df$Label,
               file.path(pdf_dir, sprintf("%s.mixture.bkmr_%s%s.%s.pdf",
                                          project, ny, name_token, date_str)),
               panel_width, font, nrow(pip_df))

      write.csv(pip_df,
                file.path(tab_dir, sprintf("%s.mixture.bkmr_%s%s.%s.csv",
                                           project, ny, name_token, date_str)),
                row.names = FALSE, quote = FALSE)
    }

    ############################################################
    #=========  QGCOMP  =========#
    ############################################################
    # quantile-discretise Z for qgcomp (q=4 quartiles)
    Z_qg <- as.data.frame(
      apply(Z_c, 2, function(x) {
        breaks <- unique(quantile(x, probs = seq(0, 1, 0.25), na.rm = TRUE))
        if (length(breaks) < 2) return(rep(0L, length(x)))
        as.integer(cut(x, breaks = breaks, include.lowest = TRUE)) - 1L
      })
    )
    qg_df       <- Z_qg
    qg_df$y_out <- y_c
    if (!is.null(covariates))
      for (cv in covariates) qg_df[[cv]] <- meta_sub[complete_idx, cv]

    cov_part <- if (!is.null(covariates))
      paste0(" + ", paste(covariates, collapse = " + ")) else ""
    safe_vars <- make.names(mix_vars)
    colnames(qg_df)[seq_along(mix_vars)] <- safe_vars
    frm_qg <- as.formula(paste0("y_out ~ ", paste(safe_vars, collapse = " + "), cov_part))

    fit_q <- tryCatch(
      qgcomp::qgcomp.noboot(frm_qg, expnms = safe_vars,
                            data = qg_df, family = gaussian(), q = NULL),
      error = function(e) { message("  QGCOMP error: ", e$message); NULL }
    )

    if (!is.null(fit_q)) {
      build_wt <- function(wts, dir) {
        if (length(wts) == 0) return(data.frame())
        data.frame(
          Variable  = mix_vars[match(names(wts), safe_vars)],
          Label     = as.character(feature_labels[mix_vars[match(names(wts), safe_vars)]]),
          Source    = as.character(sources[mix_vars[match(names(wts), safe_vars)]]),
          Weight    = if (dir == "Positive") wts else -wts,
          Direction = dir,
          row.names = NULL
        )
      }
      wt_df <- rbind(build_wt(fit_q$pos.weights, "Positive"),
                     build_wt(fit_q$neg.weights, "Negative"))

      results[[ny]]$qgcomp <- list(
        psi     = fit_q$psi,
        weights = wt_df
      )

      p_q <- ggplot2::ggplot(wt_df,
               ggplot2::aes(x = Weight, y = reorder(Label, abs(Weight)), fill = Direction)) +
        ggplot2::geom_col(color = "white") +
        ggplot2::scale_fill_manual(
          values = c("Positive" = mycols[2], "Negative" = mycols[1])
        ) +
        ggplot2::theme_bw(base_size = font) +
        ggplot2::labs(title    = sprintf("QGCOMP: Mixture → %s", ny),
                      subtitle = sprintf("Overall ψ = %.3f", fit_q$psi),
                      x        = "Component weight",
                      y        = NULL)

      auto_pdf(p_q, wt_df$Label,
               file.path(pdf_dir, sprintf("%s.mixture.qgcomp_%s%s.%s.pdf",
                                          project, ny, name_token, date_str)),
               panel_width, font, nrow(wt_df))

      write.csv(wt_df,
                file.path(tab_dir, sprintf("%s.mixture.qgcomp_%s%s.%s.csv",
                                           project, ny, name_token, date_str)),
                row.names = FALSE, quote = FALSE)
    }
  }

  invisible(results)
}
