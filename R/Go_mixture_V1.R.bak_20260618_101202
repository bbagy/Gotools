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
#' @param mycols Length-2 color vector: `[1]` negative, `[2]` positive.
#' @param pip_threshold PIP cutoff for univariate response curves (default 0.5).
#' @param panel_width Fixed panel width in inches (y-label area auto-expands;
#'   see GOTOOLS_PLOT_CONVENTIONS.md).
#' @param font Base font size in points.
#' @param name Optional report token from ensure_report_token().
#'
#' @return Invisible named list: `results[[outcome]]$bkmr`, `$bkmr_response`,
#'   and `$qgcomp` data frames.
#'
#' @export

Go_mixture <- function(psIN,
                       project,
                       mix_vars,
                       outcomes,
                       covariates    = NULL,
                       bkmr_iter     = 5000,
                       seed          = 42,
                       mycols        = c("#74ADD1", "#E31A1C"),
                       pip_threshold = 0.5,
                       panel_width   = 2.5,
                       font          = 10,
                       name          = NULL) {

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
          return(sprintf("<i>%s</i> (%s)", sp, v))
      }
      return(v)
    }
    # taxa_names: truncate long KO/PWY names
    if (nchar(v) > 45) return(paste0(substr(v, 1, 42), "..."))
    return(v)
  }
  feature_labels <- setNames(vapply(mix_vars, make_label, character(1)), mix_vars)

  # plain label (HTML tags stripped) for text repel
  plain_labels <- setNames(gsub("<[^>]+>", "", feature_labels), mix_vars)

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
    X_raw <- cov_mat[complete_idx, , drop = FALSE]
    # drop zero-variance covariates (causes Cholesky failure in BKMR)
    x_sd  <- apply(X_raw, 2, sd, na.rm = TRUE)
    x_ok  <- x_sd > 1e-8
    if (any(!x_ok))
      message(sprintf("  Dropping %d zero-variance covariate(s): %s",
                      sum(!x_ok), paste(colnames(X_raw)[!x_ok], collapse = ", ")))
    X_c <- if (sum(x_ok) > 0) X_raw[, x_ok, drop = FALSE] else NULL

    ############################################################
    #=========  BKMR  =========#
    ############################################################
    # drop near-zero-variance columns from Z before BKMR
    z_sd      <- apply(Z_c, 2, sd, na.rm = TRUE)
    z_ok      <- z_sd > 1e-8
    if (any(!z_ok))
      message(sprintf("  BKMR: dropping %d zero-variance var(s): %s",
                      sum(!z_ok), paste(colnames(Z_c)[!z_ok], collapse = ", ")))
    Z_bkmr <- Z_c[, z_ok, drop = FALSE]
    if (ncol(Z_bkmr) < 2) {
      message("  BKMR: fewer than 2 valid vars after variance filter — skipping.")
      fit_b <- NULL
    } else {
      set.seed(seed)
      fit_b <- tryCatch(
        bkmr::kmbayes(y = y_c, Z = Z_bkmr, X = X_c,
                      iter = bkmr_iter, verbose = FALSE, varsel = TRUE),
        error = function(e) { message("  BKMR error: ", e$message); NULL }
      )
    }

    if (!is.null(fit_b)) {
      bkmr_vars <- colnames(Z_bkmr)
      pip_df <- data.frame(
        Variable = bkmr_vars,
        Label    = as.character(feature_labels[bkmr_vars]),
        Source   = as.character(sources[bkmr_vars]),
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
        ggplot2::theme(axis.text.y = ggtext::element_markdown()) +
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

      # ---- BKMR univariate response curves for high-PIP taxa ----
      high_pip_vars <- pip_df$Variable[pip_df$PIP >= pip_threshold]
      if (length(high_pip_vars) >= 1) {
        message(sprintf("  BKMR response curves: %d taxa (PIP >= %.2f)",
                        length(high_pip_vars), pip_threshold))
        which_z  <- which(colnames(Z_bkmr) %in% high_pip_vars)
        burnin   <- round(bkmr_iter * 0.5)
        resp <- tryCatch(
          bkmr::PredictorResponseUnivar(
            fit     = fit_b,
            which.z = which_z,
            sel     = seq(burnin + 1L, bkmr_iter)
          ),
          error = function(e) { message("  Response curve error: ", e$message); NULL }
        )
        if (!is.null(resp)) {
          resp$label      <- feature_labels[resp$variable]
          resp$plain_label <- plain_labels[resp$variable]
          results[[ny]]$bkmr_response <- resp

          n_p     <- length(unique(resp$variable))
          ncols_p <- min(3L, n_p)
          nrows_p <- ceiling(n_p / ncols_p)

          p_resp <- ggplot2::ggplot(resp,
                     ggplot2::aes(x = z, y = est,
                                  ymin = est - 1.96 * se,
                                  ymax = est + 1.96 * se)) +
            ggplot2::geom_ribbon(fill = "steelblue", alpha = 0.2) +
            ggplot2::geom_line(color = "steelblue", linewidth = 0.8) +
            ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
            ggplot2::facet_wrap(~ label, scales = "free_x", ncol = ncols_p,
                                labeller = ggplot2::label_value) +
            ggplot2::theme_bw(base_size = font) +
            ggplot2::theme(strip.text = ggtext::element_markdown()) +
            ggplot2::labs(
              title    = sprintf("BKMR response curves → %s", ny),
              subtitle = sprintf("PIP ≥ %.2f  |  other components fixed at median", pip_threshold),
              x        = "Mixture component (CLR, scaled)",
              y        = "h(z)"
            )

          pdf_h <- max(3.5, nrows_p * 2.8 + 1.2)
          pdf_w <- max(4.0, ncols_p * (panel_width + 0.8) + 1.0)
          fname_r <- file.path(pdf_dir,
                               sprintf("%s.mixture.bkmr_response_%s%s.%s.pdf",
                                       project, ny, name_token, date_str))
          pdf(fname_r, width = pdf_w, height = pdf_h)
          print(p_resp)
          dev.off()
          message(sprintf("[Go_mixture] %s (%.2f x %.2f in)",
                          basename(fname_r), pdf_w, pdf_h))

          write.csv(resp[, c("variable","plain_label","z","est","se")],
                    file.path(tab_dir,
                              sprintf("%s.mixture.bkmr_response_%s%s.%s.csv",
                                      project, ny, name_token, date_str)),
                    row.names = FALSE, quote = FALSE)
        }
      }
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
        ggplot2::theme(axis.text.y = ggtext::element_markdown()) +
        ggplot2::labs(title    = sprintf("QGCOMP: Mixture → %s", ny),
                      subtitle = sprintf("Overall psi = %.3f", fit_q$psi),
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


#' Compare BKMR PIP across stratified groups
#'
#' Takes named list of Go_mixture() results (one per group) and produces
#' a dot plot with Grp1 PIP on x-axis and Grp2 PIP on y-axis per outcome.
#' Points above the diagonal are more important in Grp2; below in Grp1.
#'
#' @param results_list Named list of Go_mixture() return values, e.g.
#'   `list(Grp1 = res1, Grp2 = res2)`.
#' @param project Project name string (for output filenames).
#' @param outcomes Character vector of outcome names to plot. NULL = all.
#' @param pip_threshold Minimum PIP in at least one group to include a taxon.
#' @param font Base font size in points.
#' @param name Optional name token appended to filenames.
#'
#' @return Invisible NULL.
#' @export

Go_mixture_compare <- function(results_list,
                                project,
                                outcomes      = NULL,
                                pip_threshold = 0.1,
                                font          = 10,
                                name          = NULL) {

  grp_names <- names(results_list)
  if (is.null(grp_names) || length(grp_names) < 2)
    stop("[Go_mixture_compare] results_list must be a named list with >= 2 groups.")

  dir_root   <- sprintf("%s_%s", project, format(Sys.Date(), "%y%m%d"))
  pdf_dir    <- file.path(dir_root, "pdf")
  tab_dir    <- file.path(dir_root, "table")
  dir.create(pdf_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

  name_token <- if (!is.null(name) && nzchar(name)) paste0(".", name) else ""
  date_str   <- format(Sys.Date(), "%y%m%d")

  all_outcomes <- unique(unlist(lapply(results_list, names)))
  if (!is.null(outcomes)) all_outcomes <- intersect(outcomes, all_outcomes)

  for (ny in all_outcomes) {
    pip_all <- lapply(grp_names, function(g) {
      pip <- results_list[[g]][[ny]]$bkmr
      if (is.null(pip)) return(NULL)
      pip$Group <- g
      pip
    })
    pip_all <- do.call(rbind, Filter(Negate(is.null), pip_all))
    if (is.null(pip_all) || nrow(pip_all) == 0) next

    # keep taxa where at least one group exceeds threshold
    vars_show <- pip_all |>
      dplyr::group_by(Variable) |>
      dplyr::summarise(max_pip = max(PIP), .groups = "drop") |>
      dplyr::filter(max_pip >= pip_threshold) |>
      dplyr::pull(Variable)

    if (length(vars_show) == 0) {
      message(sprintf("[Go_mixture_compare] No taxa with PIP >= %.2f for %s",
                      pip_threshold, ny))
      next
    }

    pip_wide <- pip_all |>
      dplyr::filter(Variable %in% vars_show) |>
      dplyr::select(Variable, Label, Group, PIP) |>
      tidyr::pivot_wider(names_from = Group, values_from = PIP,
                         values_fill = 0)

    # plain label for text repel (strip HTML tags)
    pip_wide$PlainLabel <- gsub("<[^>]+>", "", pip_wide$Label)

    g1 <- grp_names[1]; g2 <- grp_names[2]
    pip_wide$MaxPIP <- pmax(pip_wide[[g1]], pip_wide[[g2]])
    pip_wide$Sig    <- ifelse(pip_wide[[g1]] >= 0.5 | pip_wide[[g2]] >= 0.5,
                              "PIP ≥ 0.5", "PIP < 0.5")

    p_cmp <- ggplot2::ggplot(pip_wide,
               ggplot2::aes(x = .data[[g1]], y = .data[[g2]])) +
      ggplot2::geom_abline(slope = 1, intercept = 0,
                           linetype = "dashed", color = "grey60") +
      ggplot2::geom_vline(xintercept = 0.5, linetype = "dotted", color = "grey70") +
      ggplot2::geom_hline(yintercept = 0.5, linetype = "dotted", color = "grey70") +
      ggplot2::geom_point(
        ggplot2::aes(color = MaxPIP, size = MaxPIP, shape = Sig)
      ) +
      ggrepel::geom_text_repel(
        ggplot2::aes(label = PlainLabel),
        size         = font * 0.28,
        max.overlaps = 20,
        segment.color = "grey70"
      ) +
      ggplot2::scale_color_gradient(low = "grey75", high = "#E31A1C",
                                    limits = c(0, 1), name = "Max PIP") +
      ggplot2::scale_size_continuous(range = c(1.5, 5),
                                     limits = c(0, 1), guide = "none") +
      ggplot2::scale_shape_manual(values = c("PIP ≥ 0.5" = 16,
                                             "PIP < 0.5"  = 1),
                                  name = NULL) +
      ggplot2::coord_fixed(xlim = c(0, 1), ylim = c(0, 1)) +
      ggplot2::theme_bw(base_size = font) +
      ggplot2::labs(
        title    = sprintf("BKMR PIP: %s vs %s  →  %s", g1, g2, ny),
        subtitle = sprintf("Above diagonal = more important in %s; below = %s", g2, g1),
        x        = sprintf("PIP (%s)", g1),
        y        = sprintf("PIP (%s)", g2)
      )

    fname_c <- file.path(pdf_dir,
                         sprintf("%s.mixture.pip_compare_%s%s.%s.pdf",
                                 project, ny, name_token, date_str))
    pdf(fname_c, width = 5.5, height = 5.5)
    print(p_cmp)
    dev.off()
    message(sprintf("[Go_mixture_compare] %s", basename(fname_c)))

    write.csv(pip_wide[, c("Variable","PlainLabel", g1, g2, "MaxPIP")],
              file.path(tab_dir,
                        sprintf("%s.mixture.pip_compare_%s%s.%s.csv",
                                project, ny, name_token, date_str)),
              row.names = FALSE, quote = FALSE)
  }
  invisible(NULL)
}
