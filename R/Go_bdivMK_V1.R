#' Beta-Diversity Ordination with MiRKAT / MiRKAT-LMM
#'
#' Identical structure to Go_bdivPM() — PCoA plots with kernel-based association
#' test (MiRKAT or MiRKAT-LMM) annotated in the plot corner and saved to CSV.
#'
#' @param psIN          Phyloseq object.
#' @param cate.vars     Grouping variable(s) to test.
#' @param project       Project name (used in filenames).
#' @param orders        Factor level order vector.
#' @param distance_metrics e.g. c("bray","unifrac","wunifrac").
#' @param cate.conf     Covariate column name(s) (optional).
#' @param plot          Ordination method, e.g. "PCoA".
#' @param ellipse       TRUE / FALSE / grouping column for ellipse color.
#' @param statistics    Logical. Compute and annotate MiRKAT p-values.
#' @param mycols        Manual color palette.
#' @param paired        Column name for paired-line overlay (optional).
#' @param combination   NULL or integer (2/3) for pairwise combinations.
#' @param shapes        Column name for point shape mapping (optional).
#' @param ID            Column name for text labels (optional).
#' @param facet         Facet column name (optional).
#' @param name          Extra filename suffix (optional).
#' @param addnumber     Add (n=x) to group legend labels.
#' @param height,width  PDF dimensions.
#' @param plotCols,plotRows  Multiplot grid dimensions.
#' @param strata_var    Subject/block ID column for MiRKAT-LMM (repeated measures).
#'   - NULL or confounded design  → MiRKAT() [subtitle shows fallback notice]
#'   - Valid repeated-measures ID → MiRKAT_LMM()
#'   Same logic as strata_var in Go_bdivPM().
#'
#' @return PDF(s) + CSV tables written to {project}_{date}/pdf/ and table/mirkat/.
#' @export

Go_bdivMK <- function(psIN, cate.vars, project, orders, distance_metrics,
                      cate.conf   = NULL,
                      plot        = "PCoA",
                      ellipse     = TRUE,
                      statistics  = TRUE,
                      mycols      = NULL,
                      paired      = NULL,
                      combination = NULL,
                      shapes      = NULL,
                      ID          = NULL,
                      facet       = NULL,
                      name        = NULL,
                      addnumber   = TRUE,
                      height, width,
                      plotCols    = 2, plotRows = 1,
                      strata_var  = NULL) {

  if (!is.null(dev.list())) dev.off()

  # ── out dirs ──────────────────────────────────────────────────────────────
  out        <- file.path(sprintf("%s_%s", project, format(Sys.Date(), "%y%m%d")))
  out_path   <- file.path(out, "pdf")
  out_table  <- file.path(out, "table")
  out_dist   <- file.path(out, "table", "dist")
  out_mirkat <- file.path(out, "table", "mirkat")
  for (d in c(out, out_path, out_table, out_dist, out_mirkat))
    if (!dir.exists(d)) dir.create(d, recursive = TRUE)

  # ── subtitle builder (identical to Go_bdivPM) ─────────────────────────────
  build_plot_subtitle <- function(cate.conf, strata_var, strata_confounded = FALSE) {
    parts <- character(0)
    if (!is.null(cate.conf) && length(cate.conf) > 0)
      parts <- c(parts, sprintf("covariates=%s", paste(cate.conf, collapse = "+")))
    if (!is.null(strata_var) && length(strata_var) > 0) {
      lbl <- if (strata_confounded)
        sprintf("strata=%s [fallback: unconstrained MiRKAT]", paste(strata_var, collapse = "+"))
      else
        sprintf("strata=%s", paste(strata_var, collapse = "+"))
      parts <- c(parts, lbl)
    }
    if (length(parts) == 0) return(NULL)
    paste(parts, collapse = "\n")
  }

  # ── shared helpers (identical to Go_bdivPM) ───────────────────────────────
  .safe_levels <- function(x, levs) {
    x <- factor(x); if (is.null(levs)) return(x)
    factor(x, levels = levs[levs %in% levels(x)])
  }
  .shape_levels <- function(x, levs = NULL) {
    x <- factor(x); if (nlevels(x) == 0 || is.null(levs)) return(x)
    keep <- levs[levs %in% levels(x)]
    if (length(keep) == 0) return(x)
    factor(x, levels = keep)
  }
  .shape_values <- function(n) {
    base_shapes <- c(1,16,8,0,15,2,17,11,10,12,3,4,5,6,7,9,13,14,18,19,20,21,22,23,24,25)
    if (n > length(base_shapes))
      warning(sprintf("Shape variable has %d levels but only %d shapes supported; shapes will repeat.",
                      n, length(base_shapes)))
    rep(base_shapes, length.out = n)
  }
  .has_any    <- function(x, key) !is.null(x) && length(x) >= 1 && any(x %in% key)
  .adj_label  <- function(method) if (tolower(as.character(method)) == "bh") "FDR(BH)" else sprintf("Adj(%s)", method)
  .axis_percent <- function(ordi_obj) {
    rel <- try(ordi_obj$values$Relative_eig, silent = TRUE)
    if (!inherits(rel, "try-error") && !is.null(rel) && length(rel) >= 2)
      return(as.numeric(rel[1:2]) * 100)
    eig <- try(ordi_obj$values$Eigenvalues, silent = TRUE)
    if (inherits(eig, "try-error") || is.null(eig)) return(c(NA_real_, NA_real_))
    eig <- as.numeric(eig)
    if (length(eig) < 2 || all(is.na(eig))) return(c(NA_real_, NA_real_))
    den <- sum(eig[eig > 0], na.rm = TRUE)
    if (!is.finite(den) || den <= 0) den <- sum(abs(eig), na.rm = TRUE)
    if (!is.finite(den) || den <= 0) return(c(NA_real_, NA_real_))
    as.numeric(pmax(eig[1:2], 0) / den * 100)
  }

  if (!is.null(facet) && length(facet) >= 1 && isTRUE(addnumber)) {
    addnumber <- FALSE
    message("facet detected: addnumber forced to FALSE")
  }

  # ── strata confounding check (identical to Go_bdivPM) ─────────────────────
  .is_strata_confounded <- function(id_vec, mvar_vec) {
    if (is.null(id_vec) || is.null(mvar_vec)) return(FALSE)
    id_vec <- as.character(id_vec); mvar_vec <- as.character(mvar_vec)
    if (length(unique(na.omit(id_vec))) < 2) return(FALSE)
    all(tapply(mvar_vec, id_vec, function(x) length(unique(na.omit(x)))) <= 1)
  }

  # ── MiRKAT-specific helpers ───────────────────────────────────────────────

  # Build numeric outcome from factor
  # binary → 0/1; 3+ levels → warn and use numeric codes (combination mode recommended)
  .mk_y <- function(vec, mvar_name = "") {
    f <- factor(vec)
    if (nlevels(f) < 2) return(NULL)
    if (nlevels(f) > 2) {
      warning(sprintf(
        "mvar '%s' has %d levels. MiRKAT with out_type='C' treats numeric codes ",
        mvar_name, nlevels(f)
      ), "as a continuous proxy — unordered factor coding is arbitrary. ",
      "Use combination=2 for pairwise comparisons (recommended).")
    }
    as.numeric(f) - 1L
  }

  # Build covariate matrix — preserves numeric columns as-is
  .mk_X <- function(map_df, cate.conf, mvar) {
    if (is.null(cate.conf) || length(cate.conf) == 0) return(NULL)
    cvars <- intersect(setdiff(cate.conf, c("SampleType", mvar)), names(map_df))
    if (length(cvars) == 0) return(NULL)
    do.call(cbind, lapply(cvars, function(cv) {
      x <- map_df[[cv]]
      if (is.numeric(x)) x else as.numeric(factor(x))
    }))
  }

  # Run MiRKAT or MiRKAT_LMM; returns list(pval, krv, method_tag)
  .run_mirkat <- function(map_df, mvar, dist_mat, strata_var, cate.conf) {
    y <- .mk_y(map_df[[mvar]], mvar_name = mvar)
    if (is.null(y)) return(list(pval = NA_real_, krv = NA_real_, method_tag = "MiRKAT"))

    K <- try(MiRKAT::D2K(as.matrix(dist_mat)), silent = TRUE)
    if (inherits(K, "try-error")) return(list(pval = NA_real_, krv = NA_real_, method_tag = "MiRKAT"))

    X_mat    <- .mk_X(map_df, cate.conf, mvar)
    out_type <- if (length(unique(y)) == 2) "D" else "C"

    use_lmm <- !is.null(strata_var) && strata_var %in% names(map_df) &&
      !.is_strata_confounded(map_df[[strata_var]], map_df[[mvar]])

    # MiRKAT_LMM availability check
    if (use_lmm && !("MiRKAT_LMM" %in% getNamespaceExports("MiRKAT"))) {
      warning("MiRKAT_LMM not available in your MiRKAT version; falling back to MiRKAT.")
      use_lmm <- FALSE
    }

    if (use_lmm) {
      tt <- try(
        res <- MiRKAT::MiRKAT_LMM(
          y      = y,
          X      = X_mat,
          id     = as.character(map_df[[strata_var]]),
          Ks     = list(K = K),
          method = "davies"
        ), silent = TRUE
      )
      if (inherits(tt, "try-error")) return(list(pval = NA_real_, krv = NA_real_, method_tag = "MiRKAT-LMM"))
      pval <- tryCatch(as.numeric(res$p_values[1]), error = function(e) NA_real_)
      return(list(pval = pval, krv = NA_real_, method_tag = "MiRKAT-LMM"))
    } else {
      tt <- try(
        res <- MiRKAT::MiRKAT(
          y          = y,
          Ks         = list(K = K),
          X          = X_mat,
          out_type   = out_type,
          method     = "davies",
          omnibus    = "permutation",
          returnKRV  = TRUE,
          returnR2   = FALSE
        ), silent = TRUE
      )
      if (inherits(tt, "try-error")) return(list(pval = NA_real_, krv = NA_real_, method_tag = "MiRKAT"))
      pval <- tryCatch(as.numeric(res$p_values[1]), error = function(e) NA_real_)
      krv  <- tryCatch(as.numeric(res$KRV[1]),      error = function(e) NA_real_)
      return(list(pval = pval, krv = krv, method_tag = "MiRKAT"))
    }
  }

  # Variables used for complete.cases (mvar + cate.conf + strata_var)
  .mirkat_model_vars <- function(mvar, cate.conf, strata_var = NULL) {
    conf_vars <- if (!is.null(cate.conf) && length(cate.conf) > 0)
      setdiff(cate.conf, "SampleType") else character(0)
    unique(c(mvar, conf_vars, strata_var))
  }

  # Annotation label builder
  .mk_ann_label <- function(distance_metric, mk_res) {
    pstr <- formatC(mk_res$pval, format = "f", digits = 3)
    if (mk_res$method_tag == "MiRKAT-LMM" || is.na(mk_res$krv)) {
      sprintf("%-12s\n%-18s", distance_metric,
              paste0(mk_res$method_tag, " p=", pstr))
    } else {
      sprintf("%-12s\n%-12s\n%-18s", distance_metric,
              paste0("KRV=", formatC(mk_res$krv, format = "f", digits = 3)),
              paste0(mk_res$method_tag, " p=", pstr))
    }
  }

  # ── multiplot (identical to Go_bdivPM) ────────────────────────────────────
  multiplot <- function(..., plotlist = NULL, file, cols = 1, rows = 1) {
    require(grid)
    plots    <- c(list(...), plotlist)
    numPlots <- length(plots)
    i <- 1
    while (i < numPlots + 1) {
      numToPlot <- min(numPlots - i + 1, cols * rows)
      layout    <- matrix(seq(i, i + cols * rows - 1), ncol = cols, nrow = rows, byrow = TRUE)
      if (numToPlot == 1) {
        print(plots[[i]])
      } else {
        grid.newpage()
        pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
        for (j in i:(i + numToPlot - 1)) {
          matchidx <- as.data.frame(which(layout == j, arr.ind = TRUE))
          print(plots[[j]], vp = viewport(layout.pos.row = matchidx$row,
                                          layout.pos.col = matchidx$col))
        }
      }
      i <- i + numToPlot
    }
  }

  # ── facet-level MiRKAT annotation (mirrors .perm_ann_by_facet) ───────────
  .mirkat_ann_by_facet <- function(ps_obj, pdataframe, mvar, facet,
                                   distance_metric, cate.conf, project,
                                   name, strata_var = NULL, comparison = NULL) {
    distance_metric <- as.character(distance_metric)
    levs <- if (is.factor(pdataframe[[facet]]))
      levels(droplevels(pdataframe[[facet]])) else unique(pdataframe[[facet]])

    stat_list <- lapply(levs, function(flv) {
      map_sub <- as.data.frame(phyloseq::sample_data(ps_obj))
      map_sub <- map_sub[map_sub[[facet]] == flv, , drop = FALSE]
      if (nrow(map_sub) < 3 || length(unique(map_sub[[mvar]])) < 2)
        return(data.frame(facet_val = flv, pval = NA_real_, krv = NA_real_,
                          method_tag = "MiRKAT", stringsAsFactors = FALSE))

      ps_sub <- phyloseq::prune_samples(rownames(map_sub), ps_obj)
      map_sub2 <- as.data.frame(phyloseq::sample_data(ps_sub))

      model_vars <- .mirkat_model_vars(mvar, cate.conf, strata_var)
      model_vars <- intersect(model_vars, names(map_sub2))
      if (length(model_vars) > 0) {
        keep     <- stats::complete.cases(map_sub2[, model_vars, drop = FALSE])
        map_sub2 <- map_sub2[keep, , drop = FALSE]
      }
      if (nrow(map_sub2) < 3 || length(unique(map_sub2[[mvar]])) < 2)
        return(data.frame(facet_val = flv, pval = NA_real_, krv = NA_real_,
                          method_tag = "MiRKAT", stringsAsFactors = FALSE))

      ps_sub    <- phyloseq::prune_samples(rownames(map_sub2), ps_sub)
      dist_list <- Go_dist(psIN = ps_sub, project = project, name = name,
                           cate.vars = mvar, distance_metrics = distance_metric)
      dist_mat  <- dist_list[[distance_metric]]

      mk_res <- .run_mirkat(map_sub2, mvar, dist_mat, strata_var, cate.conf)

      # save facet CSV
      cbn_tag <- if (!is.null(comparison)) paste0(".", paste(comparison, collapse = "_vs_")) else ""
      fn <- file.path(out_mirkat,
                      sprintf("%s.%s.%s.%s%s.facet=%s%s.csv",
                              mk_res$method_tag, distance_metric, project, mvar,
                              cbn_tag, flv,
                              ifelse(is.null(strata_var), "", paste0(".strata=", strata_var))))
      utils::write.csv(data.frame(facet = flv, metric = distance_metric,
                                  mvar = mvar, pval = mk_res$pval,
                                  krv = mk_res$krv, method = mk_res$method_tag),
                       fn, row.names = FALSE)

      data.frame(facet_val = flv, pval = mk_res$pval, krv = mk_res$krv,
                 method_tag = mk_res$method_tag, stringsAsFactors = FALSE)
    })

    stat_df <- dplyr::bind_rows(stat_list)
    ann_df  <- dplyr::mutate(stat_df,
      label = mapply(.mk_ann_label, distance_metric, split(stat_df, seq_len(nrow(stat_df)))),
      x = -Inf, y = -Inf
    )
    ann_df[[facet]] <- ann_df$facet_val
    ann_df[[facet]] <- factor(ann_df[[facet]], levels = levels(pdataframe[[facet]]))
    ann_df
  }

  # ── sanity ────────────────────────────────────────────────────────────────
  tt <- try(mycols, TRUE); if (inherits(tt, "try-error")) mycols <- NULL
  tt <- try(orders, TRUE); if (inherits(tt, "try-error")) orders <- NULL

  # ══════════════════════════════════════════════════════════════════════════
  # MAIN LOOP  (identical structure to Go_bdivPM)
  # ══════════════════════════════════════════════════════════════════════════
  for (mvar in cate.vars) {
    mapping <- data.frame(sample_data(psIN))
    mapping[, mvar] <- factor(mapping[, mvar])
    sample_data(psIN) <- mapping

    if (.has_any(facet, mvar))  next
    if (.has_any(shapes, mvar)) next

    # ────────────────────────────────────────────────────────────────────────
    # COMBINATION path
    # ────────────────────────────────────────────────────────────────────────
    if (!is.null(combination)) {
      mapping[, mvar]   <- .safe_levels(mapping[, mvar], orders)
      group.cbn         <- combn(x = levels(mapping[, mvar]), m = combination)
      group_comparisons <- lapply(seq_len(ncol(group.cbn)), function(i) group.cbn[, i])

      pdf(sprintf("%s/ordi.%s.MK.%s.%s.%s%s%s%s%s%s%s%s%s.pdf", out_path,
                  plot, paste(distance_metrics, collapse = "+"), project, mvar,
                  ifelse(is.null(facet),       "", paste0(facet, ".")),
                  paste0("(cbn=", combination, ")."),
                  ifelse(is.null(cate.conf),   "", "with_confounder."),
                  ifelse(is.null(paired),       "", paste0("(paired=", paired, ").")),
                  ifelse(is.null(name),         "", paste0(name, ".")),
                  ifelse(isFALSE(ellipse), "ellipse_FALSE.",
                         ifelse(isTRUE(ellipse), "", paste0("ellipse_", ellipse, "."))),
                  ifelse(is.null(strata_var),   "", paste0("(strata=", strata_var, ").")),
                  format(Sys.Date(), "%y%m%d")),
          height = height, width = width)
      plotlist <- list()

      for (i in seq_along(group_comparisons)) {
        group.combination <- unlist(group_comparisons[i])
        mapping.cbn <- subset(mapping, mapping[, mvar] %in% group.combination[seq_len(combination)])
        psIN.cbn <- psIN; sample_data(psIN.cbn) <- mapping.cbn

        for (distance_metric in distance_metrics) {
          # remove NA
          mapping.sel         <- data.frame(sample_data(psIN.cbn))
          psIN.cbn.na         <- prune_samples(rownames(mapping.sel[!is.na(mapping.sel[, mvar]), ]), psIN.cbn)
          mapping.sel.na.rem  <- data.frame(sample_data(psIN.cbn.na))
          mapping.sel.na.rem[, mvar] <- factor(mapping.sel.na.rem[, mvar])

          # ordination
          plist <- plyr::llply(as.list(plot), function(i, psIN.na, distance_metric) {
            ordi <- ordinate(psIN.na, method = i, distance = distance_metric)
            df   <- as.data.frame(ordi$vectors[, 1:2]); colnames(df) <- c("Axis_1", "Axis_2")
            pct  <- .axis_percent(ordi)
            df$Axis1_Percent <- pct[1]; df$Axis2_Percent <- pct[2]
            cbind(df, as.data.frame(sample_data(psIN.na)))
          }, psIN.cbn.na, distance_metric)
          names(plist) <- plot
          pdataframe   <- plyr::ldply(plist, identity); names(pdataframe)[1] <- "method"

          if (!is.null(facet) && facet %in% names(pdataframe))
            pdataframe[, facet] <- factor(pdataframe[, facet],
                                          levels = intersect(orders, unique(pdataframe[, facet])))
          pdataframe[, mvar] <- factor(pdataframe[, mvar],
                                       levels = intersect(orders, unique(pdataframe[, mvar])))

          if (isTRUE(addnumber)) {
            for (nm in unique(pdataframe[, mvar])) {
              tot <- sum(pdataframe[, mvar] == nm, na.rm = TRUE)
              levels(pdataframe[[mvar]])[levels(pdataframe[[mvar]]) == nm] <- paste0(nm, " (n=", tot, ")")
            }
          }

          axis1_pct <- mean(pdataframe$Axis1_Percent, na.rm = TRUE)
          axis2_pct <- mean(pdataframe$Axis2_Percent, na.rm = TRUE)

          # strata status for subtitle
          strata_confounded <- !is.null(strata_var) &&
            strata_var %in% names(mapping.sel.na.rem) &&
            .is_strata_confounded(mapping.sel.na.rem[[strata_var]], mapping.sel.na.rem[[mvar]])
          subtitle_text <- build_plot_subtitle(cate.conf, strata_var, strata_confounded)

          if (!is.null(shapes) && shapes %in% names(pdataframe)) {
            pdataframe[, shapes] <- .shape_levels(pdataframe[, shapes], orders)
            p <- ggplot(pdataframe, aes_string(x = "Axis_1", y = "Axis_2", color = mvar)) +
              geom_point(aes_string(shape = shapes), size = 0.9, alpha = 1) +
              scale_shape_manual(values = .shape_values(nlevels(pdataframe[, shapes])))
          } else {
            p <- ggplot(pdataframe, aes_string(x = "Axis_1", y = "Axis_2", color = mvar)) +
              geom_point(size = 0.9, alpha = 1)
          }

          p <- p + labs(x       = paste0("Axis 1 (", sprintf("%.2f", axis1_pct), "%)"),
                        y       = paste0("Axis 2 (", sprintf("%.2f", axis2_pct), "%)"),
                        title   = sprintf("%s (%s)", mvar, distance_metric),
                        subtitle = subtitle_text) +
            facet_wrap(~ method, scales = "free") + theme_bw() +
            theme(strip.background  = element_blank(),
                  legend.position   = "bottom", legend.title = element_blank(),
                  legend.justification = "left", legend.box = "vertical",
                  legend.box.margin = ggplot2::margin(0,0,0,-1,"cm"),
                  legend.spacing.y  = ggplot2::unit(0.02, "cm"),
                  legend.key.height = ggplot2::unit(0.25, "cm"),
                  legend.key.width  = ggplot2::unit(0.35, "cm"),
                  legend.margin     = ggplot2::margin(0,0,0,0,"cm"),
                  plot.subtitle     = element_text(size = 6, lineheight = 0.9))

          if (!is.null(mycols)) p <- p + scale_color_manual(values = mycols)
          p <- p + guides(color = ggplot2::guide_legend(ncol = 1, byrow = TRUE),
                          shape = ggplot2::guide_legend(ncol = 1, byrow = TRUE))
          if (!is.null(ID) && ID %in% names(pdataframe))
            p <- p + ggrepel::geom_text_repel(aes_string(label = ID), size = 2)
          if (isTRUE(ellipse))    p <- p + stat_ellipse(type = "norm", linetype = 2)
          else if (!is.null(ellipse) && !isFALSE(ellipse))
            p <- p + stat_ellipse(aes_string(group = ellipse, color = ellipse), type = "norm", linetype = 2)

          # ── MiRKAT statistics ──────────────────────────────────────────
          if (isTRUE(statistics)) {
            cbn_tag <- paste(group.combination, collapse = "_vs_")
            if (!is.null(facet) && facet %in% names(pdataframe)) {
              ann_df <- .mirkat_ann_by_facet(ps_obj = psIN.cbn.na, pdataframe = pdataframe,
                                             mvar = mvar, facet = facet,
                                             distance_metric = distance_metric,
                                             cate.conf = cate.conf, project = project,
                                             name = name, strata_var = strata_var,
                                             comparison = group.combination)
              p <- p + ggplot2::geom_text(data = ann_df,
                                          mapping = aes(x = x, y = y, label = label),
                                          size = 3, hjust = -0.005, vjust = -0.3,
                                          lineheight = 0.95, inherit.aes = FALSE)
            } else {
              set.seed(1)
              dist_obj  <- Go_dist(psIN = psIN.cbn.na, project = project, cate.vars = mvar,
                                   name = name, distance_metrics = distance_metric)
              x         <- as.dist(dist_obj[[distance_metric]])
              factors   <- mapping.sel.na.rem[, mvar]
              map.pair  <- subset(mapping.sel.na.rem, mapping.sel.na.rem[, mvar] %in% unique(factors))

              model_vars <- .mirkat_model_vars(mvar, cate.conf, strata_var)
              model_vars <- intersect(model_vars, names(map.pair))
              if (length(model_vars) > 0) {
                keep     <- stats::complete.cases(map.pair[, model_vars, drop = FALSE])
                map.pair <- map.pair[keep, , drop = FALSE]
              }
              if (nrow(map.pair) < 3 || length(unique(map.pair[, mvar])) < 2) next

              x1 <- as.matrix(x)[rownames(map.pair), rownames(map.pair), drop = FALSE]

              mk_res <- .run_mirkat(map.pair, mvar, x1, strata_var, cate.conf)

              # save CSV — combination info in filename to prevent overwrite
              fn <- file.path(out_mirkat,
                              sprintf("%s.%s.%s.%s.%s%s.csv",
                                      mk_res$method_tag, distance_metric, project, mvar,
                                      cbn_tag,
                                      ifelse(is.null(strata_var), "", paste0(".strata=", strata_var))))
              utils::write.csv(data.frame(comparison = cbn_tag, metric = distance_metric,
                                          mvar = mvar, pval = mk_res$pval, krv = mk_res$krv,
                                          method = mk_res$method_tag),
                               fn, row.names = FALSE)

              ann_df <- data.frame(x = -Inf, y = -Inf,
                                   label = .mk_ann_label(distance_metric, mk_res))
              p <- p + ggplot2::geom_text(data = ann_df, aes(x = x, y = y, label = label),
                                          size = 3, hjust = -0.005, vjust = -0.3,
                                          lineheight = 0.95, inherit.aes = FALSE)
            }
          }
          # ──────────────────────────────────────────────────────────────

          if (!is.null(facet) && facet %in% names(pdataframe))
            p <- p + facet_wrap(as.formula(sprintf("~ %s", facet)), scales = "free_x", ncol = 1)
          if (!is.null(paired) && paired %in% names(pdataframe))
            p <- p + geom_line(aes_string(group = paired), color = "grey", linewidth = 0.2,
                               arrow = arrow(type = "closed", length = unit(0.025, "inches")))
          p <- p +
            scale_x_continuous(expand = expansion(mult = c(0.03, 0.03))) +
            scale_y_continuous(expand = expansion(mult = c(0.06, 0.03))) +
            theme(panel.grid = element_blank(), legend.key = element_blank(),
                  panel.background = element_rect(fill = "white", colour = "Black",
                                                  linewidth = 0.7, linetype = "solid"),
                  aspect.ratio = 1,
                  plot.title = element_text(size = 8, face = "bold")) +
            geom_vline(xintercept = 0, linewidth = 0.1) +
            geom_hline(yintercept = 0, linewidth = 0.1)

          plotlist[[length(plotlist) + 1]] <- p
        }
      }
      multiplot(plotlist = plotlist, cols = plotCols, rows = plotRows)
      dev.off()

    # ────────────────────────────────────────────────────────────────────────
    # NO-COMBINATION path
    # ────────────────────────────────────────────────────────────────────────
    } else {
      pdf(sprintf("%s/ordi.%s.MK.%s.%s.%s%s%s%s%s%s%s%s.pdf", out_path,
                  plot, paste(distance_metrics, collapse = "+"), project, mvar,
                  ifelse(is.null(facet),       "", paste0(facet, ".")),
                  ifelse(is.null(cate.conf),   "", "with_confounder."),
                  ifelse(is.null(paired),       "", paste0("(paired=", paired, ").")),
                  ifelse(is.null(name),         "", paste0(name, ".")),
                  ifelse(isFALSE(ellipse), "ellipse_FALSE.",
                         ifelse(isTRUE(ellipse), "", paste0("ellipse_", ellipse, "."))),
                  ifelse(is.null(strata_var),   "", paste0("(strata=", strata_var, ").")),
                  format(Sys.Date(), "%y%m%d")),
          height = height, width = width)
      plotlist <- list()

      for (distance_metric in distance_metrics) {
        mapping.sel         <- data.frame(sample_data(psIN))
        psIN.na             <- prune_samples(rownames(mapping.sel[!is.na(mapping.sel[, mvar]), ]), psIN)
        mapping.sel.na.rem  <- data.frame(sample_data(psIN.na))
        mapping.sel.na.rem[, mvar] <- factor(mapping.sel.na.rem[, mvar])

        # ordination
        plist <- plyr::llply(as.list(plot), function(i, psIN.na, distance_metric) {
          ordi <- ordinate(psIN.na, method = i, distance = distance_metric)
          df   <- as.data.frame(ordi$vectors[, 1:2]); colnames(df) <- c("Axis_1", "Axis_2")
          pct  <- .axis_percent(ordi)
          df$Axis1_Percent <- pct[1]; df$Axis2_Percent <- pct[2]
          cbind(df, as.data.frame(sample_data(psIN.na)))
        }, psIN.na, distance_metric)
        names(plist) <- plot
        pdataframe   <- plyr::ldply(plist, identity); names(pdataframe)[1] <- "method"

        if (!is.null(facet) && facet %in% names(pdataframe))
          pdataframe[, facet] <- factor(pdataframe[, facet],
                                        levels = intersect(orders, unique(pdataframe[, facet])))
        pdataframe[, mvar] <- factor(pdataframe[, mvar],
                                     levels = intersect(orders, unique(pdataframe[, mvar])))

        if (isTRUE(addnumber)) {
          for (nm in unique(pdataframe[, mvar])) {
            tot <- sum(pdataframe[, mvar] == nm, na.rm = TRUE)
            levels(pdataframe[[mvar]])[levels(pdataframe[[mvar]]) == nm] <- paste0(nm, " (n=", tot, ")")
          }
        }

        axis1_pct <- mean(pdataframe$Axis1_Percent, na.rm = TRUE)
        axis2_pct <- mean(pdataframe$Axis2_Percent, na.rm = TRUE)

        # strata status for subtitle
        strata_confounded <- !is.null(strata_var) &&
          strata_var %in% names(mapping.sel.na.rem) &&
          .is_strata_confounded(mapping.sel.na.rem[[strata_var]], mapping.sel.na.rem[[mvar]])
        subtitle_text <- build_plot_subtitle(cate.conf, strata_var, strata_confounded)

        if (!is.null(shapes) && shapes %in% names(pdataframe)) {
          pdataframe[, shapes] <- .shape_levels(pdataframe[, shapes], orders)
          p <- ggplot(pdataframe, aes_string(x = "Axis_1", y = "Axis_2", color = mvar)) +
            geom_point(aes_string(shape = shapes), size = 0.9, alpha = 1) +
            scale_shape_manual(values = .shape_values(nlevels(pdataframe[, shapes])))
        } else {
          p <- ggplot(pdataframe, aes_string(x = "Axis_1", y = "Axis_2", color = mvar)) +
            geom_point(size = 0.9, alpha = 1)
        }

        p <- p + labs(x        = paste0("Axis 1 (", sprintf("%.2f", axis1_pct), "%)"),
                      y        = paste0("Axis 2 (", sprintf("%.2f", axis2_pct), "%)"),
                      title    = sprintf("%s (%s)", mvar, distance_metric),
                      subtitle = subtitle_text) +
          facet_wrap(~ method, scales = "free") + theme_bw() +
          theme(strip.background  = element_blank(),
                legend.position   = "bottom", legend.title = element_blank(),
                legend.justification = "left", legend.box = "vertical",
                legend.box.margin = ggplot2::margin(0,0,0,-1,"cm"),
                legend.spacing.y  = ggplot2::unit(0.02, "cm"),
                legend.key.height = ggplot2::unit(0.25, "cm"),
                legend.key.width  = ggplot2::unit(0.35, "cm"),
                legend.margin     = ggplot2::margin(0,0,0,0,"cm"),
                plot.title        = element_text(size = 8, face = "bold"),
                plot.subtitle     = element_text(size = 6, lineheight = 0.9))

        if (!is.null(mycols)) p <- p + scale_color_manual(values = mycols)
        p <- p + guides(color = ggplot2::guide_legend(ncol = 1, byrow = TRUE),
                        shape = ggplot2::guide_legend(ncol = 1, byrow = TRUE))
        if (!is.null(ID) && ID %in% names(pdataframe))
          p <- p + ggrepel::geom_text_repel(aes_string(label = ID), size = 2)
        if (isTRUE(ellipse))    p <- p + stat_ellipse(type = "norm", linetype = 2)
        else if (!is.null(ellipse) && !isFALSE(ellipse))
          p <- p + stat_ellipse(aes_string(group = ellipse, color = ellipse), type = "norm", linetype = 2)

        # ── MiRKAT statistics ────────────────────────────────────────────
        if (isTRUE(statistics)) {
          if (!is.null(facet) && facet %in% names(pdataframe)) {
            ann_df <- .mirkat_ann_by_facet(ps_obj = psIN.na, pdataframe = pdataframe,
                                           mvar = mvar, facet = facet,
                                           distance_metric = distance_metric,
                                           cate.conf = cate.conf, project = project,
                                           name = name, strata_var = strata_var)
            p <- p + ggplot2::geom_text(data = ann_df,
                                        mapping = aes(x = x, y = y, label = label),
                                        size = 3, hjust = -0.005, vjust = -0.3,
                                        lineheight = 0.95, inherit.aes = FALSE)
          } else {
            set.seed(1)
            dist_obj  <- Go_dist(psIN = psIN.na, project = project, cate.vars = mvar,
                                 name = name, distance_metrics = distance_metric)
            x         <- as.dist(dist_obj[[distance_metric]])
            factors   <- mapping.sel.na.rem[, mvar]
            map.pair  <- subset(mapping.sel.na.rem, mapping.sel.na.rem[, mvar] %in% unique(factors))

            model_vars <- .mirkat_model_vars(mvar, cate.conf, strata_var)
            model_vars <- intersect(model_vars, names(map.pair))
            if (length(model_vars) > 0) {
              keep     <- stats::complete.cases(map.pair[, model_vars, drop = FALSE])
              map.pair <- map.pair[keep, , drop = FALSE]
            }
            if (nrow(map.pair) < 3 || length(unique(map.pair[, mvar])) < 2) next

            x1 <- as.matrix(x)[rownames(map.pair), rownames(map.pair), drop = FALSE]

            mk_res <- .run_mirkat(map.pair, mvar, x1, strata_var, cate.conf)

            # save CSV
            fn <- file.path(out_mirkat,
                            sprintf("%s.%s.%s.%s%s.csv",
                                    mk_res$method_tag, distance_metric, project, mvar,
                                    ifelse(is.null(strata_var), "", paste0(".strata=", strata_var))))
            utils::write.csv(data.frame(metric = distance_metric, mvar = mvar,
                                        pval = mk_res$pval, krv = mk_res$krv,
                                        method = mk_res$method_tag),
                             fn, row.names = FALSE)

            ann_df <- data.frame(x = -Inf, y = -Inf,
                                 label = .mk_ann_label(distance_metric, mk_res))
            p <- p + ggplot2::geom_text(data = ann_df, aes(x = x, y = y, label = label),
                                        size = 3, hjust = -0.005, vjust = -0.3,
                                        lineheight = 0.95, inherit.aes = FALSE)
          }
        }
        # ────────────────────────────────────────────────────────────────

        if (!is.null(facet) && facet %in% names(pdataframe)) {
          mapping.sel.na.rem[[facet]] <- factor(mapping.sel.na.rem[[facet]], levels = orders)
          p <- p + facet_wrap(as.formula(sprintf("~ %s", facet)), scales = "free_x", ncol = 1)
        }
        if (!is.null(paired) && paired %in% names(pdataframe))
          p <- p + geom_line(aes_string(group = paired), color = "grey", linewidth = 0.2,
                             arrow = arrow(type = "closed", length = unit(0.025, "inches")))
        p <- p +
          scale_x_continuous(expand = expansion(mult = c(0.03, 0.03))) +
          scale_y_continuous(expand = expansion(mult = c(0.06, 0.03))) +
          theme(panel.grid = element_blank(), legend.key = element_blank(),
                panel.background = element_rect(fill = "white", colour = "Black",
                                                linewidth = 0.7, linetype = "solid"),
                aspect.ratio = 1) +
          geom_vline(xintercept = 0, linewidth = 0.1) +
          geom_hline(yintercept = 0, linewidth = 0.1)

        plotlist[[length(plotlist) + 1]] <- p
      }
      multiplot(plotlist = plotlist, cols = plotCols, rows = plotRows)
      dev.off()
    }
  }
}
