#' Generate Biodiversity Plots for Microbiome Data
#'
#' @param psIN Phyloseq object containing microbiome data.
#' @param cate.vars Categorical variables to be used for grouping in the plot.
#' @param project Name of the project for labeling.
#' @param orders Ordering of factors in the plot.
#' @param distance_metrics Vector of distance metrics to be used.
#' @param cate.conf Additional confounder variables, if any.
#' @param plot Type of plot to generate, e.g., "PCoA".
#' @param ellipse Logical indicating whether to include ellipses.
#' @param statistics Logical indicating whether to include statistical tests.
#' @param mycols Custom color palette.
#' @param paired Variable for paired analysis, if applicable.
#' @param combination Number of groups to combine for comparison.
#' @param shapes Shape of points in the plot.
#' @param ID Identifier for labeling points.
#' @param facet Faceting variable for the plot.
#' @param name Optional name for the output plot.
#' @param addnumber Logical indicating whether to add sample numbers to groups.
#' @param height Height of the plot.
#' @param width Width of the plot.
#' @param plotCols Number of columns in the multiplot layout.
#' @param plotRows Number of rows in the multiplot layout.
#' @param strata_var (NEW) Column name to use as permutation blocks for PERMANOVA.
#' @param p_adjust P-value adjustment method for PERMANOVA labels (e.g., "BH", "bonferroni").
#'
#' @return PDF(s) and CSV tables on disk.
#' @export
Go_bdivPM <- function(psIN, cate.vars, project, orders, distance_metrics,
                    cate.conf=NULL,
                    plot="PCoA",
                    ellipse=TRUE,
                    statistics = TRUE,
                    mycols=NULL,
                    paired = NULL,
                    combination=NULL,
                    shapes = NULL,
                    ID = NULL,
                    facet=NULL,
                    name=NULL,
                    addnumber=TRUE,
                    height, width,
                    plotCols = 2, plotRows = 1,
                    strata_var = NULL,
                    p_adjust = "BH") {

  if(!is.null(dev.list())) dev.off()

  # out dir
  out <- file.path(sprintf("%s_%s",project, format(Sys.Date(), "%y%m%d")))
  if(!dir.exists(out)) dir.create(out, recursive = TRUE)
  out_path <- file.path(sprintf("%s_%s/pdf",project, format(Sys.Date(), "%y%m%d")))
  if(!dir.exists(out_path)) dir.create(out_path, recursive = TRUE)
  out_table <- file.path(sprintf("%s_%s/table",project, format(Sys.Date(), "%y%m%d")))
  if(!dir.exists(out_table)) dir.create(out_table, recursive = TRUE)
  out_dist <- file.path(sprintf("%s_%s/table/dist",project, format(Sys.Date(), "%y%m%d")))
  if(!dir.exists(out_dist)) dir.create(out_dist, recursive = TRUE)

  # NEW: perm dir (PERMANOVA 결과 저장 경로)  # <<< NEW
  out_perm <- file.path(sprintf("%s_%s/table/perm",project, format(Sys.Date(), "%y%m%d")))
  if(!dir.exists(out_perm)) dir.create(out_perm, recursive = TRUE)

  title_suffix <- ""
  build_plot_subtitle <- function(cate.conf, strata_var) {
    subtitle_parts <- character(0)
    if (!is.null(cate.conf) && length(cate.conf) > 0) {
      subtitle_parts <- c(subtitle_parts, sprintf("covariates=%s", paste(cate.conf, collapse = "+")))
    }
    if (!is.null(strata_var) && length(strata_var) > 0) {
      subtitle_parts <- c(subtitle_parts, sprintf("strata=%s", paste(strata_var, collapse = "+")))
    }
    if (length(subtitle_parts) == 0) {
      return(NULL)
    }
    paste(subtitle_parts, collapse = "\n")
  }

  # ---------- helpers ----------
  .safe_levels <- function(x, levs) {
    x <- factor(x)
    if (is.null(levs)) return(x)
    keep <- levs[levs %in% levels(x)]
    factor(x, levels = keep)
  }
  .has_any <- function(x, key) !is.null(x) && length(x) >= 1 && any(x %in% key)
  .adj_label <- function(method) {
    if (tolower(as.character(method)) == "bh") "FDR(BH)" else sprintf("Adj(%s)", method)
  }
  .axis_percent <- function(ordi_obj) {
    # Prefer relative eigenvalues when available.
    rel <- try(ordi_obj$values$Relative_eig, silent = TRUE)
    if (!inherits(rel, "try-error") && !is.null(rel) && length(rel) >= 2) {
      return(as.numeric(rel[1:2]) * 100)
    }

    eig <- try(ordi_obj$values$Eigenvalues, silent = TRUE)
    if (inherits(eig, "try-error") || is.null(eig)) return(c(NA_real_, NA_real_))
    eig <- as.numeric(eig)
    if (length(eig) < 2 || all(is.na(eig))) return(c(NA_real_, NA_real_))

    # PCoA can include negative eigenvalues; use positive-only denominator.
    den <- sum(eig[eig > 0], na.rm = TRUE)
    if (!is.finite(den) || den <= 0) {
      den <- sum(abs(eig), na.rm = TRUE)
    }
    if (!is.finite(den) || den <= 0) return(c(NA_real_, NA_real_))

    pct <- pmax(eig[1:2], 0) / den * 100
    as.numeric(pct)
  }

  # Facet mode: disable group-count suffix because counts differ by facet panel.
  if (!is.null(facet) && length(facet) >= 1 && isTRUE(addnumber)) {
    addnumber <- FALSE
    message("facet detected: addnumber forced to FALSE")
  }

  # permutation blocks
  .mk_perm_block <- function(id_vec, nperm = 999) {
    ok <- !is.null(id_vec)
    if (ok) {
      id_vec <- as.vector(id_vec)
      ok <- length(na.omit(id_vec)) > 1 && length(unique(na.omit(id_vec))) > 1
    }
    if (!ok) return(nperm)
    ctrl <- permute::how(blocks = id_vec)
    permute::setNperm(ctrl) <- nperm
    ctrl
  }

  .perm_model_vars <- function(mvar, cate.conf, strata_var = NULL) {
    conf_vars <- if (!is.null(cate.conf) && length(cate.conf) > 0) {
      setdiff(cate.conf, "SampleType")
    } else {
      character(0)
    }
    unique(c(mvar, conf_vars, strata_var))
  }

  multiplot <- function(..., plotlist=NULL, file, cols=1, rows=1) {
    require(grid)
    plots <- c(list(...), plotlist)
    numPlots = length(plots)

    i = 1
    while (i < numPlots + 1) {
      numToPlot <- min(numPlots-i+1, cols*rows)
      layout <- matrix(seq(i, i+cols*rows-1), ncol = cols, nrow = rows, byrow = TRUE)
      if (numToPlot == 1) {
        print(plots[[i]])
      } else {
        grid.newpage()
        pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
        for (j in i:(i+numToPlot-1)) {
          matchidx <- as.data.frame(which(layout == j, arr.ind = TRUE))
          print(plots[[j]], vp = viewport(layout.pos.row = matchidx$row,
                                          layout.pos.col = matchidx$col))
        }
      }
      i <- i+numToPlot
    }
  }

  # facet 있을 때: facet별 PERMANOVA + 코너 고정 주석
  .perm_ann_by_facet <- function(ps_obj, pdataframe, mvar, facet, distance_metric,
                                 cate.conf, project, name, strata_var = NULL) {
    distance_metric <- as.character(distance_metric)

    levs <- if (is.factor(pdataframe[[facet]])) {
      base::levels(base::droplevels(pdataframe[[facet]]))
    } else {
      base::unique(pdataframe[[facet]])
    }

    stat_list <- base::lapply(levs, function(flv){
      map_sub <- base::as.data.frame(phyloseq::sample_data(ps_obj))
      map_sub <- map_sub[map_sub[[facet]] == flv, , drop = FALSE]
      if (base::nrow(map_sub) < 3 || base::length(base::unique(map_sub[[mvar]])) < 2) {
        return(base::data.frame(facet_val = flv, R2 = NA_real_, p = NA_real_))
      }
      ps_sub <- phyloseq::prune_samples(base::rownames(map_sub), ps_obj)

      map_sub2 <- base::data.frame(phyloseq::sample_data(ps_sub))
      model_vars <- .perm_model_vars(mvar, cate.conf, strata_var)
      model_vars <- intersect(model_vars, names(map_sub2))
      if (length(model_vars) > 0) {
        keep <- stats::complete.cases(map_sub2[, model_vars, drop = FALSE])
        map_sub2 <- map_sub2[keep, , drop = FALSE]
      }
      if (base::nrow(map_sub2) < 3 || base::length(base::unique(map_sub2[[mvar]])) < 2) {
        return(base::data.frame(facet_val = flv, R2 = NA_real_, p = NA_real_))
      }

      ps_sub <- phyloseq::prune_samples(base::rownames(map_sub2), ps_sub)
      dist_list <- Go_dist(psIN = ps_sub, project = project, name = name,
                           cate.vars = mvar, distance_metrics = distance_metric)
      x <- stats::as.dist(dist_list[[distance_metric]])

      perm <- 999
      if (!is.null(strata_var) && strata_var %in% names(map_sub2)) {
        perm <- .mk_perm_block(map_sub2[[strata_var]], nperm = 999)
      }

      if (!base::is.null(cate.conf) && base::length(cate.conf) > 0) {
        for (conf in cate.conf) map_sub2[[conf]] <- base::factor(map_sub2[[conf]])
        form <- stats::as.formula(base::sprintf("x ~ %s + %s", mvar,
                                                base::paste(base::setdiff(cate.conf, "SampleType"), collapse = " + ")))
      } else {
        form <- stats::as.formula(base::sprintf("x ~ %s", mvar))
      }
      set.seed(123)
      ad <- vegan::adonis2(form, data = map_sub2, permutations = perm, by = "terms")

      # SAVE facet PERMANOVA table   # <<< NEW
      fn <- file.path(
        out_perm,
        sprintf("PERMANOVA.%s.%s.%s.facet=%s%s.csv",
                distance_metric, project, mvar, flv,
                ifelse(is.null(strata_var), "", paste0(".strata=", strata_var)))
      )
      utils::write.csv(as.data.frame(ad), fn, row.names = TRUE)

      base::data.frame(
        facet_val = flv,
        R2 = base::suppressWarnings(base::round(ad[1,"R2"], 3)),
        p  = base::suppressWarnings(ad[1,"Pr(>F)"]),
        stringsAsFactors = FALSE
      )
    })

    stat_df <- dplyr::bind_rows(stat_list)

    ann_df <- dplyr::mutate(
      stat_df,
      label = base::sprintf(
        "%-12s\n%-12s\n%-18s",
        distance_metric,
        base::paste0("R2=", base::formatC(R2, format="f", digits=3)),
        base::paste0("PERMANOVA p=", base::formatC(p, format="f", digits=3))
      ),
      x = -Inf,
      y = -Inf
    )

    ann_df[[facet]] <- ann_df$facet_val
    ann_df[[facet]] <- base::factor(ann_df[[facet]], levels = base::levels(pdataframe[[facet]]))
    ann_df
  }
  # -----------------------------

  # sanity
  tt <- try(mycols, TRUE); if(inherits(tt, "try-error")) mycols <- NULL
  tt <- try(orders, TRUE); if(inherits(tt, "try-error")) orders <- NULL

  for (mvar in cate.vars) {
    mapping <- data.frame(sample_data(psIN))
    mapping[,mvar] <- factor(mapping[,mvar])
    sample_data(psIN) <- mapping

    if (.has_any(facet, mvar)) next
    if (.has_any(shapes, mvar)) next

    #------------------------------#
    # for group combination or not #
    #------------------------------#
    if (!is.null(combination)){
      mapping[,mvar] <- .safe_levels(mapping[,mvar], orders)
      group.cbn <- combn(x = levels(mapping[,mvar]), m = combination)
      group_comparisons <- lapply(seq_len(ncol(group.cbn)), function(i) group.cbn[,i])

      ord_meths = plot
      pdf(sprintf("%s/ordi.%s.%s.%s.%s%s%s%s%s%s%s%s%s.pdf", out_path,
                  ord_meths, paste(distance_metrics, collapse = "+"), project, mvar,
                  ifelse(is.null(facet), "", paste(facet, ".", sep = "")),
                  ifelse(is.null(combination), "", paste("(cbn=",combination, ").", sep = "")),
                  ifelse(is.null(cate.conf), "", paste("with_confounder", ".", sep = "")),
                  ifelse(is.null(paired), "", paste("(paired=",paired, ").", sep = "")),
                  ifelse(is.null(name), "", paste(name, ".", sep = "")),
                  ifelse(ellipse == FALSE, "ellipse_FALSE.",
                         ifelse(ellipse == TRUE, "", paste("ellipse_", ellipse, ".", sep = ""))),
                  ifelse(is.null(strata_var), "", paste("(strata=", strata_var, ").", sep = "")),
                  format(Sys.Date(), "%y%m%d")), height = height, width = width)
      plotlist <- list()

      for(i in seq_along(group_comparisons)){
        group.combination <- unlist(group_comparisons[i])

        mapping.cbn <- subset(mapping, mapping[,mvar] %in% group.combination[seq_len(combination)])

        psIN.cbn <- psIN; sample_data(psIN.cbn) <- mapping.cbn

        for(distance_metric in distance_metrics){
          # remove NA
          mapping.sel <- data.frame(sample_data(psIN.cbn))
          mapping.sel.na <- mapping.sel[!is.na(mapping.sel[,mvar]), ]
          psIN.cbn.na <- prune_samples(rownames(mapping.sel[!is.na(mapping.sel[,mvar]), ]), psIN.cbn)
          mapping.sel.na.rem <- data.frame(sample_data(psIN.cbn.na ))
          mapping.sel.na.rem[,mvar] <- factor(mapping.sel.na.rem[,mvar])

          # global perm for this subset
          perm_global <- 999
          if (!is.null(strata_var) && strata_var %in% names(mapping.sel.na.rem)) {
            perm_global <- .mk_perm_block(mapping.sel.na.rem[[strata_var]], nperm = 999)
          }

          # ordination
          ord_meths = plot
          plist = plyr::llply(as.list(ord_meths), function(i, psIN.na, distance_metric){
            ordi = ordinate(psIN.na, method=i, distance=distance_metric)
            df = as.data.frame(ordi$vectors[, 1:2]); colnames(df) = c("Axis_1", "Axis_2")
            var_explained <- .axis_percent(ordi)
            df$Axis1_Percent = var_explained[1]; df$Axis2_Percent = var_explained[2]
            metadata = as.data.frame(sample_data(psIN.na))
            cbind(df, metadata)
          }, psIN.cbn.na, distance_metric)

          names(plist) <- ord_meths
          pdataframe = plyr::ldply(plist, identity); names(pdataframe)[1] = "method"

          if (!is.null(facet) && facet %in% names(pdataframe)) {
            pdataframe[, facet] <- factor(pdataframe[, facet], levels = intersect(orders, unique(pdataframe[, facet])))
          }

          pdataframe[, mvar] <- factor(pdataframe[, mvar], levels = intersect(orders, unique(pdataframe[, mvar])))




          if(addnumber==TRUE){
            for (Name in unique(pdataframe[,mvar])) {
              total <- sum(pdataframe[,mvar] == Name, na.rm=TRUE)
              levels(pdataframe[[mvar]])[levels(pdataframe[[mvar]])== Name] <- paste0(Name, " (n=", total, ")")
            }
          }

          axis1_percent_avg <- mean(pdataframe$Axis1_Percent, na.rm = TRUE)
          axis2_percent_avg <- mean(pdataframe$Axis2_Percent, na.rm = TRUE)

          if (!is.null(shapes) && shapes %in% names(pdataframe)) {
            pdataframe[,shapes] <- .safe_levels(pdataframe[,shapes], orders)
            p = ggplot(pdataframe, aes_string(x="Axis_1", y="Axis_2", color=mvar)) +
              geom_point(aes_string(shape=shapes), size=0.9, alpha=1) +
              scale_shape_manual(values = c(1,16,8,0,15,2,17,11,10,12,3,4,5,6,7,8,9,13,14))
          } else {
            p = ggplot(pdataframe, aes_string(x="Axis_1", y="Axis_2", color=mvar)) +
              geom_point(size=0.9, alpha=1)
          }

          subtitle_text <- build_plot_subtitle(cate.conf, strata_var)
          p = p +
            labs(x = paste0("Axis 1 (", sprintf("%.2f", axis1_percent_avg),"%)"),
                 y = paste0("Axis 2 (", sprintf("%.2f", axis2_percent_avg),"%)"),
                 title = sprintf("%s (%s)%s", mvar, distance_metric, title_suffix),
                 subtitle = subtitle_text) +
            facet_wrap(~ method, scales="free") + theme_bw() +
            theme(strip.background = element_blank(),
                  legend.position = "bottom",
                  legend.title = element_blank(),
                  legend.justification="left",
                  legend.box = "vertical",
                  legend.box.margin = ggplot2::margin(0,0,0,-1,"cm"),
                  legend.spacing.y = ggplot2::unit(0.02, "cm"),
                  legend.key.height = ggplot2::unit(0.25, "cm"),
                  legend.key.width = ggplot2::unit(0.35, "cm"),
                  legend.margin = ggplot2::margin(0, 0, 0, 0, "cm"),
                  plot.subtitle = element_text(size = 6, lineheight = 0.9))

          if(!is.null(mycols)) p <- p + scale_color_manual(values = mycols)
          p <- p + guides(
            color = ggplot2::guide_legend(ncol = 1, byrow = TRUE),
            shape = ggplot2::guide_legend(ncol = 1, byrow = TRUE)
          )
          if (!is.null(ID) && ID %in% names(pdataframe)) p <- p + ggrepel::geom_text_repel(aes_string(label = ID), size = 2)

          if (ellipse == TRUE) p <- p + stat_ellipse(type="norm", linetype=2)
          else if (!is.null(ellipse) && ellipse != TRUE) p <- p + stat_ellipse(aes_string(group=ellipse, color=ellipse), type="norm", linetype=2)

          # ======= PERMANOVA (코너 고정) ======= #
          if (statistics){
            if (!is.null(facet) && facet %in% names(pdataframe)) {
              ann_df <- .perm_ann_by_facet(ps_obj = psIN.cbn.na, pdataframe = pdataframe,
                                           mvar = mvar, facet = facet,
                                           distance_metric = distance_metric,
                                           cate.conf = cate.conf,
                                           project = project, name = name,
                                           strata_var = strata_var)
              p <- p + ggplot2::geom_text(
                data = ann_df,
                mapping = aes(x = x, y = y, label = label),
                size = 3, hjust = -0.005, vjust = -0.3,
                lineheight = 0.95,
                inherit.aes = FALSE
              )
            } else {
              set.seed(1)
              distance <- Go_dist(psIN = psIN.cbn.na, project = project, cate.vars = mvar, name=name, distance_metrics = distance_metric)
              x <- as.dist(distance[[distance_metric]])
              factors <- mapping.sel.na.rem[,mvar]
              x1 <- as.matrix(x)[factors %in% unique(factors), factors %in% unique(factors)]
              map.pair <- subset(mapping.sel.na.rem, mapping.sel.na.rem[,mvar] %in% unique(factors))

              model_vars <- .perm_model_vars(mvar, cate.conf, strata_var)
              model_vars <- intersect(model_vars, names(map.pair))
              if (length(model_vars) > 0) {
                keep <- stats::complete.cases(map.pair[, model_vars, drop = FALSE])
                map.pair <- map.pair[keep, , drop = FALSE]
              }
              if (nrow(map.pair) < 3 || length(unique(map.pair[, mvar])) < 2) {
                next
              }
              x1 <- as.matrix(x)[rownames(map.pair), rownames(map.pair), drop = FALSE]

              perm_use <- 999
              if (!is.null(strata_var) && strata_var %in% names(map.pair)) {
                perm_use <- .mk_perm_block(map.pair[[strata_var]], nperm = 999)
              }

              if (!is.null(cate.conf) && length(cate.conf)>0) {
                for(conf in cate.conf) map.pair[,conf] <- factor(map.pair[,conf])
                form <- as.formula(sprintf("x1 ~ %s + %s", mvar, paste(setdiff(cate.conf,"SampleType"), collapse="+")))
              } else {
                form <- as.formula(sprintf("x1 ~ %s", mvar))
              }
              ad <- vegan::adonis2(form, data = map.pair, permutations = perm_use, by="terms")
              R2 <- round(ad[1,3], 3); p_perm <- ad[1,5]

              # SAVE non-facet PERMANOVA table   # <<< NEW
              fn <- file.path(
                out_perm,
                sprintf("PERMANOVA.%s.%s.%s%s.csv",
                        distance_metric, project, mvar,
                        ifelse(is.null(strata_var), "", paste0(".strata=", strata_var)))
              )
              utils::write.csv(as.data.frame(ad), fn, row.names = TRUE)

              ann_df <- data.frame(
                x = -Inf,
                y = -Inf,
                label = sprintf("%-12s\n%-12s\n%-18s",
                                distance_metric,
                                paste0("R2=", formatC(R2, format="f", digits=3)),
                                paste0("PERMANOVA p=", formatC(p_perm, format="f", digits=3)))
              )
              p <- p + ggplot2::geom_text(
                data = ann_df,
                aes(x = x, y = y, label = label),
                size = 3, hjust = -0.005, vjust = -0.3,
                lineheight = 0.95,
                inherit.aes = FALSE
              )
            }
          }
          # ===================================== #

          if (!is.null(facet) && facet %in% names(pdataframe)) {
            p <- p + facet_wrap(as.formula(sprintf("~ %s", facet)), scales="free_x", ncol = 1)
          }

          if(!is.null(paired) && paired %in% names(pdataframe)){
              p <- p + geom_line(aes_string(group = paired), color="grey", linewidth=0.2,
                                 arrow = arrow(type="closed", length=unit(0.025,"inches")))
          }

          p <- p +
            scale_x_continuous(expand = expansion(mult = c(0.03, 0.03))) +
            scale_y_continuous(expand = expansion(mult = c(0.06, 0.03))) +
            theme(panel.grid = element_blank(),
                  legend.key = element_blank(),
                  panel.background = element_rect(fill = "white", colour = "Black", linewidth = 0.7, linetype = "solid"),
                  aspect.ratio = 1,
                  plot.title=element_text(size=8,face="bold")) +
            geom_vline(xintercept = 0, linewidth = 0.1) + geom_hline(yintercept = 0, linewidth = 0.1)

          plotlist[[length(plotlist)+1]] <- p
        }
      }
      multiplot(plotlist = plotlist, cols = plotCols, rows = plotRows)
      dev.off()

    } else {
      # no combination
      ord_meths = plot
      pdf(sprintf("%s/ordi.%s.%s.%s.%s%s%s%s%s%s%s%s%s.pdf", out_path,
                  ord_meths, paste(distance_metrics, collapse = "+"), project, mvar,
                  ifelse(is.null(facet), "", paste(facet, ".", sep = "")),
                  ifelse(is.null(combination), "", paste("(cbn=",combination, ").", sep = "")),
                  ifelse(is.null(cate.conf), "", paste("with_confounder", ".", sep = "")),
                  ifelse(is.null(paired), "", paste("(paired=",paired, ").", sep = "")),
                  ifelse(is.null(name), "", paste(name, ".", sep = "")),
                  ifelse(ellipse == FALSE, "ellipse_FALSE.",
                         ifelse(ellipse == TRUE, "", paste("ellipse_", ellipse, ".", sep = ""))),
                  ifelse(is.null(strata_var), "", paste("(strata=", strata_var, ").", sep = "")),
                  format(Sys.Date(), "%y%m%d")), height = height, width = width)
      plotlist <- list()
      for(distance_metric in distance_metrics){
        mapping.sel <- data.frame(sample_data(psIN))
        mapping.sel.na <- mapping.sel[!is.na(mapping.sel[,mvar]), ]
        psIN.na <- prune_samples(rownames(mapping.sel[!is.na(mapping.sel[,mvar]), ]), psIN)
        mapping.sel.na.rem <- data.frame(sample_data(psIN.na ))
        mapping.sel.na.rem[,mvar] <- factor(mapping.sel.na.rem[,mvar])

        perm_global <- 999
        if (!is.null(strata_var) && strata_var %in% names(mapping.sel.na.rem)) {
          perm_global <- .mk_perm_block(mapping.sel.na.rem[[strata_var]], nperm = 999)
        }

        ord_meths = plot
        plist = plyr::llply(as.list(ord_meths), function(i, psIN.na, distance_metric){
          ordi = ordinate(psIN.na, method=i, distance=distance_metric)
          df = as.data.frame(ordi$vectors[, 1:2]); colnames(df) = c("Axis_1", "Axis_2")
          var_explained <- .axis_percent(ordi)
          df$Axis1_Percent = var_explained[1]; df$Axis2_Percent = var_explained[2]
          metadata = as.data.frame(sample_data(psIN.na))
          cbind(df, metadata)
        }, psIN.na, distance_metric)

        names(plist) <- ord_meths
        pdataframe = plyr::ldply(plist, identity); names(pdataframe)[1] = "method"

        if (!is.null(facet) && facet %in% names(pdataframe)) {
          pdataframe[, facet] <- factor(pdataframe[, facet],
                                        levels = intersect(orders, unique(pdataframe[, facet])))
        }

        pdataframe[, mvar] <- factor(pdataframe[, mvar],
                                     levels = intersect(orders, unique(pdataframe[, mvar])))

        if(addnumber==TRUE){
          for (Name in unique(pdataframe[,mvar])) {
            total <- sum(pdataframe[,mvar] == Name, na.rm=TRUE)
            levels(pdataframe[[mvar]])[levels(pdataframe[[mvar]])== Name] <- paste0(Name, " (n=", total, ")")
          }
        }

        axis1_percent_avg <- mean(pdataframe$Axis1_Percent, na.rm = TRUE)
        axis2_percent_avg <- mean(pdataframe$Axis2_Percent, na.rm = TRUE)

        if (!is.null(shapes) && shapes %in% names(pdataframe)) {
          pdataframe[,shapes] <- .safe_levels(pdataframe[,shapes], orders)
          p = ggplot(pdataframe, aes_string(x="Axis_1", y="Axis_2", color=mvar)) +
            geom_point(aes_string(shape=shapes), size=0.9, alpha=1) +
            scale_shape_manual(values = c(1,16,8,0,15,2,17,11,10,12,3,4,5,6,7,8,9,13,14))
        } else {
          p = ggplot(pdataframe, aes_string(x="Axis_1", y="Axis_2", color=mvar)) +
            geom_point(size=0.9, alpha=1)
        }

        subtitle_text <- build_plot_subtitle(cate.conf, strata_var)
        p = p +
          labs(x = paste0("Axis 1 (", sprintf("%.2f", axis1_percent_avg),"%)"),
               y = paste0("Axis 2 (", sprintf("%.2f", axis2_percent_avg),"%)"),
               title = sprintf("%s (%s)%s", mvar, distance_metric, title_suffix),
               subtitle = subtitle_text) +
          facet_wrap(~ method, scales="free") + theme_bw() +
          theme(strip.background = element_blank(),
                legend.position = "bottom",
                legend.title = element_blank(),
                legend.justification="left",
                legend.box = "vertical",
                legend.box.margin = ggplot2::margin(0,0,0,-1,"cm"),
                legend.spacing.y = ggplot2::unit(0.02, "cm"),
                legend.key.height = ggplot2::unit(0.25, "cm"),
                legend.key.width = ggplot2::unit(0.35, "cm"),
                legend.margin = ggplot2::margin(0, 0, 0, 0, "cm"),
                plot.title=element_text(size=8,face="bold"),
                plot.subtitle = element_text(size = 6, lineheight = 0.9))

        if(!is.null(mycols)) p <- p + scale_color_manual(values = mycols)
        p <- p + guides(
          color = ggplot2::guide_legend(ncol = 1, byrow = TRUE),
          shape = ggplot2::guide_legend(ncol = 1, byrow = TRUE)
        )
        if (!is.null(ID) && ID %in% names(pdataframe)) p <- p + ggrepel::geom_text_repel(aes_string(label = ID), size = 2)

        if (ellipse == TRUE) p <- p + stat_ellipse(type="norm", linetype=2)
        else if (!is.null(ellipse) && ellipse != TRUE) p <- p + stat_ellipse(aes_string(group=ellipse, color=ellipse), type="norm", linetype=2)

        # ======= PERMANOVA (코너 고정) ======= #
        if (statistics){
          if (!is.null(facet) && facet %in% names(pdataframe)) {
            ann_df <- .perm_ann_by_facet(ps_obj = psIN.na, pdataframe = pdataframe,
                                         mvar = mvar, facet = facet,
                                         distance_metric = distance_metric,
                                         cate.conf = cate.conf,
                                         project = project, name = name,
                                         strata_var = strata_var)
            p <- p + ggplot2::geom_text(
              data = ann_df,
              mapping = aes(x = x, y = y, label = label),
              size = 3, hjust = -0.005, vjust = -0.3,
              lineheight = 0.95,
              inherit.aes = FALSE
            )
          } else {
            set.seed(1)
            distance <- Go_dist(psIN = psIN.na, project = project, name=name, cate.vars = mvar, distance_metrics = distance_metric)
            x <- as.dist(distance[[distance_metric]])
            factors <- mapping.sel.na.rem[,mvar]
            x1 <- as.matrix(x)[factors %in% unique(factors), factors %in% unique(factors)]
            map.pair <- subset(mapping.sel.na.rem, mapping.sel.na.rem[,mvar] %in% unique(factors))

            model_vars <- .perm_model_vars(mvar, cate.conf, strata_var)
            model_vars <- intersect(model_vars, names(map.pair))
            if (length(model_vars) > 0) {
              keep <- stats::complete.cases(map.pair[, model_vars, drop = FALSE])
              map.pair <- map.pair[keep, , drop = FALSE]
            }
            if (nrow(map.pair) < 3 || length(unique(map.pair[, mvar])) < 2) {
              next
            }
            x1 <- as.matrix(x)[rownames(map.pair), rownames(map.pair), drop = FALSE]

            perm_use <- 999
            if (!is.null(strata_var) && strata_var %in% names(map.pair)) {
              perm_use <- .mk_perm_block(map.pair[[strata_var]], nperm = 999)
            }

            if (!is.null(cate.conf) && length(cate.conf)>0) {
              for(conf in cate.conf) map.pair[,conf] <- factor(map.pair[,conf])
              form <- as.formula(sprintf("x1 ~ %s + %s", mvar, paste(setdiff(cate.conf,"SampleType"), collapse="+")))
            } else {
              form <- as.formula(sprintf("x1 ~ %s", mvar))
            }
            ad <- vegan::adonis2(form, data = map.pair, permutations = perm_use, by="terms")
            R2 <- round(ad[1,3], 3); p_perm <- ad[1,5]

            # SAVE non-facet PERMANOVA table   # <<< NEW
            fn <- file.path(
              out_perm,
              sprintf("PERMANOVA.%s.%s.%s%s.csv",
                      distance_metric, project, mvar,
                      ifelse(is.null(strata_var), "", paste0(".strata=", strata_var)))
            )
            utils::write.csv(as.data.frame(ad), fn, row.names = TRUE)

            ann_df <- data.frame(
              x = -Inf,
              y = -Inf,
              label = sprintf(
                "%s\nR2=%.3f\nPERMANOVA p=%.3f",
                distance_metric,
                R2,
                p_perm
              ),
              stringsAsFactors = FALSE
            )
            p <- p + ggplot2::geom_text(
              data = ann_df,
              aes(x = x, y = y, label = label),
              size = 3, hjust = -0.005, vjust = -0.3,
              lineheight = 0.95,
              inherit.aes = FALSE
            )
          }
        }
        # ===================================== #

        if (!is.null(facet) && facet %in% names(pdataframe)) {
          mapping.sel.na.rem[[facet]] <- factor(mapping.sel.na.rem[[facet]], levels = orders)
          p <- p + facet_wrap(as.formula(sprintf("~ %s", facet)), scales="free_x", ncol = 1)
        }

        if(!is.null(paired) && paired %in% names(pdataframe)){
          p <- p + geom_line(aes_string(group = paired),color="grey", linewidth = 0.2,
                             arrow = arrow(type="closed", length=unit(0.025,"inches")))
        }

        p <- p +
          scale_x_continuous(expand = expansion(mult = c(0.03, 0.03))) +
          scale_y_continuous(expand = expansion(mult = c(0.06, 0.03))) +
          theme(panel.grid = element_blank(),
                legend.key = element_blank(),
                panel.background = element_rect(fill = "white", colour = "Black", linewidth = 0.7, linetype = "solid"),
                aspect.ratio = 1) +
          geom_vline(xintercept = 0, linewidth = 0.1) + geom_hline(yintercept = 0, linewidth = 0.1)

        plotlist[[length(plotlist)+1]] <- p
      }
      multiplot(plotlist = plotlist, cols = plotCols, rows = plotRows)
      dev.off()
    }
  }
}
