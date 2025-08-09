
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
#'
#' @details
#' This function creates biodiversity plots such as PCoA for microbiome data, facilitating comparisons across different conditions or factors. The function supports various distance metrics and provides options for customization.
#'
#' @return
#' A PDF file containing the biodiversity plot(s).
#'
#' @examples
#' Go_bdiv(psIN = ps_object,
#'         cate.vars = c("Condition", "Treatment"),
#'         project = "MyMicrobiomeStudy",
#'         orders = c("Condition1", "Condition2"),
#'         distance_metrics = c("bray", "unifrac"),
#'         cate.conf = "AgeGroup",
#'         plot = "PCoA",
#'         ellipse = TRUE, or group names
#'         statistics = TRUE,
#'         mycols = c("blue", "red"),
#'         paired = "PatientID",
#'         combination = 2,
#'         shapes = 16,
#'         ID = "SampleID",
#'         facet = "Group",
#'         name = "BiodivPlot",
#'         addnumber = TRUE,
#'         height = 6,
#'         width = 8)
#'
#' @export

Go_bdiv <- function(psIN, cate.vars, project, orders, distance_metrics,
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
                    height, width){

  if(!is.null(dev.list())) dev.off()

  # out dir
  out <- file.path(sprintf("%s_%s",project, format(Sys.Date(), "%y%m%d")))
  if(!file_test("-d", out)) dir.create(out)
  out_path <- file.path(sprintf("%s_%s/pdf",project, format(Sys.Date(), "%y%m%d")))
  if(!file_test("-d", out_path)) dir.create(out_path)
  out_dist <- file.path(sprintf("%s_%s/table/dist",project, format(Sys.Date(), "%y%m%d")))
  if(!file_test("-d", out_dist)) dir.create(out_dist)

  # ---------- helpers ----------
  .safe_levels <- function(x, levs) {
    x <- factor(x)
    if (!is.null(levs)) x <- factor(x, levels = intersect(levs, levels(x)))
    x
  }

  # facet 있을 때: facet별 PERMANOVA + 코너 고정 좌표(-Inf, -Inf)
  .perm_ann_by_facet <- function(ps_obj, pdataframe, mvar, facet, distance_metric, cate.conf, project, name){
    levs <- unique(pdataframe[[facet]])
    stat_list <- lapply(levs, function(flv){
      map_sub <- as.data.frame(sample_data(ps_obj))
      map_sub <- map_sub[map_sub[[facet]] == flv, , drop = FALSE]
      if (nrow(map_sub) < 3 || length(unique(map_sub[[mvar]])) < 2) {
        return(data.frame(facet_val = flv, R2 = NA_real_, p = NA_real_))
      }
      ps_sub <- prune_samples(rownames(map_sub), ps_obj)

      dist_list <- Go_dist(psIN = ps_sub, project = project, name = name,
                           cate.vars = mvar, distance_metrics = distance_metric)
      x <- as.dist(dist_list[[distance_metric]])

      map_sub2 <- data.frame(sample_data(ps_sub))
      if (!is.null(cate.conf) && length(cate.conf) > 0) {
        for (conf in cate.conf) map_sub2[[conf]] <- factor(map_sub2[[conf]])
        form <- as.formula(sprintf("x ~ %s + %s", mvar,
                                   paste(setdiff(cate.conf, "SampleType"), collapse = " + ")))
      } else {
        form <- as.formula(sprintf("x ~ %s", mvar))
      }

      ad <- vegan::adonis2(form, data = map_sub2, permutations = 999, by = "terms")
      data.frame(facet_val = flv,
                 R2 = suppressWarnings(round(ad[1,"R2"], 3)),
                 p  = suppressWarnings(ad[1,"Pr(>F)"]))
    })

    stat_df <- do.call(rbind, stat_list)
    stat_df$padj <- p.adjust(stat_df$p, method = "bonferroni")

    ann_df <- stat_df |>
      dplyr::mutate(
        label = paste0(distance_metric,
                       "\nR2=", ifelse(is.na(R2), "NA", format(R2, nsmall = 3)),
                       "\nPERMANOVA p=", ifelse(is.na(padj), "NA", signif(padj,3))),
        x = -Inf,  # 코너 고정
        y = -Inf
      )

    ann_df[[facet]] <- ann_df$facet_val
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

    if (length(facet) >= 1 && facet == mvar) next
    if (length(shapes) >= 1 && shapes == mvar) next

    #------------------------------#
    # for group combination or not #
    #------------------------------#
    if (!is.null(combination)){
      mapping[,mvar] <- .safe_levels(mapping[,mvar], orders)
      group.cbn <- combn(x = levels(mapping[,mvar]), m = combination)
      group_comparisons <- lapply(seq_len(ncol(group.cbn)), function(i) group.cbn[,i])

      ord_meths = plot
      pdf(sprintf("%s/ordi.%s.%s.%s.%s%s%s%s%s%s%s%s.pdf", out_path,
                  ord_meths, "distance_metric", project, mvar,
                  ifelse(is.null(facet), "", paste(facet, ".", sep = "")),
                  ifelse(is.null(combination), "", paste("(cbn=",combination, ").", sep = "")),
                  ifelse(is.null(cate.conf), "", paste("with_confounder", ".", sep = "")),
                  ifelse(is.null(paired), "", paste("(paired=",paired, ").", sep = "")),
                  ifelse(is.null(name), "", paste(name, ".", sep = "")),
                  ifelse(ellipse == FALSE, "ellipse_FALSE.",
                         ifelse(ellipse == TRUE, "", paste("ellipse_", ellipse, ".", sep = ""))),
                  format(Sys.Date(), "%y%m%d")), height = height, width = width)

      for(i in seq_along(group_comparisons)){
        group.combination <- unlist(group_comparisons[i])

        if(combination == 2){
          basline <- group.combination[1]; smvar <- group.combination[2]
          mapping.cbn <- subset(mapping, mapping[,mvar] %in% c(basline, smvar))
        } else if(combination == 3){
          basline <- group.combination[1]; smvar1 <- group.combination[2]; smvar2 <- group.combination[3]
          mapping.cbn <- subset(mapping, mapping[,mvar] %in% c(basline, smvar1, smvar2))
        } else if(combination == 4){
          basline <- group.combination[1]; smvar1 <- group.combination[2]; smvar2 <- group.combination[3]; smvar3 <- group.combination[4]
          mapping.cbn <- subset(mapping, mapping[,mvar] %in% c(basline, smvar1, smvar2, smvar3))
        } else { break }

        psIN.cbn <- psIN; sample_data(psIN.cbn) <- mapping.cbn

        for(distance_metric in distance_metrics){
          # remove NA
          mapping.sel <- data.frame(sample_data(psIN.cbn))
          mapping.sel.na <- mapping.sel[!is.na(mapping.sel[,mvar]), ]
          psIN.cbn.na <- prune_samples(rownames(mapping.sel[!is.na(mapping.sel[,mvar]), ]), psIN.cbn)
          mapping.sel.na.rem <- data.frame(sample_data(psIN.cbn.na ))
          mapping.sel.na.rem[,mvar] <- factor(mapping.sel.na.rem[,mvar])

          # ordination
          ord_meths = plot
          plist = plyr::llply(as.list(ord_meths), function(i, psIN.na, distance_metric){
            ordi = ordinate(psIN.na, method=i, distance=distance_metric)
            df = as.data.frame(ordi$vectors[, 1:2]); colnames(df) = c("Axis_1", "Axis_2")
            if ("Eigenvalues" %in% names(ordi$values)) {
              var_explained = ordi$values$Eigenvalues / sum(ordi$values$Eigenvalues) * 100
              df$Axis1_Percent = var_explained[1]; df$Axis2_Percent = var_explained[2]
            }
            metadata = as.data.frame(sample_data(psIN.na))
            cbind(df, metadata)
          }, psIN.cbn.na, distance_metric)

          names(plist) <- ord_meths
          pdataframe = plyr::ldply(plist, identity); names(pdataframe)[1] = "method"

          if (!is.null(facet) && facet %in% names(pdataframe)) {
            pdataframe[,facet] <- .safe_levels(pdataframe[,facet], orders)
          }
          pdataframe[,mvar] <- .safe_levels(pdataframe[,mvar], orders)

          # n 표시(전체 기준; facet별 n표시가 필요하면 말해줘서 바꿀 수 있음)
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

          p = p +
            labs(x = paste0("Axis 1 (", sprintf("%.2f", axis1_percent_avg),"%)"),
                 y = paste0("Axis 2 (", sprintf("%.2f", axis2_percent_avg),"%)")) +
            ggtitle(sprintf("%s (%s)", mvar, distance_metric)) +
            facet_wrap(~ method, scales="free") + theme_bw() +
            theme(strip.background = element_blank(),
                  legend.position = "bottom",
                  legend.title = element_blank(),
                  legend.justification="left",
                  legend.box = "vertical",
                  legend.box.margin = ggplot2::margin(0,0,0,-1,"cm"))

          if(!is.null(mycols)) p <- p + scale_color_manual(values = mycols)
          if (!is.null(ID) && ID %in% names(pdataframe)) p <- p + ggrepel::geom_text_repel(aes_string(label = ID), size = 2)

          # ellipse
          if (ellipse == TRUE) p <- p + stat_ellipse(type="norm", linetype=2)
          else if (!is.null(ellipse) && ellipse != TRUE) p <- p + stat_ellipse(aes_string(group=ellipse, color=ellipse), type="norm", linetype=2)

          # ======= PERMANOVA (코너 고정) ======= #
          if (statistics){
            if (!is.null(facet) && facet %in% names(pdataframe)) {
              ann_df <- .perm_ann_by_facet(ps_obj = psIN.cbn.na, pdataframe = pdataframe,
                                           mvar = mvar, facet = facet,
                                           distance_metric = distance_metric,
                                           cate.conf = cate.conf,
                                           project = project, name = name)

              p <- p + ggplot2::geom_text(
                data = ann_df,
                mapping = aes(x = x, y = y, label = label),
                size = 3, hjust = -0.1, vjust = -0.6,
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
              if (!is.null(cate.conf) && length(cate.conf)>0) {
                for(conf in cate.conf) map.pair[,conf] <- factor(map.pair[,conf])
                form <- as.formula(sprintf("x1 ~ %s + %s", mvar, paste(setdiff(cate.conf,"SampleType"), collapse="+")))
              } else {
                form <- as.formula(sprintf("x1 ~ %s", mvar))
              }
              ad <- vegan::adonis2(form, data = map.pair, permutations=999, by="terms")
              R2 <- round(ad[1,3], 3); padj <- ad[1,5]
              ann_df <- data.frame(
                x = -Inf,
                y = -Inf,
                label = paste0(distance_metric, "\nR2=",R2,"\nPERMANOVA p=", signif(padj,3))
              )
              p <- p + ggplot2::geom_text(
                data = ann_df,
                aes(x = x, y = y, label = label),
                size = 3, hjust = -0.1, vjust = -0.6,
                lineheight = 0.95,
                inherit.aes = FALSE
              )
            }
          }
          # ===================================== #

          if (!is.null(facet) && facet %in% names(pdataframe)) {
            ncol <- length(unique(mapping.sel.na.rem[,facet]))
            p <- p + facet_wrap(as.formula(sprintf("~ %s", facet)), scales="free_x", ncol = ncol)
          }

          if(!is.null(paired) && paired %in% names(pdataframe)){
            p <- p + geom_line(aes_string(group = paired), color="grey", size=0.2,
                               arrow = arrow(type="closed", length=unit(0.025,"inches")))
          }

          # 코너 고정 텍스트가 잘 보이도록 여백
          p <- p +
            scale_x_continuous(expand = expansion(mult = c(0.03, 0.03))) +
            scale_y_continuous(expand = expansion(mult = c(0.06, 0.03))) +
            theme(panel.grid = element_blank(),
                  legend.key = element_blank(),
                  panel.background = element_rect(fill = "white", colour = "Black", size = 0.7, linetype = "solid"),
                  aspect.ratio = 1,
                  plot.title=element_text(size=8,face="bold")) +
            geom_vline(xintercept = 0, size = 0.1) + geom_hline(yintercept = 0, size = 0.1)

          print(p)
        }
      }
      dev.off()

    } else {
      # no combination
      for(distance_metric in distance_metrics){
        mapping.sel <- data.frame(sample_data(psIN))
        mapping.sel.na <- mapping.sel[!is.na(mapping.sel[,mvar]), ]
        psIN.na <- prune_samples(rownames(mapping.sel[!is.na(mapping.sel[,mvar]), ]), psIN)
        mapping.sel.na.rem <- data.frame(sample_data(psIN.na ))
        mapping.sel.na.rem[,mvar] <- factor(mapping.sel.na.rem[,mvar])

        ord_meths = plot
        plist = plyr::llply(as.list(ord_meths), function(i, psIN.na, distance_metric){
          ordi = ordinate(psIN.na, method=i, distance=distance_metric)
          df = as.data.frame(ordi$vectors[, 1:2]); colnames(df) = c("Axis_1", "Axis_2")
          if ("Eigenvalues" %in% names(ordi$values)) {
            var_explained = ordi$values$Eigenvalues / sum(ordi$values$Eigenvalues) * 100
            df$Axis1_Percent = var_explained[1]; df$Axis2_Percent = var_explained[2]
          }
          metadata = as.data.frame(sample_data(psIN.na))
          cbind(df, metadata)
        }, psIN.na, distance_metric)

        names(plist) <- ord_meths
        pdataframe = plyr::ldply(plist, identity); names(pdataframe)[1] = "method"

        if (!is.null(facet) && facet %in% names(pdataframe)) {
          pdataframe[,facet] <- .safe_levels(pdataframe[,facet], orders)
        }
        pdataframe[,mvar] <- .safe_levels(pdataframe[,mvar], orders)

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

        p = p +
          labs(x = paste0("Axis 1 (", sprintf("%.2f", axis1_percent_avg),"%)"),
               y = paste0("Axis 2 (", sprintf("%.2f", axis2_percent_avg),"%)")) +
          ggtitle(sprintf("%s (%s)", mvar, distance_metric)) +
          facet_wrap(~ method, scales="free") + theme_bw() +
          theme(strip.background = element_blank(),
                legend.position = "bottom",
                legend.title = element_blank(),
                legend.justification="left",
                legend.box = "vertical",
                legend.box.margin = ggplot2::margin(0,0,0,-1,"cm"),
                plot.title=element_text(size=8,face="bold"))

        if(!is.null(mycols)) p <- p + scale_color_manual(values = mycols)
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
                                         project = project, name = name)
            p <- p + ggplot2::geom_text(
              data = ann_df,
              mapping = aes(x = x, y = y, label = label),
              size = 3, hjust = -0.1, vjust = -0.6,
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
            if (!is.null(cate.conf) && length(cate.conf)>0) {
              for(conf in cate.conf) map.pair[,conf] <- factor(map.pair[,conf])
              form <- as.formula(sprintf("x1 ~ %s + %s", mvar, paste(setdiff(cate.conf,"SampleType"), collapse="+")))
            } else {
              form <- as.formula(sprintf("x1 ~ %s", mvar))
            }
            ad <- vegan::adonis2(form, data = map.pair, permutations=999, by="terms")
            R2 <- round(ad[1,3], 3); padj <- ad[1,5]
            ann_df <- data.frame(
              x = -Inf,
              y = -Inf,
              label = paste0(distance_metric, "\nR2=",R2,"\nPERMANOVA p=", signif(padj,3))
            )
            p <- p + ggplot2::geom_text(
              data = ann_df,
              aes(x = x, y = y, label = label),
              size = 3, hjust = -0.1, vjust = -0.6,
              lineheight = 0.95,
              inherit.aes = FALSE
            )
          }
        }
        # ===================================== #

        if (!is.null(facet) && facet %in% names(pdataframe)) {
          ncol <- length(unique(mapping.sel.na.rem[,facet]))
          p <- p + facet_wrap(as.formula(sprintf("~ %s", facet)), scales="free_x", ncol = ncol)
        }

        if(!is.null(paired) && paired %in% names(pdataframe)){
          p <- p + geom_line(aes_string(group = paired),color="grey", size = 0.2,
                             arrow = arrow(type="closed", length=unit(0.025,"inches")))
        }

        # 코너 고정 텍스트가 잘 보이도록 여백
        p <- p +
          scale_x_continuous(expand = expansion(mult = c(0.03, 0.03))) +
          scale_y_continuous(expand = expansion(mult = c(0.06, 0.03))) +
          theme(panel.grid = element_blank(),
                legend.key = element_blank(),
                panel.background = element_rect(fill = "white", colour = "Black", size = 0.7, linetype = "solid"),
                aspect.ratio = 1) +
          geom_vline(xintercept = 0, size = 0.1) + geom_hline(yintercept = 0, size = 0.1)

        pdf(sprintf("%s/ordi.%s.%s.%s.%s%s%s%s%s%s%s%s.pdf", out_path,
                    ord_meths, distance_metric, project, mvar,
                    ifelse(is.null(facet), "", paste(facet, ".", sep = "")),
                    ifelse(is.null(combination), "", paste("(cbn=",combination, ").", sep = "")),
                    ifelse(is.null(cate.conf), "", paste("with_confounder", ".", sep = "")),
                    ifelse(is.null(paired), "", paste("(paired=",paired, ").", sep = "")),
                    ifelse(is.null(name), "", paste(name, ".", sep = "")),
                    ifelse(ellipse == FALSE, "ellipse_FALSE.",
                           ifelse(ellipse == TRUE, "", paste("ellipse_", ellipse, ".", sep = ""))),
                    format(Sys.Date(), "%y%m%d")), height = height, width = width)
        print(p)
        dev.off()
      }
    }
  }
}
