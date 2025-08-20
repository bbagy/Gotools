#' Beta-diversity with MiRKAT / MiRKAT_LMM (facet annotation + CSV saving)
#'
#' @param psIN      phyloseq object
#' @param cate.vars grouping variables to test (vector)
#' @param project   project name (used in filenames)
#' @param orders    factor order vector (optional)
#' @param distance_metrics e.g. c("bray","unifrac","wunifrac",...)
#' @param cate.conf covariates (character vector) [currently not used inside MiRKAT core test]
#' @param plot      ordination method(s), e.g. "PCoA"
#' @param ellipse   TRUE/FALSE or grouping name for ellipse color
#' @param statistics TRUE to compute MiRKAT/MiRKAT_LMM p-values
#' @param mycols    manual color vector
#' @param rand.eff  random-effect column for MiRKAT_LMM (e.g., "UID"); NULL → MiRKAT
#' @param combination NULL or 2/3/4 (subset pages; single PDF)
#' @param shapes    shape mapping column (optional)
#' @param ID        text label column (optional)
#' @param facet     facet column (optional)
#' @param name      extra name suffix (optional)
#' @param addnumber add n to group labels
#' @param height,width PDF size
#'
#' @return PDFs and CSVs written to disk
#' @export
Go_bdivMK <- function(psIN, cate.vars, project, orders, distance_metrics,
                      cate.conf=NULL,
                      plot="PCoA",
                      ellipse=TRUE,
                      statistics = TRUE,
                      mycols=NULL,
                      rand.eff = NULL,     # <- switch: NULL=MiRKAT, not NULL=MiRKAT_LMM (if available)
                      combination=NULL,
                      shapes = NULL,
                      ID = NULL,
                      facet=NULL,
                      name=NULL,
                      addnumber=TRUE,
                      height, width){

  if(!is.null(dev.list())) dev.off()

  # out dirs
  out_root  <- sprintf("%s_%s", project, format(Sys.Date(), "%y%m%d"))
  out_path  <- file.path(out_root, "pdf");        if(!file_test("-d", out_path))  dir.create(out_path, recursive = TRUE)
  out_dist  <- file.path(out_root, "table", "dist");   if(!file_test("-d", out_dist)) dir.create(out_dist, recursive = TRUE)
  out_table <- file.path(out_root, "table", "mirkat"); if(!file_test("-d", out_table)) dir.create(out_table, recursive = TRUE)

  title_suffix <- if (!is.null(rand.eff)) sprintf(" | rand=%s", rand.eff) else ""

  # ---------- helpers ----------
  .safe_levels <- function(x, levs) {
    x <- factor(x)
    if (is.null(levs)) return(x)
    keep <- levs[levs %in% levels(x)]
    factor(x, levels = keep)
  }

  .mk_kernel <- function(dist_obj){
    # dist_obj: class 'dist' or dist matrix
    D <- if (inherits(dist_obj,"dist")) as.matrix(dist_obj) else as.matrix(dist_obj)
    # MiRKAT::D2K(D)  (no extra args; earlier error was from passing an extra "Gaussian")
    K <- try(MiRKAT::D2K(D), silent = TRUE)
    if (inherits(K, "try-error")) return(NULL) else return(K)
  }

  .mk_y_from_factor <- function(vec){
    f <- factor(vec)
    k <- nlevels(f)
    if (k < 2) return(NULL)
    if (k == 2) {
      # binary: 0/1
      as.numeric(f) - 1L
    } else {
      # 3+ levels: simple numeric coding (1,2,3,...) — continuous approx.
      # If you prefer multi-df joint test, say the word; we can expand.
      as.numeric(f)
    }
  }

  .mk_label_block <- function(metric, pval, method_tag){
    ptxt <- ifelse(is.na(pval), "NA", formatC(pval, format="f", digits=3))
    sprintf("%s\n%s p=%s", metric, method_tag, ptxt)
  }

  # facet 주석 데이터프레임: 각 facet 레벨별 한 행
  .mk_ann_facet <- function(pdata, facet, label_text){
    levs <- unique(pdata[[facet]])
    ann <- data.frame(
      x = -Inf,
      y = -Inf,
      label = rep(label_text, length(levs)),
      stringsAsFactors = FALSE
    )
    ann[[facet]] <- levs
    ann
  }

  # facet별 MiRKAT / MiRKAT_LMM
  .mirkat_by_facet <- function(ps_obj, mvar, facet, metric, rand.eff){
    map_all <- as.data.frame(phyloseq::sample_data(ps_obj))
    # facet를 문자 → factor로 안전 변환
    if (!is.factor(map_all[[facet]])) map_all[[facet]] <- factor(map_all[[facet]])
    levs <- levels(droplevels(map_all[[facet]]))
    out_list <- lapply(levs, function(flv){
      map_sub <- as.data.frame(phyloseq::sample_data(ps_obj))
      map_sub <- map_sub[map_sub[[facet]] == flv, , drop = FALSE]
      if (nrow(map_sub) < 3 || length(unique(map_sub[[mvar]])) < 2)
        return(data.frame(facet_val = flv, p = NA_real_))

      ps_sub <- phyloseq::prune_samples(rownames(map_sub), ps_obj)
      # distance
      dist_list <- Go_dist(psIN = ps_sub, project = project, name = name,
                           cate.vars = mvar, distance_metrics = metric)
      D <- dist_list[[metric]]
      K <- .mk_kernel(D)
      if (is.null(K))
        return(data.frame(facet_val = flv, p = NA_real_))

      map2 <- as.data.frame(phyloseq::sample_data(ps_sub))
      y    <- .mk_y_from_factor(map2[[mvar]])
      if (is.null(y))
        return(data.frame(facet_val = flv, p = NA_real_))

      pval <- NA_real_
      if (is.null(rand.eff)) {
        # MiRKAT
        p_try <- try(MiRKAT::MiRKAT(y = y, Ks = list(K), out_type = if (nlevels(factor(map2[[mvar]]))==2) "D" else "C"),
                     silent = TRUE)
        if (!inherits(p_try,"try-error") && is.list(p_try) && length(p_try) > 0) {
          # First kernel p-value
          pv <- suppressWarnings(as.numeric(p_try[[1]]))
          if (is.finite(pv)) pval <- pv
        }
      } else {
        # MiRKAT_LMM (only if available)
        if (isTRUE(requireNamespace("MiRKAT", quietly = TRUE)) &&
            "MiRKAT_LMM" %in% getNamespaceExports("MiRKAT")) {
          Z <- as.factor(map2[[rand.eff]])
          p_try <- try(MiRKAT::MiRKAT_LMM(y = y, K = K, id = Z),
                       silent = TRUE)
          if (!inherits(p_try,"try-error")) {
            pv <- suppressWarnings(as.numeric(p_try))
            if (is.finite(pv)) pval <- pv
          }
        } else {
          message("[MiRKAT_LMM] not available in your MiRKAT version; returning NA p.")
        }
      }

      data.frame(facet_val = flv, p = pval, stringsAsFactors = FALSE)
    })
    dplyr::bind_rows(out_list)
  }

  # 비-facet MiRKAT / MiRKAT_LMM (한 번 계산)
  .mirkat_simple <- function(ps_obj, mvar, metric, rand.eff){
    mapping <- as.data.frame(phyloseq::sample_data(ps_obj))
    if (nrow(mapping) < 3 || length(unique(mapping[[mvar]])) < 2)
      return(NA_real_)
    dist_list <- Go_dist(psIN = ps_obj, project = project, name = name,
                         cate.vars = mvar, distance_metrics = metric)
    D <- dist_list[[metric]]
    K <- .mk_kernel(D)
    if (is.null(K)) return(NA_real_)

    y <- .mk_y_from_factor(mapping[[mvar]])
    if (is.null(y)) return(NA_real_)

    if (is.null(rand.eff)) {
      p_try <- try(MiRKAT::MiRKAT(y = y, Ks = list(K), out_type = if (nlevels(factor(mapping[[mvar]]))==2) "D" else "C"),
                   silent = TRUE)
      if (!inherits(p_try,"try-error") && is.list(p_try) && length(p_try) > 0) {
        pv <- suppressWarnings(as.numeric(p_try[[1]]))
        if (is.finite(pv)) return(pv)
      }
      return(NA_real_)
    } else {
      if (isTRUE(requireNamespace("MiRKAT", quietly = TRUE)) &&
          "MiRKAT_LMM" %in% getNamespaceExports("MiRKAT")) {
        Z <- as.factor(mapping[[rand.eff]])
        p_try <- try(MiRKAT::MiRKAT_LMM(y = y, K = K, id = Z),
                     silent = TRUE)
        if (!inherits(p_try,"try-error")) {
          pv <- suppressWarnings(as.numeric(p_try))
          if (is.finite(pv)) return(pv)
        }
        return(NA_real_)
      } else {
        message("[MiRKAT_LMM] not available; returning NA p.")
        return(NA_real_)
      }
    }
  }
  # ------------------------------------------

  # sanity
  tt <- try(mycols, TRUE); if(inherits(tt, "try-error")) mycols <- NULL
  tt <- try(orders, TRUE); if(inherits(tt, "try-error")) orders <- NULL

  # 한 번의 pdf()로 여러 페이지 찍는 구조(= Go_bdivPM 스타일)
  # combination=NULL인 경우: metric마다 1 PDF
  # combination!=NULL인 경우: metric마다 1 PDF 안에 그룹조합을 페이지로
  for (metric in distance_metrics) {
    ord_meths <- plot
    pdf(file = sprintf("%s/ordi.%s.%s.%s%s%s%s%s%s.pdf", out_path,
                       ord_meths, metric, project,
                       ifelse(is.null(facet), "", paste0(".", facet)),
                       ifelse(is.null(combination), "", paste0(".(cbn=", combination, ")")),
                       ifelse(is.null(name), "", paste0(".", name)),
                       ifelse(is.null(rand.eff), "", paste0(".(RE=", rand.eff, ")")),
                       format(Sys.Date(), "%y%m%d")),
        height = height, width = width)

    for (mvar in cate.vars) {
      mapping <- data.frame(sample_data(psIN))
      mapping[,mvar] <- factor(mapping[,mvar])
      sample_data(psIN) <- mapping

      if (!is.null(facet) && facet == mvar) next
      if (!is.null(shapes) && shapes == mvar) next

      # -------- 조합 모드 --------
      if (!is.null(combination)) {
        mapping[,mvar] <- .safe_levels(mapping[,mvar], orders)
        group.cbn <- combn(x = levels(mapping[,mvar]), m = combination)
        group_comparisons <- lapply(seq_len(ncol(group.cbn)), function(i) group.cbn[,i])

        for (i in seq_along(group_comparisons)) {
          sel <- unlist(group_comparisons[i])
          map_cbn <- subset(mapping, mapping[,mvar] %in% sel)
          ps_cbn  <- psIN; sample_data(ps_cbn) <- map_cbn

          # NA 제거
          map_na  <- data.frame(sample_data(ps_cbn))
          map_na  <- map_na[!is.na(map_na[,mvar]), ]
          ps_na   <- prune_samples(rownames(map_na), ps_cbn)
          map_rem <- data.frame(sample_data(ps_na))
          map_rem[,mvar] <- factor(map_rem[,mvar])

          # ordination 좌표
          plist <- plyr::llply(as.list(ord_meths), function(i, ps_obj, metric){
            ordi = ordinate(ps_obj, method=i, distance=metric)
            df = as.data.frame(ordi$vectors[, 1:2]); colnames(df) = c("Axis_1", "Axis_2")
            if ("Eigenvalues" %in% names(ordi$values)) {
              var_explained = ordi$values$Eigenvalues / sum(ordi$values$Eigenvalues) * 100
              df$Axis1_Percent = var_explained[1]; df$Axis2_Percent = var_explained[2]
            }
            metadata = as.data.frame(sample_data(ps_obj))
            cbind(df, metadata)
          }, ps_na, metric)
          names(plist) <- ord_meths
          pdata <- plyr::ldply(plist, identity); names(pdata)[1] = "method"

          # facet/levels
          if (!is.null(facet) && facet %in% names(pdata)) {
            pdata[, facet] <- factor(pdata[, facet], levels = intersect(orders, unique(pdata[, facet])))
          }
          pdata[, mvar] <- factor(pdata[, mvar], levels = intersect(orders, unique(pdata[, mvar])))

          # n 라벨
          if (isTRUE(addnumber)) {
            for (nm in unique(pdata[,mvar])) {
              total <- sum(pdata[,mvar] == nm, na.rm=TRUE)
              levels(pdata[[mvar]])[levels(pdata[[mvar]])== nm] <- paste0(nm, " (n=", total, ")")
            }
          }




          # 플롯
          axis1_percent_avg <- mean(pdata$Axis1_Percent, na.rm = TRUE)
          axis2_percent_avg <- mean(pdata$Axis2_Percent, na.rm = TRUE)
          p <- ggplot(pdata, aes_string(x="Axis_1", y="Axis_2", color=mvar)) +
            geom_point(size=0.9, alpha=1) +
            labs(x = paste0("Axis 1 (", sprintf("%.2f", axis1_percent_avg),"%)"),
                 y = paste0("Axis 2 (", sprintf("%.2f", axis2_percent_avg),"%)")) +
            ggtitle(sprintf("%s (%s)%s", mvar, metric, title_suffix)) +
            facet_wrap(~ method, scales="free") + theme_bw() +
            theme(strip.background = element_blank(),
                  legend.position = "bottom",
                  legend.title = element_blank(),
                  legend.justification="left",
                  legend.box = "vertical",
                  legend.box.margin = ggplot2::margin(0,0,0,-1,"cm"),
                  plot.title=element_text(size=8,face="bold"))

          if(!is.null(mycols)) p <- p + scale_color_manual(values = mycols)
          if (!is.null(ellipse)) {
            if (identical(ellipse, TRUE)) p <- p + stat_ellipse(type="norm", linetype=2)
            else if (is.character(ellipse) && ellipse %in% names(pdata))
              p <- p + stat_ellipse(aes_string(group=ellipse, color=ellipse), type="norm", linetype=2)
          }
          if (!is.null(ID) && ID %in% names(pdata)) p <- p + ggrepel::geom_text_repel(aes_string(label = ID), size = 2)

          # 통계: MiRKAT / MiRKAT_LMM
          if (isTRUE(statistics)) {
            method_tag <- if (is.null(rand.eff)) "MiRKAT" else "MiRKAT-LMM"

            if (!is.null(facet) && facet %in% names(pdata)) {
              # facet별
              res <- .mirkat_by_facet(ps_obj = ps_na, mvar = mvar, facet = facet, metric = metric, rand.eff = rand.eff)

              # 저장 (facet별)
              for (j in seq_len(nrow(res))) {
                flev <- res$facet_val[j]
                fn <- file.path(out_table,
                                sprintf("%s.%s.%s.%s.facet=%s%s.csv",
                                        method_tag, metric, project, mvar, flev,
                                        ifelse(is.null(rand.eff), "", paste0(".rand=", rand.eff))))
                utils::write.csv(res[j, , drop=FALSE], fn, row.names = FALSE)
              }

              # 주석(모든 facet에 동일 문자열을 찍는게 아니라, **각 facet마다 해당 p**를 찍고 싶으면 facet-join 필요)
              # 간단히: facet 각각에 "metric\nmethod p=..."을 같은 텍스트로 찍는 대신,
              # 여기서는 레벨별로 붙이도록 합칩니다.
              ann <- data.frame(
                x = -Inf, y = -Inf,
                label = sprintf("%s\n%s p=%s",
                                metric, method_tag,
                                ifelse(is.na(res$p), "NA", formatC(res$p, format="f", digits=3))),
                stringsAsFactors = FALSE
              )
              ann[[facet]] <- res$facet_val
              ann[[facet]] <- factor(ann[, facet], levels = intersect(orders, unique(ann[, facet])))
              p <- p + ggplot2::geom_text(
                data = ann,
                mapping = aes(x = x, y = y, label = label),
                size = 3, hjust = -0.005, vjust = -0.3,
                lineheight = 0.95,
                inherit.aes = FALSE
              )
            } else {
              # 비-facet
              pval <- .mirkat_simple(ps_obj = ps_na, mvar = mvar, metric = metric, rand.eff = rand.eff)

              # 저장
              fn <- file.path(out_table,
                              sprintf("%s.%s.%s.%s%s.csv",
                                      method_tag, metric, project, mvar,
                                      ifelse(is.null(rand.eff), "", paste0(".rand=", rand.eff))))
              utils::write.csv(data.frame(metric=metric, mvar=mvar, p=pval, method=method_tag),
                               fn, row.names = FALSE)

              # 주석
              p <- p + ggplot2::geom_text(
                data = data.frame(x=-Inf, y=-Inf,
                                  label=.mk_label_block(metric, pval, method_tag)),
                aes(x=x,y=y,label=label),
                size=3, hjust=-0.005, vjust=-0.3, lineheight=0.95,
                inherit.aes = FALSE
              )
            }
          }

          # facet wrap 유지 (열 갯수는 facet 레벨 수로 자동 — 원본 유지)
          if (!is.null(facet) && facet %in% names(pdata)) {
            ncol <- length(unique(pdata[,facet]))
            pdata[[facet]] <- factor(pdata[[facet]], intersect(orders, unique(pdata[, facet])))
            p <- p + facet_wrap(as.formula(sprintf("~ %s", facet)), scales = "free_x", ncol = ncol)
          }

          # 마진/축
          p <- p +
            scale_x_continuous(expand = expansion(mult = c(0.03, 0.03))) +
            scale_y_continuous(expand = expansion(mult = c(0.06, 0.03))) +
            theme(panel.grid = element_blank(),
                  legend.key = element_blank(),
                  panel.background = element_rect(fill = "white", colour = "Black", linewidth = 0.7, linetype = "solid"),
                  aspect.ratio = 1) +
            geom_vline(xintercept = 0, linewidth = 0.1) +
            geom_hline(yintercept = 0, linewidth = 0.1)

          print(p)  # <- 한 PDF 안에 페이지로 누적
        }

        # -------- 전체 모드 --------
      } else {
        # 전체 데이터
        mapping.sel <- data.frame(sample_data(psIN))
        mapping.sel.na <- mapping.sel[!is.na(mapping.sel[,mvar]), ]
        psIN.na <- prune_samples(rownames(mapping.sel.na), psIN)
        mapping.sel.na.rem <- data.frame(sample_data(psIN.na))
        mapping.sel.na.rem[,mvar] <- factor(mapping.sel.na.rem[,mvar])

        # ordination 좌표
        plist <- plyr::llply(as.list(ord_meths), function(i, ps_obj, metric){
          ordi = ordinate(ps_obj, method=i, distance=metric)
          df = as.data.frame(ordi$vectors[, 1:2]); colnames(df) = c("Axis_1", "Axis_2")
          if ("Eigenvalues" %in% names(ordi$values)) {
            var_explained = ordi$values$Eigenvalues / sum(ordi$values$Eigenvalues) * 100
            df$Axis1_Percent = var_explained[1]; df$Axis2_Percent = var_explained[2]
          }
          metadata = as.data.frame(sample_data(ps_obj))
          cbind(df, metadata)
        }, psIN.na, metric)
        names(plist) <- ord_meths
        pdata <- plyr::ldply(plist, identity); names(pdata)[1] = "method"

        if (!is.null(facet) && facet %in% names(pdata)) {
          pdata[, facet] <- factor(pdata[, facet], levels = intersect(orders, unique(pdata[, facet])))
        }
        pdata[, mvar] <- factor(pdata[, mvar], levels = intersect(orders, unique(pdata[, mvar])))

        if (isTRUE(addnumber)) {
          for (nm in unique(pdata[,mvar])) {
            total <- sum(pdata[,mvar] == nm, na.rm=TRUE)
            levels(pdata[[mvar]])[levels(pdata[[mvar]])== nm] <- paste0(nm, " (n=", total, ")")
          }
        }

        axis1_percent_avg <- mean(pdata$Axis1_Percent, na.rm = TRUE)
        axis2_percent_avg <- mean(pdata$Axis2_Percent, na.rm = TRUE)

        p <- ggplot(pdata, aes_string(x="Axis_1", y="Axis_2", color=mvar)) +
          geom_point(size=0.9, alpha=1) +
          labs(x = paste0("Axis 1 (", sprintf("%.2f", axis1_percent_avg),"%)"),
               y = paste0("Axis 2 (", sprintf("%.2f", axis2_percent_avg),"%)")) +
          ggtitle(sprintf("%s (%s)%s", mvar, metric, title_suffix)) +
          facet_wrap(~ method, scales="free") + theme_bw() +
          theme(strip.background = element_blank(),
                legend.position = "bottom",
                legend.title = element_blank(),
                legend.justification="left",
                legend.box = "vertical",
                legend.box.margin = ggplot2::margin(0,0,0,-1,"cm"),
                plot.title=element_text(size=8,face="bold"))
        if(!is.null(mycols)) p <- p + scale_color_manual(values = mycols)
        if (!is.null(ellipse)) {
          if (identical(ellipse, TRUE)) p <- p + stat_ellipse(type="norm", linetype=2)
          else if (is.character(ellipse) && ellipse %in% names(pdata))
            p <- p + stat_ellipse(aes_string(group=ellipse, color=ellipse), type="norm", linetype=2)
        }
        if (!is.null(ID) && ID %in% names(pdata)) p <- p + ggrepel::geom_text_repel(aes_string(label = ID), size = 2)

        # 통계
        if (isTRUE(statistics)) {
          method_tag <- if (is.null(rand.eff)) "MiRKAT" else "MiRKAT-LMM"
          if (!is.null(facet) && facet %in% names(pdata)) {
            res <- .mirkat_by_facet(ps_obj = psIN.na, mvar = mvar, facet = facet, metric = metric, rand.eff = rand.eff)
            # 저장
            for (j in seq_len(nrow(res))) {
              flev <- res$facet_val[j]
              fn <- file.path(out_table,
                              sprintf("%s.%s.%s.%s.facet=%s%s.csv",
                                      method_tag, metric, project, mvar, flev,
                                      ifelse(is.null(rand.eff), "", paste0(".rand=", rand.eff))))
              utils::write.csv(res[j, , drop=FALSE], fn, row.names = FALSE)
            }
            # 주석(각 facet별 p)
            ann <- data.frame(
              x = -Inf, y = -Inf,
              label = sprintf("%s\n%s p=%s", metric, method_tag,
                              ifelse(is.na(res$p), "NA", formatC(res$p, format="f", digits=3))),
              stringsAsFactors = FALSE
            )
            ann[[facet]] <- res$facet_val

            ann[[facet]] <- factor(ann[, facet], levels = intersect(orders, unique(ann[, facet])))

            p <- p + ggplot2::geom_text(
              data = ann,
              mapping = aes(x = x, y = y, label = label),
              size = 3, hjust = -0.005, vjust = -0.3,
              lineheight = 0.95,
              inherit.aes = FALSE
            )
          } else {
            pval <- .mirkat_simple(ps_obj = psIN.na, mvar = mvar, metric = metric, rand.eff = rand.eff)
            fn <- file.path(out_table,
                            sprintf("%s.%s.%s.%s%s.csv",
                                    method_tag, metric, project, mvar,
                                    ifelse(is.null(rand.eff), "", paste0(".rand=", rand.eff))))
            utils::write.csv(data.frame(metric=metric, mvar=mvar, p=pval, method=method_tag),
                             fn, row.names = FALSE)

            p <- p + ggplot2::geom_text(
              data = data.frame(x=-Inf, y=-Inf,
                                label=.mk_label_block(metric, pval, method_tag)),
              aes(x=x,y=y,label=label),
              size=3, hjust=-0.005, vjust=-0.3, lineheight=0.95,
              inherit.aes = FALSE
            )
          }
        }

        if (!is.null(facet) && facet %in% names(pdata)) {
          ncol <- length(unique(pdata[,facet]))
          pdata[[facet]] <- factor(pdata[[facet]], intersect(orders, unique(pdata[, facet])))
          p <- p + facet_wrap(as.formula(sprintf("~ %s", facet)), scales="free_x", ncol = ncol)
        }

        p <- p +
          scale_x_continuous(expand = expansion(mult = c(0.03, 0.03))) +
          scale_y_continuous(expand = expansion(mult = c(0.06, 0.03))) +
          theme(panel.grid = element_blank(),
                legend.key = element_blank(),
                panel.background = element_rect(fill = "white", colour = "Black", linewidth = 0.7, linetype = "solid"),
                aspect.ratio = 1) +
          geom_vline(xintercept = 0, linewidth = 0.1) +
          geom_hline(yintercept = 0, linewidth = 0.1)

        print(p)
      }
    }
    dev.off()
  }
}
