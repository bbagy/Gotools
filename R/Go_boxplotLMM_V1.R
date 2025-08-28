
#' Generate Boxplots for Multiple Variables
#'
#' @importFrom lmerTest lmer
#' @importFrom emmeans emmeans contrast
#' @importFrom ggpubr stat_pvalue_manual
#' @importFrom stats anova as.formula complete.cases
#' @param df Data frame containing the data to be plotted.
#' @param cate.vars Categorical variables to be used for the boxplot's x-axis.
#' @param project Project name used for output file naming.
#' @param outcomes Numeric variables to be plotted on the y-axis.
#' @param orders Vector of ordered factor levels for the categorical variables.
#' @param mycols Custom color palette for the plots.
#' @param combination Number of group combinations to display in the boxplot.
#' @param ylim Y-axis limits for the plot.
#' @param title Title of the plot.
#' @param facet Optional variable for creating facetted plots.
#' @param paired Indicates if the data points are paired.
#' @param name Optional name for saving the plot.
#' @param ncol Number of columns for facet wrapping in the plot.
#' @param addnumber Boolean to add the number of samples in each group.
#' @param standardsize Boolean to control the size of the plot based on group size.
#' @param statistics Whether to perform statistical tests.
#' @param parametric Whether to use parametric tests.
#' @param star Whether to show significance levels as stars.
#' @param xangle Angle of x-axis labels.
#' @param cutoff Significance level for statistical tests.
#' @param height Height of the plot.
#' @param width Width of the plot.
#' @param plotCols Number of columns in the multiplot layout.
#' @param plotRows Number of rows in the multiplot layout.
#'
#' @details
#' This function creates boxplots for multiple categorical variables against one or more numeric outcomes. It supports faceting, custom color palettes, and statistical testing with options for parametric or non-parametric methods. The function can handle paired data and allows customization of plot size and layout.
#'
#' @return
#' A PDF file containing the boxplots.
#'
#' @examples
#' Go_boxplot(df = my_data, cate.vars = c("Group", "Condition"), project = "MyProject",
#'            outcomes = c("Variable1", "Variable2"), orders = c("Control", "Treatment"),
#'            height = 10, width = 8, plotCols = 2, plotRows = 1)
#'
#' @export

Go_boxplotLMM <- function(df, cate.vars, project, outcomes,
                       orders=NULL,
                       mycols=NULL,
                       combination=NULL,
                       ylim =NULL,
                       title= NULL,
                       facet= NULL,
                       paired=NULL,
                       name= NULL,
                       ncol=NULL,
                       addnumber=TRUE,
                       standardsize=TRUE,
                       statistics = TRUE,
                       parametric= FALSE,
                       star=TRUE,
                       xangle=90,
                       cutoff = 0.1,
                       height, width, plotCols, plotRows){

  suppressPackageStartupMessages({
    library(ggplot2)
    library(ggpubr)
    library(dplyr)
    library(lme4)
    library(lmerTest)
    library(emmeans)
  })

  emmeans::emm_options(lmer.df = "asymptotic")
  use_lmm_pairwise    <- !is.null(paired)  # paired 있으면 LMM만
  use_ggpubr_pairwise <-  is.null(paired)  # paired 없으면 ggpubr만
  .lmm_summary_mode   <- "compact"
  .lmm_adjust_method  <- "BH"   # 다중비교 보정(BH)

  # --------------------------------------------
  # LMM 보조 함수들 (최소 침습 추가)
  # --------------------------------------------

  # 1) LMM + emmeans 요약 콘솔 출력 (paired일 때만 호출)
  multiplot <- function(..., plotlist=NULL, file, cols=1, rows=1) {
    require(grid)
    # Make a list from the ... arguments and plotlist
    plots <- c(list(...), plotlist)
    numPlots = length(plots)

    i = 1
    while (i < numPlots) {
      numToPlot <- min(numPlots-i+1, cols*rows)
      # Make the panel
      # ncol: Number of columns of plots
      # nrow: Number of rows needed, calculated from # of cols
      layout <- matrix(seq(i, i+cols*rows-1), ncol = cols, nrow = rows, byrow=T)
      if (numToPlot==1) {
        print(plots[[i]])
      } else {
        # Set up the page
        grid.newpage()
        pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
        # Make each plot, in the correct location
        for (j in i:(i+numToPlot-1)) {
          # Get the i,j matrix positions of the regions that contain this subplot
          matchidx <- as.data.frame(which(layout == j, arr.ind = TRUE))
          print(plots[[j]], vp = viewport(layout.pos.row = matchidx$row,
                                          layout.pos.col = matchidx$col))
        }
      }
      i <- i+numToPlot
    }
  }


  .print_lmm_emm <- function(DAT, mvar, oc, id_col, emm_adjust = "none", digits = 4) {
    if (!all(c(mvar, oc, id_col) %in% names(DAT))) return(invisible(NULL))
    keep <- stats::complete.cases(DAT[, c(mvar, oc, id_col)])
    DAT  <- DAT[keep, , drop = FALSE]
    if (!nrow(DAT)) return(invisible(NULL))

    # (n=) 라벨 제거 후 모델에 투입
    DAT[[mvar]] <- as.character(DAT[[mvar]])
    DAT[[mvar]] <- sub(" \\(n=.*\\)$", "", DAT[[mvar]])
    DAT[[mvar]] <- droplevels(factor(DAT[[mvar]]))
    if (nlevels(DAT[[mvar]]) < 2) return(invisible(NULL))

    fm <- stats::as.formula(sprintf("%s ~ %s + (1|%s)", oc, mvar, id_col))
    mod <- try(lme4::lmer(fm, data = DAT), silent = TRUE)
    if (inherits(mod, "try-error")) {
      message(sprintf("[LMM] fit failed for %s ~ %s | id=%s", oc, mvar, id_col))
      return(invisible(NULL))
    }

    cat("\n==================== LMM summary ====================\n")
    cat(sprintf("Outcome: %s | Fixed: %s | Random: (1|%s)\n", oc, mvar, id_col))
    print(summary(mod), digits = digits)

    cat("\n-------------------- ANOVA(mod) ---------------------\n")
    suppressMessages(print(anova(mod), digits = digits))  # lmerTest 로드 시 Pr(>F) 포함

    emm <- try(emmeans::emmeans(mod, stats::as.formula(paste0("~ ", mvar))), silent = TRUE)
    if (!inherits(emm, "try-error")) {
      cat("\n================= emmeans summary ===================\n")
      print(summary(emm, infer = c(TRUE, TRUE)), digits = digits)

      contr <- try(emmeans::contrast(emm, method = "pairwise", adjust = emm_adjust), silent = TRUE)
      if (!inherits(contr, "try-error")) {
        cat("\n============= emmeans pairwise contrasts ============\n")
        print(summary(contr, infer = c(TRUE, TRUE)), digits = digits)
      }
    } else {
      cat("\n[emmeans] failed to compute marginal means.\n")
    }
    invisible(NULL)
  }
  # ---- dispatcher: 모드에 따라 compact/full/none 선택 ----
  .print_lmm <- function(DAT, mvar, oc, id_col, digits = 3) {
    if (identical(.lmm_summary_mode, "compact")) {
      .print_lmm_compact(DAT, mvar, oc, id_col,
                         adjust = .lmm_adjust_method,
                         digits = digits)
    } else if (identical(.lmm_summary_mode, "full")) {
      .print_lmm_emm(DAT, mvar, oc, id_col,
                     emm_adjust = "none",
                     digits = max(4, digits))
    } else {
      invisible(NULL)
    }
  }

  # ---- LMM compact printer: one concise table per plot ----
  # contrast_mode: "trt.vs.ctrl" (baseline vs others) | "pairwise" (모든 페어)
  .print_lmm_compact <- function(DAT, mvar, oc, id_col,
                                 adjust = "BH",
                                 digits = 3) {
    if (!all(c(mvar, oc, id_col) %in% names(DAT))) return(invisible(NULL))
    keep <- stats::complete.cases(DAT[, c(mvar, oc, id_col)])
    DAT  <- DAT[keep, , drop = FALSE]
    if (!nrow(DAT)) return(invisible(NULL))

    # (n=) 제거
    DAT[[mvar]] <- sub(" \\(n=.*\\)$", "", as.character(DAT[[mvar]]))
    DAT[[mvar]] <- droplevels(factor(DAT[[mvar]]))
    if (nlevels(DAT[[mvar]]) < 2) return(invisible(NULL))

    # LMM 적합
    fm  <- stats::as.formula(sprintf("%s ~ %s + (1|%s)", oc, mvar, id_col))
    mod <- try(lme4::lmer(fm, data = DAT), silent = TRUE)
    if (inherits(mod, "try-error")) return(invisible(NULL))

    # Omnibus p
    atab <- suppressMessages(stats::anova(mod))
    p_omni <- NA_real_
    if ("Pr(>F)" %in% colnames(atab)) {
      rn <- rownames(atab)
      hit <- which(rn == mvar)
      if (length(hit) >= 1) p_omni <- as.numeric(atab[hit[1], "Pr(>F)"])
    }

    # emmeans & 모든 pairwise 대비
    emm <- try(emmeans::emmeans(mod, stats::as.formula(paste0("~ ", mvar)), lmer.df = "asymptotic"), silent = TRUE)
    if (inherits(emm, "try-error")) return(invisible(NULL))
    contr <- try(emmeans::contrast(emm, method = "pairwise", adjust = "none"), silent = TRUE)
    if (inherits(contr, "try-error")) return(invisible(NULL))
    cdf <- as.data.frame(contr)
    parts <- strsplit(cdf$contrast, " - ", fixed = TRUE)
    g1 <- vapply(parts, `[`, character(1), 1)
    g2 <- vapply(parts, `[`, character(1), 2)
    cdf$padj <- stats::p.adjust(cdf$p.value, method = adjust)

    out <- data.frame(
      outcome  = oc,
      group1   = g1,
      group2   = g2,
      estimate = round(cdf$estimate, digits),
      p        = signif(cdf$p.value, digits),
      padj     = signif(cdf$padj, digits),
      stringsAsFactors = FALSE
    )

    cat(sprintf("\n[LMM compact] %s ~ %s + (1|%s) | omnibus p = %s | adjust=%s | mode=pairwise\n",
                oc, mvar, id_col,
                ifelse(is.na(p_omni), "NA", signif(p_omni, digits)),
                adjust))
    print(out, row.names = FALSE)
    invisible(NULL)
  }

  # 2) LMM pairwise p-value 주석 생성 (방향 무관 매칭 보강 + 실패시 윌콕슨 fallback)
  .lmm_pvals <- function(df, mvar, oc, id_col, comparisons, y_nudge = 1.08, digits = 3) {
    # 0) 입력/결측 가드
    need_cols <- c(mvar, oc, id_col)
    need_cols <- need_cols[!is.null(need_cols)]
    keep <- stats::complete.cases(df[, need_cols, drop = FALSE])
    df <- df[keep, , drop = FALSE]
    if (!nrow(df)) return(data.frame())

    # 1) (n=) 제거된 "원라벨" 기준 factor 구성 (모델/emmeans는 항상 원라벨로)
    raw_vec <- sub(" \\(n=.*\\)$", "", as.character(df[[mvar]]))
    df[[mvar]] <- droplevels(factor(raw_vec))
    if (nlevels(df[[mvar]]) < 2) return(data.frame())
    lvl_order <- levels(df[[mvar]])

    # 2) 비교쌍 사전 정리: 원라벨 기준으로 정리 + 최소 표본/ID 가드
    #    (비교쌍은 방향 무관하게 키를 만들어 매칭)
    .pair_key <- function(a, b, lvl = NULL) {
      if (!is.null(lvl)) {
        ia <- match(a, lvl); ib <- match(b, lvl)
        ord <- order(c(ia, ib), na.last = TRUE)
        paste(c(a, b)[ord], collapse = "__")
      } else paste(sort(c(a, b)), collapse = "__")
    }

    .keep_pair <- function(a, b, d, id_col = NULL) {
      raw <- as.character(d[[mvar]])
      nA <- sum(raw == a); nB <- sum(raw == b)
      ok <- (nA >= 2 && nB >= 2)             # 쌍별 최소 표본수
      if (!is.null(id_col)) {
        idsA <- unique(d[raw == a, id_col])
        idsB <- unique(d[raw == b, id_col])
        # 두 군에 걸쳐 최소 2명 이상 존재 (너무 빡세지 않게 완화)
        ok <- ok && length(unique(na.omit(c(idsA, idsB)))) >= 2
      }
      ok
    }

    comp_raw <- lapply(comparisons, function(x) {
      a <- sub(" \\(n=.*\\)$", "", x[1])
      b <- sub(" \\(n=.*\\)$", "", x[2])
      c(a, b)
    })
    comp_raw <- Filter(function(z) all(z %in% lvl_order) && .keep_pair(z[1], z[2], df, id_col), comp_raw)
    if (!length(comp_raw)) return(data.frame())

    # 3) LMM 적합
    fm  <- stats::as.formula(sprintf("%s ~ %s + (1|%s)", oc, mvar, id_col))
    mod <- try(lme4::lmer(fm, data = df), silent = TRUE)
    if (inherits(mod, "try-error")) return(data.frame())

    # 4) emmeans (asymptotic df로 NA 감소)
    emm <- try(emmeans::emmeans(mod, stats::as.formula(paste0("~ ", mvar)), lmer.df = "asymptotic"), silent = TRUE)
    if (inherits(emm, "try-error")) return(data.frame())

    contr <- try(emmeans::contrast(emm, method = "pairwise", adjust = "none"), silent = TRUE)
    if (inherits(contr, "try-error")) return(data.frame())

    # 5) emmeans 결과를 (a,b) 쌍 키로 매칭
    contr_df <- as.data.frame(contr)
    parse_con <- function(s) strsplit(s, " - ", fixed = TRUE)[[1]]
    contr_df$g1 <- vapply(contr_df$contrast, function(z) parse_con(z)[1], character(1))
    contr_df$g2 <- vapply(contr_df$contrast, function(z) parse_con(z)[2], character(1))
    contr_df$key <- mapply(.pair_key, contr_df$g1, contr_df$g2, MoreArgs = list(lvl = lvl_order))

    want <- data.frame(
      g1 = vapply(comp_raw, `[`, character(1), 1),
      g2 = vapply(comp_raw, `[`, character(1), 2)
    )
    want$key <- mapply(.pair_key, want$g1, want$g2, MoreArgs = list(lvl = lvl_order))

    out <- merge(contr_df[, c("key", "g1", "g2", "p.value")], unique(want["key"]), by = "key", all.y = TRUE)
    names(out)[names(out) == "p.value"] <- "p"

    # 6) NA 보수: 해당 쌍만 윌콕슨으로 대체 (paired 가능하면 paired=TRUE)
    if (any(is.na(out$p))) {
      get_p_fallback <- function(a, b) {
        raw <- as.character(df[[mvar]])
        x <- df[raw == a, oc]; y <- df[raw == b, oc]
        do_paired <- FALSE
        if (!is.null(id_col)) {
          idsA <- unique(df[raw == a, id_col])
          idsB <- unique(df[raw == b, id_col])
          if (length(intersect(idsA, idsB)) >= 2) do_paired <- TRUE
        }
        p <- try(stats::wilcox.test(x, y, paired = do_paired, exact = FALSE)$p.value, silent = TRUE)
        if (inherits(p, "try-error")) NA_real_ else p
      }
      for (i in which(is.na(out$p))) {
        out$p[i] <- get_p_fallback(out$g1[i], out$g2[i])
      }
    }

    # 7) 주석용 좌표/라벨
    ymax <- max(df[[oc]], na.rm = TRUE) * y_nudge
    out$y.position <- ymax * (1 + 0.06 * (seq_len(nrow(out)) - 1))
    out$label <- ifelse(is.na(out$p), "NA", sprintf(paste0("%.", digits, "f"), out$p))
    names(out)[names(out) %in% c("g1", "g2")] <- c("group1", "group2")

    # (n=) 라벨 매핑은 바깥에서 .map_ann_levels()로 처리
    out
  }

  # 3) (n=) 라벨로 바뀐 축 레벨에 p주석 레벨 매핑
  .map_ann_levels <- function(ann, old2new) {
    if (!nrow(ann)) return(ann)
    ann$group1 <- ifelse(ann$group1 %in% names(old2new), old2new[ann$group1], ann$group1)
    ann$group2 <- ifelse(ann$group2 %in% names(old2new), old2new[ann$group2], ann$group2)
    ann
  }

  # 4) (선택) LMM omnibus p-value (제목에 쓰고자 할 때 사용 가능)
  .lmm_omnibus <- function(df, mvar, oc, id_col) {
    keep <- stats::complete.cases(df[, c(mvar, oc, id_col)])
    df <- df[keep, , drop = FALSE]
    if (!nrow(df)) return(list(name="LMM", pval=NA_real_))

    df[[mvar]] <- as.character(df[[mvar]])
    df[[mvar]] <- sub(" \\(n=.*\\)$", "", df[[mvar]])
    df[[mvar]] <- droplevels(factor(df[[mvar]]))
    if (nlevels(df[[mvar]]) < 2) return(list(name="LMM", pval=NA_real_))

    fm  <- stats::as.formula(sprintf("%s ~ %s + (1|%s)", oc, mvar, id_col))
    mod <- lmerTest::lmer(fm, data = df)
    atab <- suppressMessages(stats::anova(mod))
    rn <- rownames(atab)
    hit <- which(rn == mvar)
    if (length(hit) == 0) hit <- grep(sprintf("^%s$", gsub("([\\W])", "\\\\\\1", mvar)), rn)
    p <- if (length(hit)>=1 && "Pr(>F)" %in% colnames(atab)) atab[hit[1], "Pr(>F)"] else NA_real_
    list(name="LMM", pval=ifelse(is.na(p), NA_real_, round(as.numeric(p), 4)))
  }



  if(!is.null(dev.list())) dev.off()

  # out dir
  out <- file.path(sprintf("%s_%s",project, format(Sys.Date(), "%y%m%d")))
  if(!file_test("-d", out)) dir.create(out)
  out_path <- file.path(sprintf("%s_%s/pdf",project, format(Sys.Date(), "%y%m%d")))
  if(!file_test("-d", out_path)) dir.create(out_path)
  set.seed(151)

  # out file
  if (class(name) == "function"){
    name <- NULL
  }

  tt <- try(mycols,T)
  if(class(tt) == "try-error"){
    print("mycols is not defined.")
    mycols <- NULL
  }

  tt <- try(orders,T)
  if(class(tt) == "try-error"){
    print("orders is not defined.")
    orders <- NULL
  }

  print("box.tickness1")

  # plot design
  if (height*width <= 6){
    dot.size = 0.7
    box.tickness = 0.3
  }else if (height*width > 6 & height*width < 10){
    dot.size = 1
    box.tickness = 0.4
  }else{
    dot.size = 1.5
    box.tickness = 0.5
  }

  pdf(sprintf("%s/box.%s.%s%s%s%s%s.pdf", out_path,
              project,
              ifelse(is.null(facet), "", paste(facet, ".", sep = "")),
              ifelse(is.null(paired), "", paste("(paired=",paired, ").", sep = "")),
              ifelse(is.null(combination), "", paste("(cbn=",combination, ").", sep = "")),
              ifelse(is.null(name), "", paste(name, ".", sep = "")),
              format(Sys.Date(), "%y%m%d")), height = height, width = width)

  # plot
  plotlist <- list()
  for (mvar in cate.vars) {
    if (length(unique(df[,mvar])) < 2){
      next
    }

    if (length(facet) >= 1){
      if (facet == mvar){
        next
      }
    } else {}

    # remove Na
    print("Control NA")
    df <- data.frame(df)
    df[,mvar] <- as.character(df[,mvar]);df[,mvar]
    df[,mvar][df[,mvar]==""] <- "NA";df[,mvar]
    df.na <- subset(df, df[,mvar] != "NA");df.na[,mvar]  # subset 를 사용한 NA 삭제
    df.na[,mvar] <- as.factor(df.na[,mvar]);df.na[,mvar]

    df.na[,mvar] <- factor(df.na[,mvar], levels = intersect(orders, df.na[,mvar]))

    # Add number of samples in the group
    print("Add sample number informations")
    if(addnumber==TRUE){
      renamed_levels <- as.character(levels(df.na[,mvar]));renamed_levels
      oldNames <- unique(df.na[,mvar]);oldNames
      if (length(renamed_levels) == 0) {
        renamed_levels <- oldNames
      }
      for (name in oldNames) {
        total <- length(which(df.na[,mvar] == name));total
        new_n <- paste(name, " (n=", total, ")", sep="");new_n
        levels(df.na[[mvar]])[levels(df.na[[mvar]])== name] <- new_n
        renamed_levels <- replace(renamed_levels, renamed_levels == name, new_n);renamed_levels
      }
    }else{
      df.na <- df.na
    }

    print(sprintf("##-- %s (total without NA: %s/%s) --##",
                  mvar, dim(df.na)[1], dim(df)[1]))

    if (length(unique(df.na[,mvar])) ==1) {
      next
    }

    summary.df.na <- summary(df.na[,mvar])



    #------------------------------#
    # for group combination or not #
    #------------------------------#

    if (!is.null(combination)){
      print(sprintf("Combination n=%s", combination))
      group.cbn <- combn(x = levels(df.na[,mvar]), m = combination)

      group_comparisons <- {}
      for(i in 1:ncol(group.cbn)){
        x <- group.cbn[,i]
        group_comparisons[[i]] <- x
      };group_comparisons

      print(1)
      for (i in seq_along(group_comparisons)) {
        group.combination <- unlist(group_comparisons[i])

        # 가드: 선택 가능한 길이로 제한
        k <- min(combination, length(group.combination))
        if (k < 2) next

        baseline <- group.combination[1]
        sel <- group.combination[seq_len(k)]

        df.cbn <- subset(df.na, df.na[[mvar]] %in% sel)

        lvl_new <- c(baseline, setdiff(sel, baseline))
        df.cbn[[mvar]] <- factor(df.cbn[[mvar]], levels = lvl_new)

        # baseline vs others
        if (k == 2) {
          my_comparisons <- list(c(lvl_new[1], lvl_new[2]))
        } else {
          my_comparisons <- lapply(lvl_new[-1], function(g) c(lvl_new[1], g))
        }

        for(oc in outcomes){
          df_full <- df.na
          df_full[[mvar]] <- sub(" \\(n=.*\\)$", "", as.character(df_full[[mvar]]))
          df_full[[mvar]] <- droplevels(factor(df_full[[mvar]]))
          if (!is.null(orders) && length(orders) > 0) {
            df_full[[mvar]] <- factor(df_full[[mvar]], levels = intersect(orders, levels(df_full[[mvar]])))
          }
          # remove NA for facet
          if (length(facet) >= 1) {
            for (fc in facet){
              df.cbn[,fc] <- as.character(df.cbn[,fc]);df.cbn[,fc]
              df.cbn[,fc][df.cbn[,fc] == ""] <- "NA"
              df.cbn.sel <- df.cbn[!is.na(df.cbn[,fc]), ]
              df.cbn <- df.cbn.sel
              # facet or not
              df.cbn[,fc] <- factor(df.cbn[,fc], levels = orders)
            }
          }



          print(oc)

          # check statistics method (원본 유지)
          if (statistics){
            if (parametric){
              if (nlevels(factor(df.cbn[,mvar])) > 2) {
                test <- aov(as.formula(sprintf("%s ~ %s", oc, mvar)), df.cbn)
                pval <- round(summary(test)[[1]][["Pr(>F)"]][1],4)
                test.name <- "ANOVA"
                testmethod <-  "t.test"
              } else {
                testmethod <-  "t.test"
                pval <- NULL
                test.name <- "Pairwise T-Test"
              }
            }else{
              if (nlevels(factor(df.cbn[,mvar])) > 2) {
                test <- kruskal.test(as.formula(sprintf("%s ~ %s", oc, mvar)), df.cbn)
                pval <- round(test$p.value, 4)
                test.name <- "KW"
                testmethod <- "wilcox.test"
              } else {
                testmethod <- "wilcox.test"
                pval <- NULL
                test.name <- "Pairwise Wilcoxon"
              }
            }
          }else{
            test.name<-NULL
            pval <- NULL
          }

          # ----- PATCH: paired이면 제목용 전역 테스트를 LMM으로 교체 -----
          if (!is.null(paired)) {
            om <- .lmm_omnibus(df = df.cbn, mvar = mvar, oc = oc, id_col = paired)
            test.name <- om$name        # "LMM"
            pval <- om$pval             # LMM omnibus p-value
            if (!is.na(pval)) pval <- round(pval, 3)  # 3자리로 표시 원하면
          }
          # ----------------------------------------------------------------


          # x축은 반드시 *_lab 사용
          xvar_lab <- paste0(mvar, "_lab")




          p1 <- ggplot(df.cbn, aes_string(x=mvar, y=oc))  + labs(y=oc, x=NULL) +
            theme(strip.background = element_blank()) +
            theme(text=element_text(size=8), axis.text.x=element_text(angle=xangle,hjust=1,vjust=0.5,size=8),
                  plot.title=element_text(size=8))



          if (!is.null(title)) {
            p1 <- p1 + ggtitle(sprintf("%s%s%s%s", title,
                                       ifelse(is.null(test.name), "", paste("\n",test.name, " ")),
                                       ifelse(is.null(pval), "", paste("p=", " ")),
                                       ifelse(is.null(pval), "", paste(pval, " "))))
          } else{
            p1 <- p1 + ggtitle(sprintf("%s%s%s%s", mvar,
                                       ifelse(is.null(test.name), "", paste("\n",test.name, " ")),
                                       ifelse(is.null(pval), "", paste("p=", " ")),
                                       ifelse(is.null(pval), "", paste(pval, " "))))
          }

          # control statistic on the plot (원본 유지)
          if (use_ggpubr_pairwise && statistics) {
            if (test.name %in% c("KW", "ANOVA")) {
              # pval이 NA일 수도 있으니 isTRUE()로 방어
              if (isTRUE(pval < cutoff)) {
                if (star) {
                  p1 <- p1 + stat_compare_means(
                    method = testmethod, label = "p.signif",
                    comparisons = my_comparisons, hide.ns = TRUE, size = 3
                  )
                } else {
                  p1 <- p1 + stat_compare_means(
                    method = testmethod, label = "p.format",
                    comparisons = my_comparisons, size = 2
                  )
                }
              }
            } else if (testmethod %in% c("wilcox.test", "t.test")) {
              if (star) {
                p1 <- p1 + stat_compare_means(
                  method = testmethod, label = "p.signif",
                  comparisons = my_comparisons, hide.ns = TRUE, size = 3
                )
              } else {
                p1 <- p1 + stat_compare_means(
                  method = testmethod, label = "p.format",
                  comparisons = my_comparisons, size = 2
                )
              }
            }
          }

          if(!is.null(ylim)){
            if(oc == "Chao1"){
              p1 = p1
            }else{
              p1 = p1 + ylim(ylim[1] , ylim[2])
            }
          }

          # paired plot type (원본 유지 + LMM 주석/요약 추가)
          if (!is.null(paired)) {

            if(!is.null(mycols)){
              p1 <- p1 + scale_color_manual(values = mycols)
            }else{
              p1 <- p1
            }

            p1 = p1 + geom_boxplot(aes_string(colour=mvar),outlier.shape = NA,lwd=box.tickness)  + theme(legend.position="none")
            p1 = p1 + geom_point(aes_string(colour=mvar,group=paired),alpha = 0.8, size = dot.size, position = position_dodge(0.3),show.legend = F)
            p1 = p1 + geom_line(aes_string(group=paired), color="grey50", size=0.3,position = position_dodge(0.3))
            p1 = p1 + theme(legend.title = element_blank(), legend.position="bottom", legend.justification="left",legend.box.margin = ggplot2::margin(0,0,0,-1,"cm"))

            # ---- (추가) LMM pairwise 주석 + 요약 출력 ----
            # 통계용 원본 레벨로 비교쌍 구성 (baseline-only 규칙 유지)
            tmp_lmm <- df.cbn
            tmp_lmm[[mvar]] <- sub(" \\(n=.*\\)$", "", as.character(tmp_lmm[[mvar]]))
            tmp_lmm[[mvar]] <- droplevels(factor(tmp_lmm[[mvar]]))
            tmp_lmm[[mvar]] <- factor(tmp_lmm[[mvar]],
                                      intersect(orders, tmp_lmm[,mvar]))

            lv_raw <- levels(tmp_lmm[[mvar]])

            if (length(lv_raw) >= 2) {
              # 2) 비교쌍 구성: combination이면 baseline(첫 레벨) vs others, 아니면 모든 페어
              if (!is.null(combination) && combination != 2) {
                my_comp0 <- lapply(lv_raw[-1], function(g) c(lv_raw[1], g))
              } else {
                cbn0 <- combn(lv_raw, 2)
                my_comp0 <- lapply(seq_len(ncol(cbn0)), function(k) cbn0[, k])
              }

              # 3) LMM pairwise p 추출
              ann <- .lmm_pvals(df = tmp_lmm, mvar = mvar, oc = oc, id_col = paired,
                                comparisons = my_comp0, y_nudge = 1.08, digits = 3)
              ## ----- [PVAL FORMAT PATCH for combination == 전체그룹수] -----
              ## ann: .lmm_pvals(...) 결과 (group1, group2, p, y.position, label ...)
              ## dfX: 지금 그리는 데이터프레임 (combination이면 df.cbn, 아니면 df.na)
              ## mvar, star, cutoff 변수는 기존 그대로 사용

              # 0) 전체 그룹 수 = combination 인지 확인
              n_groups <- nlevels(df.cbn[[mvar]])
              if (!is.null(combination) && combination == n_groups && nrow(ann) > 0) {

                # 1) 다중비교 보정 (BH). NA는 그대로 유지되므로 먼저 제거 후 보정하고 다시 병합
                keep_idx <- which(!is.na(ann$p))
                if (length(keep_idx)) {
                  padj <- stats::p.adjust(ann$p[keep_idx], method = "BH")
                  ann$padj <- NA_real_
                  ann$padj[keep_idx] <- padj
                } else {
                  ann$padj <- ann$p
                }

                # 2) 표기 라벨: star 또는 숫자 포맷
                if (isTRUE(star)) {
                  # cutpoints는 필요 시 조정 가능
                  ann$label <- ifelse(is.na(ann$p), "NA",
                                      ifelse(ann$p < 0.001, "***",
                                             ifelse(ann$p < 0.01, "**",
                                                    ifelse(ann$p < 0.05, "*",
                                                           "ns"))))
                } else {
                  ann$label <- ifelse(is.na(ann$p), "NA",
                                      sprintf("%.3g", ann$p))
                }

                # 3) cutoff 적용: 유의한 것만 남김 (겹침 방지)
                ann <- ann[!is.na(ann$p) & ann$p < cutoff, , drop = FALSE]

                # 유의한 게 하나도 없으면 주석 스킵
                if (nrow(ann) > 0) {
                  # 4) y 위치: 현재 데이터 범위 기반으로 자동 계단 배치(겹침 최소화)
                  ymax <- max(df.cbn[[oc]], na.rm = TRUE)
                  base <- ymax * 1.08
                  step <- 0.06 * ymax
                  ann <- ann[order(ann$padj, decreasing = FALSE), , drop = FALSE]  # 더 유의한 것 위쪽
                  ann$y.position <- base + step * seq_len(nrow(ann))

                  # 5) (좌표 방식) x축 인덱스 매핑
                  x_levels   <- levels(df.cbn[[mvar]])                  # 라벨 유무 무관
                  raw_from_x <- sub(" \\(n=.*\\)$", "", x_levels)    # 원라벨로 통일
                  raw2pos    <- setNames(seq_along(x_levels), raw_from_x)

                  ann$xmin <- unname(raw2pos[ann$group1])
                  ann$xmax <- unname(raw2pos[ann$group2])
                  ann <- ann[!is.na(ann$xmin) & !is.na(ann$xmax), , drop = FALSE]

                  if (nrow(ann) > 0) {
                    p1 <- p1 + ggpubr::stat_pvalue_manual(
                      ann,
                      xmin       = "xmin",
                      xmax       = "xmax",
                      y.position = "y.position",
                      label      = "label",
                      tip.length = 0.01,
                      size       = 3
                    )
                  }
                }
              }
              ## ----- [END PATCH] -----

              # 4) 좌표 기반으로 p 주석 올리기 (문자열 매칭 X → x축 NA 완전 차단)
              if (nrow(ann) > 0) {
                x_levels   <- levels(df.cbn[[mvar]])                    # 현재 x축 레벨(라벨 포함/미포함 상관없음)
                raw_from_x <- sub(" \\(n=.*\\)$", "", x_levels)         # 원라벨 추출
                raw2pos    <- setNames(seq_along(x_levels), raw_from_x) # 원라벨 -> 위치

                ann$xmin <- unname(raw2pos[ann$group1])
                ann$xmax <- unname(raw2pos[ann$group2])
                ann <- ann[!is.na(ann$xmin) & !is.na(ann$xmax), , drop = FALSE]

                if (nrow(ann)) {
                  p1 <- p1 + ggpubr::stat_pvalue_manual(
                    ann,
                    xmin       = "xmin",
                    xmax       = "xmax",
                    y.position = "y.position",
                    label      = "label",
                    tip.length = 0.01,
                    size       = 3
                  )
                }
              }

              # 5) 콘솔 요약
              .print_lmm(DAT = df_full, mvar = mvar, oc = oc, id_col = paired, digits = 3)
            }

          }  else{

            # count or table for number of variable (원본)
            if (max(table(df.cbn[,mvar])) > 150){

              if(!is.null(mycols)){
                p1 <- p1 + scale_fill_manual(values = mycols)
              }else{
                p1 <- p1
              }
              p1 = p1 + geom_boxplot(aes_string(fill=mvar),outlier.shape = NULL,lwd=box.tickness)   + theme(legend.position="none")
            } else {

              if(!is.null(mycols)){
                p1 <- p1 + scale_color_manual(values = mycols)
              }else{
                p1 <- p1
              }

              p1 = p1 + geom_boxplot(aes_string(colour=mvar),outlier.shape = NA,lwd=box.tickness) + theme(legend.position="none")
              p1 = p1 + geom_jitter(aes_string(colour=mvar),shape=16, alpha = 0.8, size = dot.size, position=position_jitter(0.2))

            }
          }

          # facet (원본)
          if (length(facet) >= 1) {
            if(is.null(ncol)){
              ncol <- length(unique(df[,facet]))
            }
            p1 = p1 + facet_wrap(as.formula(sprintf("~ %s" , paste(setdiff(facet, "SocpleType"), collapse="+"))),
                                 scales="free_x", ncol = ncol)
            p1 = p1 + guides(color = "none", size = "none", shape= "none")
          } else {
            p1 = p1 + guides(color = "none", size = "none", shape= "none")
          }

          # plot size ratio (원본)
          if (length(unique(df.cbn[,mvar])) < 5){
            if(standardsize==TRUE){
              num.subgroup <- length(unique(df.cbn[,mvar]))*0.1
            }else{
              num.subgroup <- 0.9
            }
          }else{
            num.subgroup <- length(unique(df.cbn[,mvar]))*0.1
          }

          p1  <- p1  + theme(panel.grid = element_blank(),
                             panel.background = element_rect(fill = "white", colour = "Black",size = 0.5, linetype = "solid"),
                             aspect.ratio = 1/num.subgroup)

          plotlist[[length(plotlist)+1]] <- p1
        }
      }
    }else{
      # ------------ combination 미사용 (원본) ------------
      print("Check combination for statistics")
      cbn <- combn(x = levels(df.na[,mvar]), m = 2)

      my_comparisons <- {}
      for(i in 1:ncol(cbn)){
        x <- cbn[,i]
        my_comparisons[[i]] <- x
      };my_comparisons

      for(oc in outcomes){
        df_full <- df.na
        df_full[[mvar]] <- sub(" \\(n=.*\\)$", "", as.character(df_full[[mvar]]))
        df_full[[mvar]] <- droplevels(factor(df_full[[mvar]]))
        if (!is.null(orders) && length(orders) > 0) {
          df_full[[mvar]] <- factor(df_full[[mvar]], levels = intersect(orders, levels(df_full[[mvar]])))
        }
        # remove NA for facet
        if (!is.null(facet)) {
          for (fc in facet){
            df.na[,fc] <- as.character(df.na[,fc]);df.na[,fc]
            df.na[,fc][df.na[,fc] == ""] <- "NA"
            df.na.sel <- df.na[!is.na(df.na[,fc]), ]
            df.na <- df.na.sel
            # facet or not
            df.na[,fc] <- factor(df.na[,fc], levels = orders)
          }
        }

        print(oc)

        if (statistics){
          if (parametric){
            if (nlevels(factor(df.na[,mvar])) > 2) {
              test <- aov(as.formula(sprintf("%s ~ %s", oc, mvar)), df.na)
              pval <- round(summary(test)[[1]][["Pr(>F)"]][1],4)
              test.name <- "ANOVA"
              testmethod <-  "t.test"
            } else {
              testmethod <-  "t.test"
              pval <- NULL
              test.name <- "Pairwise T-Test"
            }
          }else{
            if (nlevels(factor(df.na[,mvar])) > 2) {
              test <- kruskal.test(as.formula(sprintf("%s ~ %s", oc, mvar)), df.na)
              pval <- round(test$p.value, 4)
              test.name <- "KW"
              testmethod <- "wilcox.test"
            } else {
              testmethod <- "wilcox.test"
              pval <- NULL
              test.name <- "Pairwise Wilcoxon"
            }
          }
        }else{
          test.name<-NULL
          pval <- NULL
        }

        # ----- PATCH: paired이면 제목용 전역 테스트를 LMM으로 교체 -----
        if (!is.null(paired)) {
          om <- .lmm_omnibus(df = df.na, mvar = mvar, oc = oc, id_col = paired)
          test.name <- om$name
          pval <- om$pval
          if (!is.na(pval)) pval <- round(pval, 3)
        }
        # ----------------------------------------------------------------
        p1 <- ggplot(df.na, aes_string(x=mvar, y=oc))  + labs(y=oc, x=NULL) +
          theme(strip.background = element_blank()) +
          theme(text=element_text(size=8), axis.text.x=element_text(angle=xangle,hjust=1,vjust=0.5,size=8),
                plot.title=element_text(size=8))


        # paired plot type (원본 유지 + LMM 주석/요약 추가)
        if (!is.null(paired)) {

          if(!is.null(mycols)){
            p1 <- p1 + scale_color_manual(values = mycols)
          }else{
            p1 <- p1
          }

          p1 = p1 + geom_boxplot(aes_string(colour=mvar),outlier.shape = NA,lwd=box.tickness)  + theme(legend.position="none")
          p1 = p1 + geom_point(aes_string(colour=mvar,group=paired),alpha = 0.8, size = dot.size, position = position_dodge(0.3), show.legend = F)
          p1 = p1 + geom_line(aes_string(group=paired), color="grey50", size=0.3,position = position_dodge(0.3))
          p1 = p1 + theme(legend.title = element_blank(), legend.position="bottom", legend.justification="left",legend.box.margin = ggplot2::margin(0,0,0,-1,"cm"))

          # ---- (추가) LMM pairwise 주석 + 요약 출력 ----
          tmp_lmm <- df.na
          tmp_lmm[[mvar]] <- sub(" \\(n=.*\\)$", "", as.character(tmp_lmm[[mvar]]))
          tmp_lmm[[mvar]] <- droplevels(factor(tmp_lmm[[mvar]]))

          tmp_lmm[[mvar]] <- factor(tmp_lmm[[mvar]],
                                    intersect(orders, tmp_lmm[,mvar]))
          lv_raw <- levels(tmp_lmm[[mvar]])



          if (length(lv_raw) >= 2) {
            # 2) 비교쌍 구성: non-combination에서는 모든 페어
            cbn0 <- combn(lv_raw, 2)
            my_comp0 <- lapply(seq_len(ncol(cbn0)), function(k) cbn0[, k])

            # 3) LMM pairwise p 추출
            ann <- .lmm_pvals(df = tmp_lmm, mvar = mvar, oc = oc, id_col = paired,
                              comparisons = my_comp0, y_nudge = 1.08, digits = 3)

            ## ----- [PVAL FORMAT PATCH for combination == 전체그룹수] -----
            ## ann: .lmm_pvals(...) 결과 (group1, group2, p, y.position, label ...)
            ## dfX: 지금 그리는 데이터프레임 (combination이면 df.cbn, 아니면 df.na)
            ## mvar, star, cutoff 변수는 기존 그대로 사용

            # 0) 전체 그룹 수 = combination 인지 확인
            n_groups <- nlevels(df.na[[mvar]])
            if (!is.null(combination) && combination == n_groups && nrow(ann) > 0) {

              # 1) 다중비교 보정 (BH). NA는 그대로 유지되므로 먼저 제거 후 보정하고 다시 병합
              keep_idx <- which(!is.na(ann$p))
              if (length(keep_idx)) {
                padj <- stats::p.adjust(ann$p[keep_idx], method = "BH")
                ann$padj <- NA_real_
                ann$padj[keep_idx] <- padj
              } else {
                ann$padj <- ann$p
              }

              # 2) 표기 라벨: star 또는 숫자 포맷
              if (isTRUE(star)) {
                # cutpoints는 필요 시 조정 가능
                ann$label <- ifelse(is.na(ann$p), "NA",
                                    ifelse(ann$p < 0.001, "***",
                                           ifelse(ann$p < 0.01, "**",
                                                  ifelse(ann$p < 0.05, "*",
                                                         "ns"))))
              } else {
                ann$label <- ifelse(is.na(ann$p), "NA",
                                    sprintf("%.3g", ann$p))
              }

              # 3) cutoff 적용: 유의한 것만 남김 (겹침 방지)
              ann <- ann[!is.na(ann$p) & ann$p < cutoff, , drop = FALSE]

              # 유의한 게 하나도 없으면 주석 스킵
              if (nrow(ann) > 0) {
                # 4) y 위치: 현재 데이터 범위 기반으로 자동 계단 배치(겹침 최소화)
                ymax <- max(df.na[[oc]], na.rm = TRUE)
                base <- ymax * 1.08
                step <- 0.06 * ymax
                ann <- ann[order(ann$padj, decreasing = FALSE), , drop = FALSE]  # 더 유의한 것 위쪽
                ann$y.position <- base + step * seq_len(nrow(ann))

                # 5) (좌표 방식) x축 인덱스 매핑
                x_levels   <- levels(df.na[[mvar]])                  # 라벨 유무 무관
                raw_from_x <- sub(" \\(n=.*\\)$", "", x_levels)    # 원라벨로 통일
                raw2pos    <- setNames(seq_along(x_levels), raw_from_x)

                ann$xmin <- unname(raw2pos[ann$group1])
                ann$xmax <- unname(raw2pos[ann$group2])
                ann <- ann[!is.na(ann$xmin) & !is.na(ann$xmax), , drop = FALSE]

                if (nrow(ann) > 0) {
                  p1 <- p1 + ggpubr::stat_pvalue_manual(
                    ann,
                    xmin       = "xmin",
                    xmax       = "xmax",
                    y.position = "y.position",
                    label      = "label",
                    tip.length = 0.01,
                    size       = 3
                  )
                }
              }
            }
            ## ----- [END PATCH] -----

            # 4) 좌표 기반으로 p 주석 올리기 (문자열 매칭 X → x축 NA 완전 차단)
            if (nrow(ann) > 0) {
              x_levels   <- levels(df.na[[mvar]])
              raw_from_x <- sub(" \\(n=.*\\)$", "", x_levels)
              raw2pos    <- setNames(seq_along(x_levels), raw_from_x)

              ann$xmin <- unname(raw2pos[ann$group1])
              ann$xmax <- unname(raw2pos[ann$group2])
              ann <- ann[!is.na(ann$xmin) & !is.na(ann$xmax), , drop = FALSE]

              if (nrow(ann)) {
                p1 <- p1 + ggpubr::stat_pvalue_manual(
                  ann,
                  xmin       = "xmin",
                  xmax       = "xmax",
                  y.position = "y.position",
                  label      = "label",
                  tip.length = 0.01,
                  size       = 3
                )
              }
            }

            # 5) 콘솔 요약
            .print_lmm(DAT = df_full, mvar = mvar, oc = oc, id_col = paired, digits = 3)
          }

        } else{
          # count or table for number of variable
          if (max(table(df.na[,mvar])) > 150){

            if(!is.null(mycols)){
              p1 <- p1 + scale_fill_manual(values = mycols)
            }else{
              p1 <- p1
            }
            p1 = p1 + geom_boxplot(aes_string(fill=mvar),outlier.shape = NULL,lwd=box.tickness)   + theme(legend.position="none")

          } else {

            if(!is.null(mycols)){
              p1 <- p1 + scale_color_manual(values = mycols)
            }else{
              p1 <- p1
            }

            p1 = p1 + geom_boxplot(aes_string(colour=mvar),outlier.shape = NA,lwd=box.tickness) + theme(legend.position="none")
            p1 = p1 + geom_jitter(aes_string(colour=mvar),shape=16, alpha = 0.8, size = dot.size, position=position_jitter(0.2))

          }
        }

        # control statistic on the plot (원본 유지)
        if (use_ggpubr_pairwise && statistics) {
          # 지금 플롯에 실제로 쓰이는 데이터프레임 선택
          dfX <- if (!is.null(combination)) df.cbn else df.na

          # x축 레벨(라벨)과 원라벨
          xlab <- levels(dfX[[mvar]])                    # 예: "Pre (n=26)", "wk1_Post (n=50)" ...
          xraw <- sub(" \\(n=.*\\)$", "", xlab)          # 예: "Pre", "wk1_Post" ...

          # baseline: orders[1] 우선, 없으면 x축 첫 항목
          if (!is.null(orders) && length(orders) > 0) {
            cand <- orders[orders %in% xraw]
            base_raw <- if (length(cand) >= 1) cand[1] else xraw[1]
          } else {
            base_raw <- xraw[1]
          }
          base_lab <- xlab[match(base_raw, xraw)]

          # 비교쌍: baseline vs others (라벨 기준으로 만들어 ggpubr에 넘김)
          if (length(xlab) >= 3) {
            my_comparisons_lab <- lapply(setdiff(xlab, base_lab), function(g) c(base_lab, g))
          } else {
            # 레벨이 2개면 그대로 한 쌍
            my_comparisons_lab <- list(xlab)
          }

          # 라벨 모드
          label_use <- if (isTRUE(star)) "p.signif" else "p.format"
          size_use  <- if (isTRUE(star)) 3 else 2
          hide_ns_use <- isTRUE(star)

          # KW/ANOVA는 전역 p가 cutoff 미만일 때만 pairwise 표시
          if (test.name %in% c("KW", "ANOVA")) {
            if (!is.null(pval) && pval < cutoff) {
              p1 <- p1 + ggpubr::stat_compare_means(
                method      = testmethod,
                label       = label_use,
                comparisons = my_comparisons_lab,
                hide.ns     = hide_ns_use,
                size        = size_use
              )
            }
            # 2군 비교(t.test/wilcox)는 전역 p 없이도 표시
          } else if (testmethod %in% c("wilcox.test", "t.test")) {
            p1 <- p1 + ggpubr::stat_compare_means(
              method      = testmethod,
              label       = label_use,
              comparisons = my_comparisons_lab,
              hide.ns     = hide_ns_use,
              size        = size_use
            )
          }
        }

        # Close an image (원본)


        if (!is.null(title)) {
          p1 <- p1 + ggtitle(sprintf("%s%s%s%s", title,
                                     ifelse(is.null(test.name), "", paste("\n",test.name, " ")),
                                     ifelse(is.null(pval), "", paste("p=", " ")),
                                     ifelse(is.null(pval), "", paste(pval, " "))))
        } else{
          p1 <- p1 + ggtitle(sprintf("%s%s%s%s", mvar,
                                     ifelse(is.null(test.name), "", paste("\n",test.name, " ")),
                                     ifelse(is.null(pval), "", paste("p=", " ")),
                                     ifelse(is.null(pval), "", paste(pval, " "))))

        }

        # y axis limit (원본)
        if(!is.null(ylim)){
          if(oc == "Chao1"){
            p1 = p1
          }else{
            p1 = p1 + ylim(ylim[1] , ylim[2])
          }
        }
        # facet (원본)
        if (length(facet) >= 1) {
          if(is.null(ncol)){
            ncol <- length(unique(df[,facet]))
          }
          p1 = p1 + facet_wrap(as.formula(sprintf("~ %s" , paste(setdiff(facet, "SocpleType"), collapse="+"))), scales="free_x", ncol = ncol)
          p1 = p1 + guides(color = "none", size = "none", shape= "none")
        } else {
          p1 = p1 + guides(color = "none", size = "none", shape= "none")
        }

        # plot size ratio (원본)
        if (length(unique(df.na[,mvar])) < 5){
          if(standardsize==TRUE){
            num.subgroup <- length(unique(df.na[,mvar]))*0.1
          }else{
            num.subgroup <- 0.9
          }
        }else{
          num.subgroup <- length(unique(df.na[,mvar]))*0.1
        }

        p1  <- p1  + theme(panel.grid = element_blank(),
                           panel.background = element_rect(fill = "white", colour = "Black",size = 0.5, linetype = "solid"),
                           aspect.ratio = 1/num.subgroup)

        plotlist[[length(plotlist)+1]] <- p1
      }
    }
  }
  multiplot(plotlist=plotlist, cols=plotCols, rows=plotRows)
  dev.off()
}
