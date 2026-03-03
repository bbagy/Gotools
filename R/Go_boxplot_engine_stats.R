clean_boxplot_group_label <- function(x) {
  sub(" \\(n=.*\\)$", "", as.character(x))
}

make_pair_key <- function(a, b) {
  paste(sort(c(as.character(a), as.character(b))), collapse = "||")
}

Go_boxplot_stats_engine <- function(df, mvar, oc, comparisons,
                                    model = NULL,
                                    parametric = FALSE,
                                    covariates = NULL,
                                    paired = NULL,
                                    p_adjust = "BH") {
  model <- if (is.null(model)) {
    if (isTRUE(parametric)) "parametric" else "nonparametric"
  } else {
    tolower(as.character(model))
  }

  use_lmm <- identical(model, "lmm")
  use_param <- identical(model, "parametric")
  use_nonparam <- identical(model, "nonparametric")

  cols <- unique(c(mvar, oc, covariates, paired))
  dat <- df[, cols, drop = FALSE]
  dat <- dat[stats::complete.cases(dat), , drop = FALSE]
  if (!nrow(dat)) {
    return(list(test.name = NULL, pval = NULL, testmethod = NULL, annotation = NULL))
  }

  dat[[mvar]] <- factor(dat[[mvar]])
  levels_display <- levels(dat[[mvar]])
  levels_raw <- clean_boxplot_group_label(levels_display)
  levels(dat[[mvar]]) <- levels_raw
  dat[[mvar]] <- droplevels(dat[[mvar]])

  if (nlevels(dat[[mvar]]) < 2) {
    return(list(test.name = NULL, pval = NULL, testmethod = NULL, annotation = NULL))
  }

  cov_term <- if (length(covariates) > 0) paste(covariates, collapse = " + ") else NULL
  fixed_term <- paste(c(mvar, cov_term), collapse = " + ")
  form <- stats::as.formula(sprintf("%s ~ %s", oc, fixed_term))

  if (use_nonparam) {
    if (nlevels(dat[[mvar]]) > 2) {
      test <- stats::kruskal.test(form, dat)
      return(list(test.name = "KW", pval = round(test$p.value, 4), testmethod = "wilcox.test", annotation = NULL))
    }
    return(list(test.name = "Pairwise Wilcoxon", pval = NULL, testmethod = "wilcox.test", annotation = NULL))
  }

  if (use_param && length(covariates) == 0) {
    if (nlevels(dat[[mvar]]) > 2) {
      test <- stats::aov(form, dat)
      pval <- summary(test)[[1]][["Pr(>F)"]][1]
      return(list(test.name = "ANOVA", pval = round(pval, 4), testmethod = "t.test", annotation = NULL))
    }
    return(list(test.name = "Pairwise T-Test", pval = NULL, testmethod = "t.test", annotation = NULL))
  }

  if (use_lmm) {
    if (is.null(paired) || !paired %in% names(dat)) {
      return(list(test.name = NULL, pval = NULL, testmethod = NULL, annotation = NULL))
    }
    lmm_form <- stats::as.formula(sprintf("%s ~ %s + (1|%s)", oc, fixed_term, paired))
    fit <- try(lmerTest::lmer(lmm_form, data = dat), silent = TRUE)
    if (inherits(fit, "try-error")) {
      return(list(test.name = NULL, pval = NULL, testmethod = NULL, annotation = NULL))
    }
    atab <- suppressMessages(stats::anova(fit))
    pval <- if ("Pr(>F)" %in% colnames(atab) && mvar %in% rownames(atab)) atab[mvar, "Pr(>F)"] else NA_real_
    emm <- try(emmeans::emmeans(fit, stats::as.formula(paste0("~", mvar))), silent = TRUE)
    if (inherits(emm, "try-error")) {
      return(list(test.name = "LMM", pval = ifelse(is.na(pval), NULL, round(pval, 4)), testmethod = NULL, annotation = NULL))
    }
    ctr <- try(emmeans::contrast(emm, method = "pairwise", adjust = p_adjust), silent = TRUE)
    if (inherits(ctr, "try-error")) {
      return(list(test.name = "LMM", pval = ifelse(is.na(pval), NULL, round(pval, 4)), testmethod = NULL, annotation = NULL))
    }
  } else {
    fit <- try(stats::lm(form, data = dat), silent = TRUE)
    if (inherits(fit, "try-error")) {
      return(list(test.name = NULL, pval = NULL, testmethod = NULL, annotation = NULL))
    }
    atab <- stats::anova(fit)
    pval <- if ("Pr(>F)" %in% colnames(atab) && mvar %in% rownames(atab)) atab[mvar, "Pr(>F)"] else NA_real_
    emm <- try(emmeans::emmeans(fit, stats::as.formula(paste0("~", mvar))), silent = TRUE)
    if (inherits(emm, "try-error")) {
      return(list(test.name = "ANCOVA", pval = ifelse(is.na(pval), NULL, round(pval, 4)), testmethod = NULL, annotation = NULL))
    }
    ctr <- try(emmeans::contrast(emm, method = "pairwise", adjust = p_adjust), silent = TRUE)
    if (inherits(ctr, "try-error")) {
      return(list(test.name = "ANCOVA", pval = ifelse(is.na(pval), NULL, round(pval, 4)), testmethod = NULL, annotation = NULL))
    }
  }

  ctr_df <- as.data.frame(summary(ctr))
  parts <- strsplit(as.character(ctr_df$contrast), " - ", fixed = TRUE)
  ctr_df$group1_raw <- vapply(parts, `[`, character(1), 1)
  ctr_df$group2_raw <- vapply(parts, `[`, character(1), 2)
  ctr_df$key <- mapply(make_pair_key, ctr_df$group1_raw, ctr_df$group2_raw)

  display_keys <- vapply(comparisons, function(comp) make_pair_key(clean_boxplot_group_label(comp[1]), clean_boxplot_group_label(comp[2])), character(1))
  p_map <- setNames(ctr_df$p.value, ctr_df$key)

  pvals <- p_map[display_keys]
  if ("p.value" %in% colnames(ctr_df)) {
    pvals_adj <- stats::p.adjust(pvals, method = p_adjust)
  } else {
    pvals_adj <- pvals
  }

  y_base <- max(df[[oc]], na.rm = TRUE)
  y_span <- diff(range(df[[oc]], na.rm = TRUE))
  if (!is.finite(y_span) || y_span == 0) y_span <- 1
  y_step <- y_span * 0.08

  ann <- data.frame(
    group1 = vapply(comparisons, `[`, character(1), 1),
    group2 = vapply(comparisons, `[`, character(1), 2),
    p = as.numeric(pvals),
    p.adj = as.numeric(pvals_adj),
    y.position = y_base + y_step * seq_along(comparisons),
    stringsAsFactors = FALSE
  )
  ann <- ann[is.finite(ann$p), , drop = FALSE]
  if (nrow(ann) > 0) {
    ann$label <- paste0("p=", signif(ann$p.adj, 3))
  }

  list(
    test.name = if (use_lmm) "LMM" else "ANCOVA",
    pval = ifelse(is.na(pval), NULL, round(as.numeric(pval), 4)),
    testmethod = NULL,
    annotation = ann
  )
}

Go_boxplot_add_stats_layer <- function(p1, stat_res, my_comparisons,
                                       paired = NULL, cutoff = 0.1, star = TRUE) {
  if (is.null(stat_res$test.name)) return(p1)

  if (!is.null(stat_res$annotation)) {
    if (!is.null(stat_res$pval) && stat_res$pval >= cutoff) return(p1)
    if (nrow(stat_res$annotation) == 0) return(p1)
    return(
      p1 + ggpubr::stat_pvalue_manual(
        stat_res$annotation,
        label = "label",
        xmin = "group1",
        xmax = "group2",
        y.position = "y.position",
        inherit.aes = FALSE,
        size = 2
      )
    )
  }

  label_type <- if (isTRUE(star)) "p.signif" else "p.format"
  size_val <- if (isTRUE(star)) 3 else 2

  if (stat_res$test.name %in% c("KW", "ANOVA")) {
    if (!is.null(stat_res$pval) && stat_res$pval >= cutoff) return(p1)
    if (is.null(stat_res$testmethod)) return(p1)
    return(
      p1 + ggpubr::stat_compare_means(
        method = stat_res$testmethod,
        label = label_type,
        comparisons = my_comparisons,
        hide.ns = isTRUE(star),
        size = size_val
      )
    )
  }

  if (stat_res$testmethod %in% c("wilcox.test", "t.test")) {
    if (is.null(paired)) {
      return(
        p1 + ggpubr::stat_compare_means(
          method = stat_res$testmethod,
          label = label_type,
          comparisons = my_comparisons,
          hide.ns = isTRUE(star),
          size = size_val
        )
      )
    }
    return(
      p1 + ggpubr::stat_compare_means(
        method = stat_res$testmethod,
        label = "p.format",
        comparisons = my_comparisons,
        size = 2,
        paired = TRUE
      )
    )
  }

  p1
}
