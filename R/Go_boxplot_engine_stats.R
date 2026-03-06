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
                                    facet = NULL,
                                    p_adjust = "BH") {
  has_covariates <- !is.null(covariates) && length(covariates) > 0
  model <- if (is.null(model)) {
    if (!is.null(paired)) {
      "lmm"
    } else if (has_covariates) {
      "parametric"
    } else if (isTRUE(parametric)) {
      "parametric"
    } else {
      "nonparametric"
    }
  } else {
    tolower(as.character(model))
  }
  if (!model %in% c("nonparametric", "parametric", "lmm")) {
    stop("`model` must be one of: 'nonparametric', 'parametric', 'lmm'.")
  }

  use_lmm <- identical(model, "lmm")
  use_param <- identical(model, "parametric")
  use_nonparam <- identical(model, "nonparametric")

  facet_vars <- if (!is.null(facet)) {
    intersect(as.character(facet), colnames(df))
  } else {
    character(0)
  }

  cols <- unique(c(mvar, oc, covariates, paired, facet_vars))
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

  compute_single <- function(dat_sub) {
    cov_term <- if (length(covariates) > 0) paste(covariates, collapse = " + ") else NULL
    fixed_term <- paste(c(mvar, cov_term), collapse = " + ")
    form <- stats::as.formula(sprintf("%s ~ %s", oc, fixed_term))

    if (use_nonparam) {
      form_np <- stats::as.formula(sprintf("%s ~ %s", oc, mvar))
      if (nlevels(dat_sub[[mvar]]) > 2) {
        test <- stats::kruskal.test(form_np, dat_sub)
        return(list(test.name = "KW", pval = round(test$p.value, 4), testmethod = "wilcox.test", annotation = NULL))
      }
      return(list(test.name = "Pairwise Wilcoxon", pval = NULL, testmethod = "wilcox.test", annotation = NULL))
    }

    if (use_param && length(covariates) == 0) {
      if (nlevels(dat_sub[[mvar]]) > 2) {
        test <- stats::aov(form, dat_sub)
        pval <- summary(test)[[1]][["Pr(>F)"]][1]
        return(list(test.name = "ANOVA", pval = round(pval, 4), testmethod = "t.test", annotation = NULL))
      }
      return(list(test.name = "Pairwise T-Test", pval = NULL, testmethod = "t.test", annotation = NULL))
    }

    if (use_lmm) {
      if (is.null(paired) || !paired %in% names(dat_sub)) {
        return(list(test.name = NULL, pval = NULL, testmethod = NULL, annotation = NULL))
      }
      lmm_form <- stats::as.formula(sprintf("%s ~ %s + (1|%s)", oc, fixed_term, paired))
      fit <- try(lmerTest::lmer(lmm_form, data = dat_sub), silent = TRUE)
      if (inherits(fit, "try-error")) {
        return(list(test.name = NULL, pval = NULL, testmethod = NULL, annotation = NULL))
      }
      atab <- suppressMessages(stats::anova(fit))
      pval <- if ("Pr(>F)" %in% colnames(atab) && mvar %in% rownames(atab)) atab[mvar, "Pr(>F)"] else NA_real_
      emm <- try(emmeans::emmeans(fit, stats::as.formula(paste0("~", mvar))), silent = TRUE)
      if (inherits(emm, "try-error")) {
        return(list(test.name = "LMM", pval = ifelse(is.na(pval), NULL, round(pval, 4)), testmethod = NULL, annotation = NULL))
      }
      ctr <- try(emmeans::contrast(emm, method = "pairwise", adjust = "none"), silent = TRUE)
      if (inherits(ctr, "try-error")) {
        return(list(test.name = "LMM", pval = ifelse(is.na(pval), NULL, round(pval, 4)), testmethod = NULL, annotation = NULL))
      }
    } else {
      fit <- try(stats::lm(form, data = dat_sub), silent = TRUE)
      if (inherits(fit, "try-error")) {
        return(list(test.name = NULL, pval = NULL, testmethod = NULL, annotation = NULL))
      }
      atab <- stats::anova(fit)
      pval <- if ("Pr(>F)" %in% colnames(atab) && mvar %in% rownames(atab)) atab[mvar, "Pr(>F)"] else NA_real_
      emm <- try(emmeans::emmeans(fit, stats::as.formula(paste0("~", mvar))), silent = TRUE)
      if (inherits(emm, "try-error")) {
        return(list(test.name = "ANCOVA", pval = ifelse(is.na(pval), NULL, round(pval, 4)), testmethod = NULL, annotation = NULL))
      }
      ctr <- try(emmeans::contrast(emm, method = "pairwise", adjust = "none"), silent = TRUE)
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

    y_base <- max(dat_sub[[oc]], na.rm = TRUE)
    y_span <- diff(range(dat_sub[[oc]], na.rm = TRUE))
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

  if (length(facet_vars) == 0) {
    return(compute_single(dat))
  }

  split_ids <- split(seq_len(nrow(dat)), interaction(dat[, facet_vars, drop = FALSE], drop = TRUE, lex.order = TRUE))
  res_list <- list()
  for (idx in seq_along(split_ids)) {
    sub_idx <- split_ids[[idx]]
    dat_sub <- droplevels(dat[sub_idx, , drop = FALSE])
    if (nlevels(dat_sub[[mvar]]) < 2) next
    sub_res <- compute_single(dat_sub)
    if (!is.null(sub_res$annotation) && nrow(sub_res$annotation) > 0) {
      for (fc in facet_vars) {
        fc_levels <- levels(dat[[fc]])
        fc_value <- as.character(dat_sub[[fc]][1])
        if (is.null(fc_levels)) {
          sub_res$annotation[[fc]] <- fc_value
        } else {
          sub_res$annotation[[fc]] <- factor(fc_value, levels = fc_levels)
        }
      }
    }
    res_list[[length(res_list) + 1]] <- sub_res
  }

  if (length(res_list) == 0) {
    return(list(test.name = NULL, pval = NULL, testmethod = NULL, annotation = NULL))
  }

  test_name <- res_list[[1]]$test.name
  ann_list <- lapply(res_list, function(x) x$annotation)
  ann_list <- ann_list[vapply(ann_list, function(x) !is.null(x) && nrow(x) > 0, logical(1))]
  ann <- if (length(ann_list) > 0) do.call(rbind, ann_list) else NULL

  list(
    test.name = test_name,
    pval = NULL,
    testmethod = NULL,
    annotation = ann
  )
}

Go_boxplot_add_stats_layer <- function(p1, stat_res, my_comparisons,
                                       paired = NULL, cutoff = 0.1) {
  if (is.null(stat_res$test.name)) return(p1)

  if (!is.null(stat_res$annotation)) {
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

  label_type <- "p.format"
  size_val <- 2

  if (stat_res$test.name %in% c("KW", "ANOVA")) {
    if (!is.null(stat_res$pval) && stat_res$pval >= cutoff) return(p1)
    if (is.null(stat_res$testmethod)) return(p1)
    return(
      p1 + ggpubr::stat_compare_means(
        method = stat_res$testmethod,
        label = label_type,
        comparisons = my_comparisons,
        hide.ns = FALSE,
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
          hide.ns = FALSE,
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

# Lightweight regression checks for engine logic.
Go_boxplot_stats_engine_smoke_test <- function() {
  set.seed(1)
  dat <- data.frame(
    grp = rep(c("A", "B", "C"), each = 8),
    y = rnorm(24),
    z = sample(c("F", "M"), 24, replace = TRUE),
    id = rep(1:8, 3)
  )
  comps <- list(c("A", "B"), c("A", "C"), c("B", "C"))

  res_np <- Go_boxplot_stats_engine(
    df = dat, mvar = "grp", oc = "y", comparisons = comps,
    model = "nonparametric", covariates = "z"
  )
  stopifnot(res_np$test.name %in% c("KW", "Pairwise Wilcoxon"))

  res_cov <- Go_boxplot_stats_engine(
    df = dat, mvar = "grp", oc = "y", comparisons = comps,
    covariates = "z"
  )
  stopifnot(identical(res_cov$test.name, "ANCOVA"))

  invisible(TRUE)
}
