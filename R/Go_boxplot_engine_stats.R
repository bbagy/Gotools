clean_boxplot_group_label <- function(x) {
  sub(" \\(n=.*\\)$", "", as.character(x))
}

make_pair_key <- function(a, b) {
  paste(sort(c(as.character(a), as.character(b))), collapse = "||")
}

is_small_scale_boxplot <- function(y) {
  y_vals <- suppressWarnings(as.numeric(y))
  y_vals <- y_vals[is.finite(y_vals)]
  if (length(y_vals) == 0) {
    return(FALSE)
  }
  max(abs(y_vals), na.rm = TRUE) <= 0.05
}

compute_visible_boxplot_bound <- function(y, groups = NULL, which = c("upper", "lower")) {
  which <- match.arg(which)
  y_vals <- suppressWarnings(as.numeric(y))
  keep <- is.finite(y_vals)
  y_vals <- y_vals[keep]
  if (!length(y_vals)) {
    return(NA_real_)
  }

  if (is.null(groups)) {
    groups <- rep("all", length(y))
  }
  groups <- as.character(groups)[keep]
  split_vals <- split(y_vals, groups, drop = TRUE)
  idx <- if (which == "upper") 5 else 1
  bounds <- vapply(split_vals, function(vals) {
    vals <- vals[is.finite(vals)]
    if (!length(vals)) return(NA_real_)
    grDevices::boxplot.stats(vals)$stats[idx]
  }, numeric(1))

  out <- if (which == "upper") max(bounds, na.rm = TRUE) else min(bounds, na.rm = TRUE)
  if (!is.finite(out)) {
    out <- if (which == "upper") max(y_vals, na.rm = TRUE) else min(y_vals, na.rm = TRUE)
  }
  out
}

assign_bracket_layers <- function(comparisons, group_levels) {
  if (length(comparisons) == 0 || length(group_levels) == 0) {
    return(integer(0))
  }

  level_index <- stats::setNames(seq_along(group_levels), group_levels)
  bracket_df <- data.frame(
    idx = seq_along(comparisons),
    group1 = vapply(comparisons, `[`, character(1), 1),
    group2 = vapply(comparisons, `[`, character(1), 2),
    stringsAsFactors = FALSE
  )

  bracket_df$start <- pmin(level_index[bracket_df$group1], level_index[bracket_df$group2])
  bracket_df$end <- pmax(level_index[bracket_df$group1], level_index[bracket_df$group2])
  bracket_df <- bracket_df[stats::complete.cases(bracket_df$start, bracket_df$end), , drop = FALSE]
  if (!nrow(bracket_df)) {
    return(rep(NA_integer_, length(comparisons)))
  }

  bracket_df$span <- bracket_df$end - bracket_df$start
  bracket_df <- bracket_df[order(bracket_df$span, bracket_df$start, bracket_df$end), , drop = FALSE]

  layer_assignments <- rep(NA_integer_, length(comparisons))
  layer_spans <- list()

  overlaps_existing <- function(start, end, spans) {
    if (length(spans) == 0) {
      return(FALSE)
    }
    any(vapply(spans, function(x) !(end < x[1] || start > x[2]), logical(1)))
  }

  for (i in seq_len(nrow(bracket_df))) {
    start_i <- bracket_df$start[i]
    end_i <- bracket_df$end[i]
    layer_id <- 1L
    current_spans <- if (layer_id <= length(layer_spans)) layer_spans[[layer_id]] else list()
    while (overlaps_existing(start_i, end_i, current_spans)) {
      layer_id <- layer_id + 1L
      current_spans <- if (layer_id <= length(layer_spans)) layer_spans[[layer_id]] else list()
    }
    layer_assignments[bracket_df$idx[i]] <- layer_id
    existing_spans <- if (layer_id <= length(layer_spans)) layer_spans[[layer_id]] else list()
    layer_spans[[layer_id]] <- c(existing_spans, list(c(start_i, end_i)))
  }

  layer_assignments
}

compute_annotation_positions <- function(y, groups = NULL, n_labels, comparisons = NULL, group_levels = NULL) {
  y_vals <- suppressWarnings(as.numeric(y))
  y_vals <- y_vals[is.finite(y_vals)]
  if (length(y_vals) == 0 || n_labels <= 0) {
    return(NULL)
  }

  y_max <- compute_visible_boxplot_bound(y, groups, which = "upper")
  y_min <- compute_visible_boxplot_bound(y, groups, which = "lower")
  if (!is.finite(y_max)) y_max <- max(y_vals, na.rm = TRUE)
  if (!is.finite(y_min)) y_min <- min(y_vals, na.rm = TRUE)
  y_span <- y_max - y_min
  small_scale <- is_small_scale_boxplot(y_vals)

  if (!is.finite(y_span) || y_span <= 0) {
    y_span <- if (small_scale) max(abs(y_max) * 0.05, 1e-5) else max(abs(y_max) * 0.15, 1e-6)
  }

  if (small_scale) {
    y_step <- max(y_span * 0.018, abs(y_max) * 0.009, 5e-5)
    y_base <- y_max + y_step
  } else {
    y_step <- max(y_span * 0.04, abs(y_max) * 0.025, 1e-6)
    y_base <- y_max
  }

  if (is.null(comparisons) || is.null(group_levels)) {
    return(y_base + y_step * seq_len(n_labels))
  }

  layers <- assign_bracket_layers(comparisons = comparisons, group_levels = group_levels)
  if (length(layers) != n_labels || all(is.na(layers))) {
    return(y_base + y_step * seq_len(n_labels))
  }

  y_base + y_step * layers
}

Go_boxplot_stats_engine <- function(df, mvar, oc, comparisons,
                                    model = NULL,
                                    covariates = NULL,
                                    paired = NULL,
                                    facet = NULL,
                                    p_adjust = "BH") {
  build_local_comparisons <- function(group_levels) {
    group_levels <- as.character(group_levels)
    n_grp <- length(group_levels)
    if (n_grp < 2) {
      return(list())
    }
    if (n_grp >= 5) {
      ref_grp <- group_levels[1]
      return(lapply(group_levels[-1], function(g) c(ref_grp, g)))
    }
    lapply(seq_len(ncol(combn(group_levels, 2))), function(i) combn(group_levels, 2)[, i])
  }

  has_covariates <- !is.null(covariates) && length(covariates) > 0
  model <- if (is.null(model)) {
    if (!is.null(paired) && length(paired) > 0) {
      "lmm"
    } else if (has_covariates) {
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
    present_levels <- levels(droplevels(dat_sub[[mvar]]))
    comparisons_use <- comparisons
    if (length(facet_vars) > 0) {
      comparisons_use <- build_local_comparisons(present_levels)
    } else if (!is.null(comparisons) && length(comparisons) > 0) {
      comparisons_use <- Filter(function(comp) all(comp %in% present_levels), comparisons)
      if (length(comparisons_use) == 0) {
        comparisons_use <- build_local_comparisons(present_levels)
      }
    } else {
      comparisons_use <- build_local_comparisons(present_levels)
    }

    cov_term <- if (length(covariates) > 0) paste(covariates, collapse = " + ") else NULL
    fixed_term <- paste(c(mvar, cov_term), collapse = " + ")
    form <- stats::as.formula(sprintf("%s ~ %s", oc, fixed_term))

    if (use_nonparam) {
      form_np <- stats::as.formula(sprintf("%s ~ %s", oc, mvar))
      kw_pval <- NULL
      if (nlevels(dat_sub[[mvar]]) > 2) {
        test <- stats::kruskal.test(form_np, dat_sub)
        kw_pval <- round(test$p.value, 4)
      }

      raw_pvals <- vapply(comparisons_use, function(comp) {
        sub_dat <- dat_sub[dat_sub[[mvar]] %in% comp, , drop = FALSE]
        sub_dat[[mvar]] <- droplevels(factor(sub_dat[[mvar]]))
        if (nlevels(sub_dat[[mvar]]) < 2) return(NA_real_)
        wt <- try(stats::wilcox.test(
          stats::as.formula(sprintf("%s ~ %s", oc, mvar)),
          data = sub_dat,
          exact = FALSE
        ), silent = TRUE)
        if (inherits(wt, "try-error")) return(NA_real_)
        as.numeric(wt$p.value)
      }, numeric(1))

      adj_pvals <- stats::p.adjust(raw_pvals, method = p_adjust)
      ann_y <- compute_annotation_positions(
        dat_sub[[oc]],
        dat_sub[[mvar]],
        length(comparisons_use),
        comparisons = comparisons_use,
        group_levels = levels(dat_sub[[mvar]])
      )
      ann <- data.frame(
        group1 = vapply(comparisons_use, `[`, character(1), 1),
        group2 = vapply(comparisons_use, `[`, character(1), 2),
        p = raw_pvals,
        p.adj = adj_pvals,
        y.position = ann_y,
        stringsAsFactors = FALSE
      )
      ann <- ann[is.finite(ann$p), , drop = FALSE]
      if (nrow(ann) > 0) {
        ann$label <- as.character(signif(ann$p.adj, 3))
      }

      return(list(test.name = "KW", pval = kw_pval, testmethod = NULL, annotation = ann))
    }

    if (use_lmm) {
      if (is.null(paired) || !paired %in% names(dat_sub)) {
        return(list(test.name = NULL, pval = NULL, testmethod = NULL, annotation = NULL))
      }
      lmm_form <- stats::as.formula(sprintf("%s ~ %s + (1|%s)", oc, fixed_term, paired))
      fit <- try(lmerTest::lmer(lmm_form, data = dat_sub), silent = TRUE)
      if (inherits(fit, "try-error")) {
        message(sprintf("[Go_boxplot] LMM fit failed for %s ~ %s: %s",
                        oc, fixed_term, as.character(fit)))
        return(list(test.name = "LMM", pval = NULL, testmethod = NULL, annotation = NULL))
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

    display_keys <- vapply(comparisons_use, function(comp) make_pair_key(clean_boxplot_group_label(comp[1]), clean_boxplot_group_label(comp[2])), character(1))
    p_map <- setNames(ctr_df$p.value, ctr_df$key)

    pvals <- p_map[display_keys]
    if ("p.value" %in% colnames(ctr_df)) {
      pvals_adj <- stats::p.adjust(pvals, method = p_adjust)
    } else {
      pvals_adj <- pvals
    }

    ann_y <- compute_annotation_positions(
      dat_sub[[oc]],
      dat_sub[[mvar]],
      length(comparisons_use),
      comparisons = comparisons_use,
      group_levels = levels(dat_sub[[mvar]])
    )

    ann <- data.frame(
      group1 = vapply(comparisons_use, `[`, character(1), 1),
      group2 = vapply(comparisons_use, `[`, character(1), 2),
      p = as.numeric(pvals),
      p.adj = as.numeric(pvals_adj),
      y.position = ann_y,
      stringsAsFactors = FALSE
    )
    ann <- ann[is.finite(ann$p), , drop = FALSE]
    if (nrow(ann) > 0) {
      ann$label <- as.character(signif(ann$p.adj, 3))
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
    test.name  = test_name,
    pval       = NULL,
    testmethod = res_list[[1]]$testmethod,
    annotation = ann
  )
}

compute_boxplot_label_y <- function(y, groups = NULL, n_labels) compute_annotation_positions(y, groups, n_labels)

Go_boxplot_add_stats_layer <- function(p1, stat_res, my_comparisons,
                                       paired = NULL, cutoff = 0.1,
                                       dat = NULL, oc = NULL, mvar = NULL,
                                       label_y_override = NULL) {
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
  stopifnot(identical(res_np$test.name, "KW"))

  res_cov <- Go_boxplot_stats_engine(
    df = dat, mvar = "grp", oc = "y", comparisons = comps,
    covariates = "z"
  )
  stopifnot(identical(res_cov$test.name, "ANCOVA"))

  invisible(TRUE)
}
