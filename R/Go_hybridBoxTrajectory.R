#' Hybrid Boxplot + Trajectory Plot for Repeated Measures
#'
#' Builds a hybrid figure with a reference boxplot panel on the left and a
#' repeated-measures trajectory panel on the right for one or more outcomes.
#' The trajectory panel is summarized with a linear mixed model (LMM) and
#' optional `emmeans` contrasts.
#'
#' @param input Primary input object. Provide either a data frame or a
#'   phyloseq object here. This is the preferred interface.
#' @param project Project name used for output folder and file naming.
#' @param df Deprecated compatibility input for a data frame. If supplied,
#'   it will be merged into \code{input} internally.
#' @param psIN Deprecated compatibility input for a phyloseq object. If
#'   supplied, it will be merged into \code{input} internally.
#' @param outcomes Numeric outcomes to plot. When \code{psIN} is provided,
#'   these should be alpha-diversity metrics such as \code{"Shannon"} or
#'   \code{"Chao1"}.
#' @param ref.var Grouping variable used for the left reference boxplot.
#' @param traj.var Grouping variable used for the right trajectory panel.
#' @param orders Optional ordering specification for the grouping variables.
#'   You can provide either a single character vector, which will be filtered
#'   with \code{intersect()} against each panel variable, or a named list such
#'   as \code{list(ref = c("NC", "CD"), traj = c("Neg", "Mild", "Mod"))}.
#'   A single vector keeps the function aligned with the broader Gotool style,
#'   while the list form allows different order sets for the two panels.
#' @param time.var Numeric time variable for the trajectory x-axis.
#' @param subject.var Subject ID used for repeated-measures lines and the
#'   random intercept in the LMM.
#' @param name Optional name tag appended to output files.
#' @param mycols Optional color vector. You can provide either a named vector
#'   such as \code{c(NC = "#4C78A8", CD = "#E45756", Neg = "#72B7B2",
#'   Mild = "#F58518", Mod = "#B279A2")} or an unnamed palette vector such as
#'   the output of \code{Go_myCols()}, which will be assigned in order to the
#'   panel levels after \code{orders} are resolved.
#' @param covariates Optional covariate column names to include in the trajectory
#'   LMM as fixed effects. These are added additively after the main
#'   \code{time.var * traj.var} terms.
#' @param statistics Logical; if \code{TRUE}, fit an LMM for each outcome in
#'   the trajectory panel and export ANOVA / emmeans tables.
#' @param p_adjust Logical; if \code{TRUE}, use BH adjustment for pairwise
#'   `emmeans` contrasts and for the displayed within-model FDR label.
#' @param log_transform One of \code{"auto"}, \code{"none"}, or \code{"log1p"}.
#'   In \code{"auto"} mode, richness-like outcomes such as Chao1, ACE, and
#'   Observed use `log1p()` in the LMM.
#' @param emmeans_at Optional named list passed to \code{emmeans(..., at = ...)}.
#'   If \code{NULL} and \code{time.var} is numeric, group means are estimated at
#'   the median observed time.
#' @param span LOESS span for the colored smooth trajectory.
#' @param loess_min_n Minimum observations per group for LOESS overlay.
#' @param lm_min_n Minimum observations per group for dashed LM overlay.
#' @param width Output PDF width in inches.
#' @param height Output PDF height in inches. If \code{NULL}, height is set from
#'   the number of outcomes.
#'
#' @return Invisibly returns a list containing plots, fitted model summaries,
#'   and output paths.
#'
#' @examples
#' \dontrun{
#' Go_hybridBoxTrajectory(
#'   input = dat,
#'   project = "MAPS",
#'   outcomes = c("Chao1", "Shannon"),
#'   ref.var = "TruncGroup",
#'   traj.var = "TruncITxOnly",
#'   orders = list(ref = c("NC", "CD"), traj = c("Neg", "Mild", "Mod")),
#'   mycols = c(
#'     NC = "#4C78A8",
#'     CD = "#E45756",
#'     Neg = "#72B7B2",
#'     Mild = "#F58518",
#'     Mod = "#B279A2"
#'   ),
#'   time.var = "POD_num",
#'   subject.var = "subject_id"
#' )
#' }
#'
#' @param patchwork Logical. If \code{TRUE}, skip saving and return the plot object(s) for use with \code{Gg_patchwork()} or the \pkg{patchwork} package. Default \code{FALSE}.
#' @export
Go_hybridBoxTrajectory <- function(input = NULL,
                                   project = NULL,
                                   outcomes,
                                   ref.var,
                                   traj.var,
                                   orders = NULL,
                                   time.var,
                                   subject.var,
                                   name = NULL,
                                   mycols = NULL,
                                   covariates = NULL,
                                   statistics = TRUE,
                                   p_adjust = TRUE,
                                   log_transform = c("auto", "none", "log1p"),
                                   emmeans_at = NULL,
                                   span = 0.85,
                                   loess_min_n = 5,
                                   lm_min_n = 3,
                                   width = 8,
                                   height = NULL,
                                   df = NULL,
                                   psIN = NULL,
                                   patchwork = FALSE) {

  log_transform <- match.arg(log_transform)

  # Backward-compatible rescue for calls like
  # Go_hybridBoxTrajectory(df = dat, project_name, ...)
  if (is.null(project) &&
      is.character(input) &&
      length(input) == 1 &&
      !is.null(df) &&
      inherits(df, "data.frame")) {
    project <- input
    input <- df
    df <- NULL
  }

  if (is.null(project) &&
      is.character(input) &&
      length(input) == 1 &&
      !is.null(psIN) &&
      inherits(psIN, "phyloseq")) {
    project <- input
    input <- psIN
    psIN <- NULL
  }

  input_candidates <- list(input = input, df = df, psIN = psIN)
  provided_idx <- which(!vapply(input_candidates, is.null, logical(1)))

  if (length(provided_idx) == 0) {
    stop("Provide `input` as either a data.frame or a phyloseq object.")
  }

  if (length(provided_idx) > 1) {
    stop("Provide only one of `input`, `df`, or `psIN`.")
  }

  if (missing(project) || is.null(project) || !nzchar(project)) {
    stop("`project` is required.")
  }
  if (length(outcomes) == 0) {
    stop("`outcomes` must contain at least one variable.")
  }

  resolve_input_df <- function(input, outcomes, project, name) {
    if (inherits(input, "data.frame")) {
      return(input)
    }

    if (!inherits(input, "phyloseq")) {
      stop(
        "`input` must be either a data.frame (or tibble/data.frame-like object) ",
        "or a phyloseq object."
      )
    }

    alpha_metrics <- c("Observed", "Chao1", "ACE", "Shannon",
                       "Simpson", "InvSimpson", "Fisher", "PD")
    bad <- setdiff(outcomes, alpha_metrics)
    if (length(bad) > 0) {
      stop("`phyloseq` input currently supports alpha-diversity outcomes only. Unsupported outcomes: ",
           paste(bad, collapse = ", "))
    }

    Go_adiv(psIN = input, project = project, alpha_metrics = outcomes, name = name)
  }

  safe_numeric <- function(x) suppressWarnings(as.numeric(x))

  format_p_value <- function(p) {
    if (length(p) == 0 || is.na(p)) {
      return("NA")
    }
    if (p < 0.001) {
      return("<0.001")
    }
    formatC(p, digits = 3, format = "f")
  }

  format_covariate_label <- function(df, covariates) {
    if (is.null(covariates) || length(covariates) == 0) {
      return(NULL)
    }

    parts <- vapply(covariates, function(covar) {
      if (!covar %in% names(df)) {
        return(sprintf("%s (?)", covar))
      }
      x <- df[[covar]]
      cov_type <- if (is.factor(x) || is.character(x)) "fct" else "num"
      sprintf("%s (%s)", covar, cov_type)
    }, character(1))

    paste(parts, collapse = ", ")
  }

  is_log_metric <- function(metric, mode) {
    if (identical(mode, "log1p")) {
      return(TRUE)
    }
    if (identical(mode, "none")) {
      return(FALSE)
    }
    metric %in% c("Observed", "Chao1", "ACE")
  }

  normalize_mycols <- function(mycols) {
    if (is.null(mycols)) {
      return(NULL)
    }

    if (!is.atomic(mycols)) {
      stop("`mycols` must be an atomic color vector.")
    }

    mycols <- as.character(mycols)
    nm <- names(mycols)
    if (is.null(nm)) {
      return(mycols)
    }

    nm <- trimws(nm)
    if (all(!nzchar(nm))) {
      names(mycols) <- NULL
      return(mycols)
    }

    if (any(!nzchar(nm))) {
      stop(
        "`mycols` with names must name every color. ",
        "Use either a fully named vector or a fully unnamed palette vector."
      )
    }

    names(mycols) <- nm
    mycols
  }

  build_palette <- function(mycols, levels_needed) {
    levels_needed <- unique(as.character(levels_needed))
    if (is.null(mycols)) {
      cols <- grDevices::hcl.colors(length(levels_needed), palette = "Dark 3")
      names(cols) <- levels_needed
      return(cols)
    }

    if (is.null(names(mycols))) {
      cols <- rep(mycols, length.out = length(levels_needed))
      names(cols) <- levels_needed
      return(cols)
    }

    missing_cols <- setdiff(levels_needed, names(mycols))
    if (length(missing_cols) > 0) {
      extra <- grDevices::hcl.colors(length(missing_cols), palette = "Dark 3")
      names(extra) <- missing_cols
      mycols <- c(mycols, extra)
    }
    mycols[levels_needed]
  }

  resolve_orders <- function(df, ref.var, traj.var, orders) {
    ref_seen <- unique(as.character(stats::na.omit(df[[ref.var]])))
    traj_seen <- unique(as.character(stats::na.omit(df[[traj.var]])))

    if (is.null(orders)) {
      ref_levels <- ref_seen
      traj_levels <- traj_seen
    } else if (is.atomic(orders) && !is.list(orders)) {
      ord_vec <- as.character(orders)
      ref_levels <- intersect(ord_vec, ref_seen)
      traj_levels <- intersect(ord_vec, traj_seen)
    } else {
      if (!is.list(orders)) {
        stop("`orders` must be NULL, a character vector, or a named list like list(ref = c(...), traj = c(...)).")
      }
      ref_levels <- intersect(as.character(orders$ref), ref_seen)
      traj_levels <- intersect(as.character(orders$traj), traj_seen)
    }

    if (is.null(ref_levels) || length(ref_levels) < 2) {
      stop("Reference panel needs at least two valid levels after applying `orders`.")
    }
    if (is.null(traj_levels) || length(traj_levels) < 2) {
      stop("Trajectory panel needs at least two valid levels after applying `orders`.")
    }

    list(
      ref = as.character(ref_levels),
      traj = as.character(traj_levels)
    )
  }

  finite_range <- function(...) {
    vals <- c(...)
    vals <- safe_numeric(vals)
    vals <- vals[is.finite(vals)]
    if (!length(vals)) {
      return(c(0, 1))
    }
    rng <- range(vals, na.rm = TRUE)
    if (!all(is.finite(rng))) {
      return(c(0, 1))
    }
    if (diff(rng) == 0) {
      rng <- c(rng[1] - 0.5, rng[2] + 0.5)
    }
    rng
  }

  build_y_scale <- function(ref_vals, traj_vals, extra_vals = NULL) {
    y_range <- finite_range(ref_vals, traj_vals, extra_vals)
    y_pad <- diff(y_range) * 0.05
    if (!is.finite(y_pad) || y_pad == 0) {
      y_pad <- 0.1
    }
    y_limits <- c(y_range[1] - y_pad, y_range[2] + y_pad)
    y_breaks <- pretty(y_limits, n = 5)
    list(limits = y_limits, breaks = y_breaks)
  }

  build_emmeans_at <- function(use_df, time.var, emmeans_at) {
    if (!is.null(emmeans_at)) {
      return(emmeans_at)
    }

    time_vals <- safe_numeric(use_df[[time.var]])
    time_vals <- time_vals[is.finite(time_vals)]
    if (!length(time_vals)) {
      return(NULL)
    }

    out <- list(stats::median(time_vals, na.rm = TRUE))
    names(out) <- time.var
    out
  }

  sanitize_tag <- function(x) gsub("[^A-Za-z0-9._-]+", "-", as.character(x))

  fit_lmm_bundle <- function(df, outcome, time.var, group.var, subject.var,
                             covariates, log_transform, p_adjust, emmeans_at) {
    needed <- unique(c(outcome, time.var, group.var, subject.var, covariates))
    use_df <- df[, needed, drop = FALSE]
    rename_map <- c(
      outcome = outcome,
      time_value = time.var,
      group_value = group.var,
      subject_value = subject.var
    )
    for (nm in names(rename_map)) {
      names(use_df)[names(use_df) == rename_map[[nm]]] <- nm
    }
    use_df$outcome <- safe_numeric(use_df$outcome)
    use_df$time_value <- safe_numeric(use_df$time_value)
    use_df$group_value <- factor(as.character(use_df$group_value))
    use_df$subject_value <- as.character(use_df$subject_value)
    covariate_terms <- character(0)
    if (!is.null(covariates) && length(covariates) > 0) {
      for (covar in covariates) {
        if (!covar %in% names(use_df)) {
          next
        }
        if (is.character(use_df[[covar]]) || is.factor(use_df[[covar]])) {
          use_df[[covar]] <- factor(as.character(use_df[[covar]]))
        } else {
          use_df[[covar]] <- safe_numeric(use_df[[covar]])
        }
        covariate_terms <- c(covariate_terms, covar)
      }
    }
    covariate_label <- format_covariate_label(use_df, covariate_terms)
    use_df <- use_df[stats::complete.cases(use_df), , drop = FALSE]

    if (nrow(use_df) < 5 ||
        length(unique(use_df$subject_value)) < 2 ||
        length(unique(use_df$group_value)) < 2 ||
        length(unique(use_df$time_value)) < 2) {
      return(list(
        label = "lmer p / FDR\ninsufficient data",
        fit = NULL,
        anova = NULL,
        emmeans = NULL,
        pairs = NULL
      ))
    }

    rhs_terms <- c("time_value * group_value", covariate_terms, "(1 | subject_value)")
    lhs_term <- if (is_log_metric(outcome, log_transform)) "log1p(outcome)" else "outcome"
    form_txt <- paste(lhs_term, "~", paste(rhs_terms, collapse = " + "))

    fit <- tryCatch(
      lmerTest::lmer(stats::as.formula(form_txt), data = use_df, REML = FALSE),
      error = function(e) NULL
    )
    if (is.null(fit)) {
      return(list(
        label = "lmer p / FDR\nmodel failed",
        fit = NULL,
        anova = NULL,
        emmeans = NULL,
        pairs = NULL
      ))
    }

    anova_tbl <- tryCatch(as.data.frame(stats::anova(fit)), error = function(e) NULL)
    if (is.null(anova_tbl)) {
      return(list(
        label = "lmer p / FDR\nmodel failed",
        fit = fit,
        anova = NULL,
        emmeans = NULL,
        pairs = NULL
      ))
    }

    anova_tbl$term <- rownames(anova_tbl)
    rownames(anova_tbl) <- NULL
    p_col <- grep("^Pr\\(>F\\)$", names(anova_tbl), value = TRUE)
    anova_tbl$p_value <- if (length(p_col) == 1) anova_tbl[[p_col]] else NA_real_

    get_term_p <- function(pattern) {
      hit <- anova_tbl$term[grepl(pattern, anova_tbl$term)]
      if (length(hit) == 0) {
        return(NA_real_)
      }
      anova_tbl$p_value[match(hit[1], anova_tbl$term)]
    }

    raw_p <- c(
      time = get_term_p("^time_value$"),
      group = get_term_p("^group_value$"),
      `time x group` = get_term_p("time_value:group_value|group_value:time_value")
    )
    fdr_p <- stats::p.adjust(raw_p, method = if (isTRUE(p_adjust)) "BH" else "none")

    label_txt <- paste0(
      "lmer p / FDR\n",
      "time ", format_p_value(raw_p["time"]), " / ", format_p_value(fdr_p["time"]), "\n",
      "group ", format_p_value(raw_p["group"]), " / ", format_p_value(fdr_p["group"]), "\n",
      "time x group ", format_p_value(raw_p["time x group"]), " / ", format_p_value(fdr_p["time x group"]),
      if (!is.null(covariate_label)) paste0("\ncovariates ", covariate_label) else ""
    )

    emm_at <- build_emmeans_at(use_df, "time_value", emmeans_at)
    emm <- tryCatch(
      emmeans::emmeans(fit, specs = stats::as.formula("~ group_value"), at = emm_at),
      error = function(e) NULL
    )

    emm_tbl <- NULL
    pair_tbl <- NULL
    if (!is.null(emm)) {
      emm_tbl <- as.data.frame(emm)
      names(emm_tbl)[names(emm_tbl) == "group_value"] <- "group"

      pair_tbl <- tryCatch(
        as.data.frame(emmeans::pairs(emm, adjust = if (isTRUE(p_adjust)) "BH" else "none")),
        error = function(e) NULL
      )
    }

    anova_out <- dplyr::transmute(
      anova_tbl,
      outcome = outcome,
      model = "lmm",
      covariates = if (length(covariate_terms) > 0) paste(covariate_terms, collapse = ";") else NA_character_,
      term = term,
      num_df = if ("NumDF" %in% names(anova_tbl)) NumDF else NA_real_,
      den_df = if ("DenDF" %in% names(anova_tbl)) DenDF else NA_real_,
      statistic = if ("F value" %in% names(anova_tbl)) `F value` else NA_real_,
      p_value = p_value,
      fdr_within_model = stats::p.adjust(p_value, method = if (isTRUE(p_adjust)) "BH" else "none")
    )

    if (!is.null(emm_tbl)) {
      emm_tbl$outcome <- outcome
      emm_tbl$model <- "lmm"
      emm_tbl$covariates <- if (length(covariate_terms) > 0) paste(covariate_terms, collapse = ";") else NA_character_
      if (!is.null(emm_at) && length(emm_at) == 1) {
        emm_tbl$time_at <- unname(emm_at[[1]])
      }
    }

    if (!is.null(pair_tbl)) {
      pair_tbl$outcome <- outcome
      pair_tbl$model <- "lmm"
      pair_tbl$covariates <- if (length(covariate_terms) > 0) paste(covariate_terms, collapse = ";") else NA_character_
      if (!is.null(emm_at) && length(emm_at) == 1) {
        pair_tbl$time_at <- unname(emm_at[[1]])
      }
    }

    list(
      label = label_txt,
      fit = fit,
      anova = anova_out,
      emmeans = emm_tbl,
      pairs = pair_tbl
    )
  }

  build_default_ref_comparisons <- function(group_levels) {
    group_levels <- as.character(group_levels)
    n_grp <- length(group_levels)
    if (n_grp < 2) {
      return(list())
    }
    if (n_grp >= 5) {
      ref_grp <- group_levels[1]
      return(lapply(group_levels[-1], function(g) c(ref_grp, g)))
    }
    cmb <- combn(group_levels, 2)
    lapply(seq_len(ncol(cmb)), function(i) cmb[, i])
  }

  compute_simple_annotation_positions <- function(y, groups, comparisons) {
    y_vals <- safe_numeric(y)
    y_vals <- y_vals[is.finite(y_vals)]
    if (!length(y_vals) || length(comparisons) == 0) {
      return(NULL)
    }

    y_max <- max(y_vals, na.rm = TRUE)
    y_min <- min(y_vals, na.rm = TRUE)
    y_span <- y_max - y_min
    if (!is.finite(y_span) || y_span <= 0) {
      y_span <- max(abs(y_max) * 0.15, 1e-6)
    }
    step <- max(y_span * 0.08, abs(y_max) * 0.04, 1e-6)
    y_max + step * seq_along(comparisons)
  }

  compute_ref_stats <- function(ref_df, outcome, p_adjust) {
    levels_now <- levels(ref_df$display_group)
    comparisons <- build_default_ref_comparisons(levels_now)
    if (length(levels_now) < 2 || length(comparisons) == 0) {
      return(list(test.name = NULL, pval = NULL, testmethod = NULL, annotation = NULL))
    }

    raw_pvals <- vapply(comparisons, function(comp) {
      sub_df <- ref_df[ref_df$display_group %in% comp, , drop = FALSE]
      sub_df$display_group <- droplevels(sub_df$display_group)
      if (nlevels(sub_df$display_group) < 2) {
        return(NA_real_)
      }
      wt <- try(stats::wilcox.test(stats::as.formula(sprintf("%s ~ display_group", outcome)),
                                   data = sub_df, exact = FALSE),
                silent = TRUE)
      if (inherits(wt, "try-error")) {
        return(NA_real_)
      }
      as.numeric(wt$p.value)
    }, numeric(1))

    adj_pvals <- stats::p.adjust(raw_pvals, method = if (isTRUE(p_adjust)) "BH" else "none")
    ann <- data.frame(
      group1 = vapply(comparisons, `[`, character(1), 1),
      group2 = vapply(comparisons, `[`, character(1), 2),
      p = raw_pvals,
      p.adj = adj_pvals,
      y.position = compute_simple_annotation_positions(ref_df[[outcome]], ref_df$display_group, comparisons),
      stringsAsFactors = FALSE
    )
    ann <- ann[is.finite(ann$p), , drop = FALSE]
    if (nrow(ann) > 0) {
      ann$label <- as.character(signif(ann$p.adj, 3))
    }

    list(
      test.name = "Pairwise Wilcoxon",
      pval = NULL,
      testmethod = NULL,
      annotation = ann
    )
  }

  add_ref_stats_layer <- function(p1, stat_res, ref_df, outcome) {
    if (is.null(stat_res$test.name)) {
      return(p1)
    }

    if (!is.null(stat_res$annotation) && nrow(stat_res$annotation) > 0) {
      if (exists("Go_boxplot_add_stats_layer", mode = "function")) {
        return(Go_boxplot_add_stats_layer(
          p1 = p1,
          stat_res = stat_res,
          my_comparisons = Map(c, stat_res$annotation$group1, stat_res$annotation$group2),
          paired = NULL,
          cutoff = 0.1,
          dat = ref_df,
          oc = outcome,
          mvar = "display_group"
        ))
      }

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

  build_reference_plot <- function(ref_df, outcome, y_limits, y_breaks, ref_palette) {
    ggplot2::ggplot(ref_df, ggplot2::aes(x = display_group, y = .data[[outcome]],
                                         fill = display_group, color = display_group)) +
      ggplot2::geom_boxplot(outlier.shape = NA, alpha = 0.22, width = 0.62) +
      ggplot2::geom_jitter(width = 0.10, size = 2.1, alpha = 0.75) +
      ggplot2::scale_fill_manual(values = ref_palette, drop = FALSE) +
      ggplot2::scale_color_manual(values = ref_palette, drop = FALSE) +
      ggplot2::scale_y_continuous(limits = y_limits, breaks = y_breaks, expand = c(0, 0)) +
      ggplot2::labs(x = NULL, y = outcome,
                    title = paste0(outcome, " in ", ref.var, " reference samples")) +
      ggplot2::theme_bw(base_size = 10) +
      ggplot2::theme(
        legend.position = "none",
        plot.title = ggplot2::element_text(face = "bold", size = 10),
        plot.margin = ggplot2::margin(8, 2, 8, 8)
      )
  }

  build_trajectory_plot <- function(traj_df, outcome, y_limits, y_breaks, traj_palette,
                                    title_txt, label_txt) {
    line_df <- dplyr::arrange(traj_df, .data$subject_id, .data$time_value)

    loess_df <- dplyr::ungroup(
      dplyr::filter(
        dplyr::group_by(line_df, display_group),
        dplyr::n() >= loess_min_n,
        dplyr::n_distinct(time_value) >= 4
      )
    )

    lm_df <- dplyr::ungroup(
      dplyr::filter(
        dplyr::group_by(line_df, display_group),
        dplyr::n() >= lm_min_n,
        dplyr::n_distinct(time_value) >= 2
      )
    )

    ggplot2::ggplot(line_df, ggplot2::aes(x = time_value, y = .data[[outcome]], group = subject_id)) +
      ggplot2::geom_line(color = "grey70", linewidth = 0.45, alpha = 0.95) +
      ggplot2::geom_smooth(
        data = loess_df,
        ggplot2::aes(color = display_group, group = display_group),
        method = "loess",
        se = FALSE,
        linewidth = 1.1,
        span = span
      ) +
      ggplot2::geom_smooth(
        data = lm_df,
        ggplot2::aes(color = display_group, group = display_group),
        method = "lm",
        se = FALSE,
        linewidth = 0.8,
        linetype = "dashed"
      ) +
      ggplot2::geom_point(ggplot2::aes(color = display_group), size = 2.2, alpha = 0.25) +
      ggplot2::annotate("text", x = Inf, y = Inf, label = label_txt,
                        hjust = 1.02, vjust = 1.1, size = 2.8) +
      ggplot2::scale_color_manual(values = traj_palette, drop = FALSE) +
      ggplot2::scale_y_continuous(limits = y_limits, breaks = y_breaks, expand = c(0, 0)) +
      ggplot2::labs(x = NULL, y = NULL, title = title_txt) +
      ggplot2::theme_bw(base_size = 10) +
      ggplot2::theme(
        legend.position = "none",
        plot.title = ggplot2::element_text(face = "bold", size = 10),
        axis.title.y = ggplot2::element_blank(),
        axis.text.y = ggplot2::element_blank(),
        axis.ticks.y = ggplot2::element_blank(),
        axis.line.y = ggplot2::element_blank(),
        plot.margin = ggplot2::margin(8, -4, 8, 0)
      )
  }

  draw_side_annotation <- function(metric, traj_levels, traj_palette) {
    label_gp <- grid::gpar(fontsize = 9, col = "#222222")

    grid::grid.text(paste0(metric, " guide"), x = 0.02, y = 0.96,
                    just = c("left", "top"),
                    gp = grid::gpar(fontsize = 10, fontface = "bold"))

    grid::grid.lines(x = grid::unit(c(0.06, 0.22), "npc"),
                     y = grid::unit(c(0.84, 0.84), "npc"),
                     gp = grid::gpar(col = "#7A7A7A", lwd = 2))
    grid::grid.text("Patient line", x = 0.27, y = 0.84,
                    just = c("left", "center"), gp = label_gp)

    grid::grid.lines(x = grid::unit(c(0.06, 0.22), "npc"),
                     y = grid::unit(c(0.76, 0.76), "npc"),
                     gp = grid::gpar(col = "#7A7A7A", lwd = 2))
    grid::grid.text("LOESS", x = 0.27, y = 0.76,
                    just = c("left", "center"), gp = label_gp)

    grid::grid.lines(x = grid::unit(c(0.06, 0.22), "npc"),
                     y = grid::unit(c(0.68, 0.68), "npc"),
                     gp = grid::gpar(col = "#7A7A7A", lwd = 1.5, lty = 2))
    grid::grid.text("LM", x = 0.27, y = 0.68,
                    just = c("left", "center"), gp = label_gp)

    ys <- seq(0.52, max(0.32, 0.52 - 0.08 * (length(traj_levels) - 1)),
              length.out = length(traj_levels))
    for (i in seq_along(traj_levels)) {
      this_level <- traj_levels[i]
      this_col <- traj_palette[[this_level]]
      grid::grid.points(x = grid::unit(0.11, "npc"), y = grid::unit(ys[i], "npc"),
                        pch = 16, size = grid::unit(3.1, "mm"),
                        gp = grid::gpar(col = this_col, fill = this_col))
      grid::grid.text(this_level, x = 0.20, y = ys[i],
                      just = c("left", "center"), gp = label_gp)
    }
  }

  active_input <- input_candidates[[provided_idx]]
  df <- resolve_input_df(input = active_input, outcomes = outcomes, project = project, name = name)
  mycols <- normalize_mycols(mycols)

  required_cols <- unique(c(outcomes, ref.var, traj.var, time.var, subject.var))
  if (!is.null(covariates) && length(covariates) > 0) {
    required_cols <- unique(c(required_cols, covariates))
  }
  missing_cols <- setdiff(required_cols, names(df))
  if (length(missing_cols) > 0) {
    stop("Missing required columns in `df`: ", paste(missing_cols, collapse = ", "))
  }

  orders <- resolve_orders(df = df, ref.var = ref.var, traj.var = traj.var, orders = orders)
  ref.levels <- orders$ref
  traj.levels <- orders$traj

  out_dirs <- Go_path(project = project, pdf = "yes", table = "yes", path = NULL)
  out_pdf <- out_dirs$pdf
  out_tab <- out_dirs$tab
  stat_dir <- file.path(out_tab, "hybridBoxTrajectory")
  if (!dir.exists(stat_dir)) {
    dir.create(stat_dir, recursive = TRUE)
  }

  all_levels <- unique(c(ref.levels, traj.levels))
  all_palette <- build_palette(mycols, all_levels)
  ref_palette <- all_palette[ref.levels]
  traj_palette <- all_palette[traj.levels]

  ref_df <- df[df[[ref.var]] %in% ref.levels, , drop = FALSE]
  ref_df$display_group <- factor(ref_df[[ref.var]], levels = ref.levels)

  traj_df <- df[df[[traj.var]] %in% traj.levels, , drop = FALSE]
  traj_df$display_group <- factor(traj_df[[traj.var]], levels = traj.levels)
  traj_df$time_value <- safe_numeric(traj_df[[time.var]])
  traj_df$subject_id <- as.character(traj_df[[subject.var]])
  traj_df <- traj_df[!is.na(traj_df$time_value) &
                       !is.na(traj_df$subject_id) &
                       nzchar(traj_df$subject_id), , drop = FALSE]

  if (nrow(ref_df) == 0) {
    stop("No rows available for the reference panel after filtering `ref.var` / `ref.levels`.")
  }
  if (nrow(traj_df) == 0) {
    stop("No rows available for the trajectory panel after filtering `traj.var` / `traj.levels`.")
  }

  plot_store <- list()
  model_store <- list()

  for (outcome in outcomes) {
    if (!outcome %in% names(df)) {
      next
    }

    ref_metric_df <- ref_df
    ref_metric_df[[outcome]] <- safe_numeric(ref_metric_df[[outcome]])
    ref_metric_df <- ref_metric_df[!is.na(ref_metric_df[[outcome]]), , drop = FALSE]

    traj_metric_df <- traj_df
    traj_metric_df[[outcome]] <- safe_numeric(traj_metric_df[[outcome]])
    traj_metric_df <- traj_metric_df[!is.na(traj_metric_df[[outcome]]), , drop = FALSE]

    ref_stat_res <- if (isTRUE(statistics)) {
      compute_ref_stats(ref_metric_df, outcome, p_adjust = p_adjust)
    } else {
      list(test.name = NULL, pval = NULL, testmethod = NULL, annotation = NULL)
    }

    ref_ann_y <- NULL
    if (!is.null(ref_stat_res$annotation) && nrow(ref_stat_res$annotation) > 0) {
      ref_ann_y <- ref_stat_res$annotation$y.position
    }

    y_scale <- build_y_scale(ref_metric_df[[outcome]], traj_metric_df[[outcome]], extra_vals = ref_ann_y)

    model_res <- if (isTRUE(statistics)) {
      fit_lmm_bundle(
        df = traj_metric_df,
        outcome = outcome,
        time.var = "time_value",
        group.var = "display_group",
        subject.var = "subject_id",
        covariates = covariates,
        log_transform = log_transform,
        p_adjust = p_adjust,
        emmeans_at = emmeans_at
      )
    } else {
      list(label = "statistics disabled", fit = NULL, anova = NULL, emmeans = NULL, pairs = NULL)
    }

    p_box <- build_reference_plot(
      ref_df = ref_metric_df,
      outcome = outcome,
      y_limits = y_scale$limits,
      y_breaks = y_scale$breaks,
      ref_palette = ref_palette
    )
    p_box <- add_ref_stats_layer(
      p1 = p_box,
      stat_res = ref_stat_res,
      ref_df = ref_metric_df,
      outcome = outcome
    )

    p_traj <- build_trajectory_plot(
      traj_df = traj_metric_df,
      outcome = outcome,
      y_limits = y_scale$limits,
      y_breaks = y_scale$breaks,
      traj_palette = traj_palette,
      title_txt = paste0(outcome, " trajectories across ", time.var),
      label_txt = model_res$label
    )

    plot_store[[outcome]] <- list(box = p_box, traj = p_traj)
    model_store[[outcome]] <- model_res
  }

  if (length(plot_store) == 0) {
    message("[Go_hybridBoxTrajectory] No plots to render.")
    return(invisible(NULL))
  }

  file_stub <- paste0(
    "hybridBoxTrajectory.", project, ".",
    "ref-", sanitize_tag(ref.var), ".",
    "traj-", sanitize_tag(traj.var), ".",
    ifelse(is.null(name), "", paste0(name, ".")),
    format(Sys.Date(), "%y%m%d")
  )
  pdf_path <- file.path(out_pdf, paste0(file_stub, ".pdf"))
  if (is.null(height)) {
    height <- max(3, 3 * length(plot_store))
  }

  if (isTRUE(patchwork)) return(invisible(plot_store))
  grDevices::pdf(pdf_path, width = width, height = height)
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(layout = grid::grid.layout(
    nrow = length(plot_store),
    ncol = 3,
    widths = grid::unit(c(0.86, 1.52, 0.62), "null"),
    heights = grid::unit(rep(1, length(plot_store)), "null")
  )))

  row_idx <- 1
  for (outcome in names(plot_store)) {
    print(plot_store[[outcome]]$box,
          vp = grid::viewport(layout.pos.row = row_idx, layout.pos.col = 1))
    print(plot_store[[outcome]]$traj,
          vp = grid::viewport(layout.pos.row = row_idx, layout.pos.col = 2))
    grid::pushViewport(grid::viewport(layout.pos.row = row_idx, layout.pos.col = 3))
    draw_side_annotation(outcome, traj_levels = traj.levels, traj_palette = traj_palette)
    grid::upViewport()
    row_idx <- row_idx + 1
  }
  grDevices::dev.off()

  anova_tbl <- dplyr::bind_rows(lapply(model_store, `[[`, "anova"))
  emmeans_tbl <- dplyr::bind_rows(lapply(model_store, `[[`, "emmeans"))
  pairwise_tbl <- dplyr::bind_rows(lapply(model_store, `[[`, "pairs"))

  if (isTRUE(statistics)) {
    utils::write.csv(anova_tbl,
                     file = file.path(stat_dir, paste0(file_stub, ".lmm_anova.csv")),
                     row.names = FALSE)
    utils::write.csv(emmeans_tbl,
                     file = file.path(stat_dir, paste0(file_stub, ".emmeans.csv")),
                     row.names = FALSE)
    utils::write.csv(pairwise_tbl,
                     file = file.path(stat_dir, paste0(file_stub, ".emmeans_pairs.csv")),
                     row.names = FALSE)
  }

  invisible(list(
    plots = plot_store,
    models = model_store,
    outputs = list(
      pdf = pdf_path,
      anova = file.path(stat_dir, paste0(file_stub, ".lmm_anova.csv")),
      emmeans = file.path(stat_dir, paste0(file_stub, ".emmeans.csv")),
      pairs = file.path(stat_dir, paste0(file_stub, ".emmeans_pairs.csv"))
    )
  ))
}
