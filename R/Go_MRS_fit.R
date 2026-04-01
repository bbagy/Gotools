#' Go_MRS_fit
#'
#' Fit an MRS-style model directly from a phyloseq object.
#'
#' @description
#' Uses \code{psIN} as the primary input. By default the microbial feature set
#' is ASV-level. Users can optionally aggregate to \code{"Species"},
#' \code{"Genus"}, or any taxonomic rank present in \code{tax_table(psIN)}.
#'
#' Outcome and optional random-effect columns are read from
#' \code{sample_data(psIN)}. The function automatically chooses the modeling
#' engine from outcome type and random-effect status:
#'
#' \itemize{
#'   \item no random effect: \code{glmnet::cv.glmnet()} with family
#'     \code{gaussian}, \code{binomial}, or \code{poisson}
#'   \item random effect present: \code{lme4::lmer()} for continuous outcomes
#'     and \code{lme4::glmer()} for binary/count outcomes
#' }
#'
#' @param psIN A \code{phyloseq} object.
#' @param outcome Character scalar; outcome column in \code{sample_data(psIN)}.
#' @param taxrank Character taxonomic rank. Defaults to \code{"ASV"}.
#' @param meta_vars Optional character vector of sample-level covariates to add
#'   to the predictor matrix.
#' @param random_effect Optional character scalar; grouping variable in
#'   \code{sample_data(psIN)} for mixed models.
#' @param outcome_type One of \code{"auto"}, \code{"continuous"},
#'   \code{"binary"}, or \code{"count"}.
#' @param family One of \code{"auto"}, \code{"gaussian"}, \code{"binomial"},
#'   or \code{"poisson"}.
#' @param transform One of \code{"log_rel"}, \code{"relative"}, or
#'   \code{"count"}.
#' @param top_n Optional integer; keep the top N taxa by mean abundance before
#'   model fitting. \code{NULL} keeps all taxa.
#' @param prevalence_min Optional prevalence filter on relative abundance
#'   presence rate. Defaults to \code{0}.
#' @param abundance_min Optional mean relative abundance filter. Defaults to
#'   \code{0}.
#' @param alpha Elastic-net mixing parameter for glmnet models.
#' @param nfolds Number of folds for \code{cv.glmnet()}.
#' @param standardize Logical; passed to \code{cv.glmnet()}.
#' @param validation One of \code{"apparent"} or \code{"oof"}. \code{NULL}
#'   is treated as \code{"apparent"}.
#' @param score_scale One of \code{"response"} or \code{"link"} for returned
#'   scores from glmnet.
#' @param score_method One of \code{"weighted_sum"} or \code{"predict"}.
#'   \code{"weighted_sum"} is the default and reproduces the
#'   \code{compute_mrs()} strategy used in
#'   \code{20260329_CCM_Phase2_Trajectory_Deep.R}.
#' @param na_action Currently only \code{"complete"} is supported.
#'
#' @return A list of class \code{"Go_MRS_fit"}.
#'
#' @export
Go_MRS_fit <- function(psIN,
                       outcome,
                       taxrank = "ASV",
                       meta_vars = NULL,
                       random_effect = NULL,
                       outcome_type = c("auto", "continuous", "binary", "count"),
                       family = c("auto", "gaussian", "binomial", "poisson"),
                       transform = c("log_rel", "relative", "count"),
                       top_n = NULL,
                       prevalence_min = 0,
                       abundance_min = 0,
                       alpha = 0.5,
                       nfolds = 5,
                       standardize = TRUE,
                       validation = NULL,
                       score_scale = c("response", "link"),
                       score_method = "weighted_sum",
                       na_action = c("complete")) {

  `%||%` <- function(x, y) if (is.null(x) || length(x) == 0) y else x

  outcome_type <- match.arg(outcome_type)
  family <- match.arg(family)
  transform <- match.arg(transform)
  validation <- if (is.null(validation)) "apparent" else match.arg(validation, c("apparent", "oof"))
  score_scale <- match.arg(score_scale)
  score_method <- match.arg(score_method, c("weighted_sum", "predict", "legacy_mrs"))
  if (identical(score_method, "legacy_mrs")) score_method <- "weighted_sum"
  na_action <- match.arg(na_action)

  if (!inherits(psIN, "phyloseq")) stop("`psIN` must be a phyloseq object.")
  if (!requireNamespace("phyloseq", quietly = TRUE)) stop("Package `phyloseq` is required.")

  meta <- data.frame(phyloseq::sample_data(psIN), check.names = FALSE, stringsAsFactors = FALSE)
  if (!is.character(outcome) || length(outcome) != 1 || !outcome %in% names(meta)) {
    stop("`outcome` must be a sample_data column name in `psIN`.")
  }
  if (!is.null(random_effect) && !random_effect %in% names(meta)) {
    stop("`random_effect` must be a sample_data column name in `psIN`.")
  }

  detect_outcome_type <- function(y) {
    y2 <- y[!is.na(y)]
    if (is.factor(y2) || is.character(y2)) {
      if (length(unique(as.character(y2))) == 2) return("binary")
      stop("Factor/character outcomes currently must be binary.")
    }
    if (is.logical(y2)) return("binary")
    if (is.numeric(y2)) {
      u <- sort(unique(y2))
      if (length(u) == 2 && all(u %in% c(0, 1))) return("binary")
      if (all(y2 >= 0) && all(abs(y2 - round(y2)) < .Machine$double.eps^0.5) && length(u) > 2) {
        return("count")
      }
      return("continuous")
    }
    stop("Unsupported outcome type.")
  }

  infer_family <- function(outcome_type_detected) {
    switch(
      outcome_type_detected,
      continuous = "gaussian",
      binary = "binomial",
      count = "poisson",
      stop("Unsupported detected outcome type.")
    )
  }

  compute_metrics <- function(y, score, outcome_type_detected) {
    out <- list()
    if (outcome_type_detected == "binary") {
      y_bin <- as.integer(y)
      if (requireNamespace("pROC", quietly = TRUE)) {
        roc_obj <- pROC::roc(response = y_bin, predictor = score, quiet = TRUE)
        out$roc <- roc_obj
        out$auc <- as.numeric(pROC::auc(roc_obj))
      } else {
        out$roc <- NULL
        out$auc <- NA_real_
      }
      pred_class <- as.integer(score >= 0.5)
      out$accuracy <- mean(pred_class == y_bin, na.rm = TRUE)
    } else {
      y_num <- as.numeric(y)
      out$rmse <- sqrt(mean((y_num - score)^2, na.rm = TRUE))
      out$cor <- suppressWarnings(stats::cor(y_num, score, use = "complete.obs"))
      out$r2 <- suppressWarnings(out$cor^2)
    }
    out
  }

  format_subtitle <- function(engine, family_used, family_mode,
                              outcome_type_detected, outcome_type_mode,
                              random_effect_used, validation_used,
                              n_obs, n_features,
                              taxrank_used, transform_used, metrics) {
    metric_txt <- switch(
      outcome_type_detected,
      binary = paste0(
        "AUC = ",
        ifelse(is.null(metrics$auc) || is.na(metrics$auc), "NA", sprintf("%.3f", metrics$auc))
      ),
      paste0(
        "R2 = ",
        ifelse(is.null(metrics$r2) || is.na(metrics$r2), "NA", sprintf("%.3f", metrics$r2))
      )
    )

    paste0(
      "Model: ", engine,
      " | family = ", family_used, " [", family_mode, "]",
      " | taxrank = ", taxrank_used,
      " | transform = ", transform_used,
      " | outcome = ", outcome_type_detected, " [", outcome_type_mode, "]",
      " | random effect = ", random_effect_used %||% "none",
      " | validation = ", validation_used,
      " | features = ", n_features,
      " | n = ", n_obs,
      " | ", metric_txt
    )
  }

  make_tax_labels <- function(ps_obj, taxrank_used) {
    taxa_ids <- phyloseq::taxa_names(ps_obj)
    if (identical(taxrank_used, "ASV") || is.null(phyloseq::tax_table(ps_obj, errorIfNULL = FALSE))) {
      return(stats::setNames(taxa_ids, taxa_ids))
    }
    tt <- as.data.frame(phyloseq::tax_table(ps_obj), stringsAsFactors = FALSE)
    if (taxrank_used %in% colnames(tt)) {
      vals <- tt[[taxrank_used]]
      vals[is.na(vals) | vals == ""] <- taxa_ids[is.na(vals) | vals == ""]
      vals <- make.unique(as.character(vals))
      return(stats::setNames(vals, taxa_ids))
    }
    stats::setNames(taxa_ids, taxa_ids)
  }

  ps_work <- psIN
  if (!identical(taxrank, "ASV")) {
    tt <- phyloseq::tax_table(ps_work, errorIfNULL = FALSE)
    if (is.null(tt) || !(taxrank %in% colnames(tt))) {
      stop(sprintf("taxrank '%s' not found in tax_table(psIN).", taxrank))
    }
    ps_work <- phyloseq::tax_glom(
      ps_work,
      taxrank = taxrank,
      NArm = !identical(score_method, "weighted_sum")
    )
  }

  otu_mat <- as(phyloseq::otu_table(ps_work), "matrix")
  if (phyloseq::taxa_are_rows(ps_work)) otu_mat <- t(otu_mat)

  if (!nrow(otu_mat) || !ncol(otu_mat)) stop("No OTU data found in `psIN`.")

  relab <- sweep(otu_mat, 1, pmax(1e-12, rowSums(otu_mat)), "/")
  keep_taxa <- rep(TRUE, ncol(relab))
  if (prevalence_min > 0) {
    keep_taxa <- keep_taxa & (colMeans(relab > 0) >= prevalence_min)
  }
  if (abundance_min > 0) {
    keep_taxa <- keep_taxa & (colMeans(relab) >= abundance_min)
  }
  relab <- relab[, keep_taxa, drop = FALSE]
  otu_mat <- otu_mat[, keep_taxa, drop = FALSE]
  if (!ncol(relab)) stop("No taxa remain after prevalence/abundance filtering.")

  if (!is.null(top_n) && !identical(score_method, "weighted_sum")) {
    top_n <- min(as.integer(top_n), ncol(relab))
    keep_top <- order(colMeans(relab), decreasing = TRUE)[seq_len(top_n)]
    relab <- relab[, keep_top, drop = FALSE]
    otu_mat <- otu_mat[, keep_top, drop = FALSE]
  }

  feat_micro <- switch(
    transform,
    log_rel = log(relab * 100 + 1),
    relative = relab,
    count = otu_mat
  )

  if (!is.null(top_n) && identical(score_method, "weighted_sum")) {
    top_n <- min(as.integer(top_n), ncol(feat_micro))
    keep_top <- order(colMeans(feat_micro), decreasing = TRUE)[seq_len(top_n)]
    feat_micro <- feat_micro[, keep_top, drop = FALSE]
    relab <- relab[, keep_top, drop = FALSE]
    otu_mat <- otu_mat[, keep_top, drop = FALSE]
  }

  label_map <- make_tax_labels(ps_work, taxrank)
  feature_names_raw <- colnames(feat_micro)
  feature_names_lab <- unname(label_map[feature_names_raw])
  feature_names_lab[is.na(feature_names_lab) | feature_names_lab == ""] <- feature_names_raw[is.na(feature_names_lab) | feature_names_lab == ""]
  feature_names_lab <- make.unique(feature_names_lab)
  colnames(feat_micro) <- feature_names_lab

  meta_use <- meta
  meta_vars <- intersect(meta_vars %||% character(0), names(meta_use))
  keep_meta_cols <- unique(c(outcome, meta_vars, random_effect))
  meta_use <- meta_use[, keep_meta_cols, drop = FALSE]
  meta_use$.SampleID <- rownames(meta_use)

  feat_df <- data.frame(feat_micro, check.names = FALSE, stringsAsFactors = FALSE)
  feat_df$.SampleID <- rownames(feat_df)

  work_df <- merge(meta_use, feat_df, by = ".SampleID", all = FALSE, sort = FALSE)
  rownames(work_df) <- work_df$.SampleID

  y_raw <- work_df[[outcome]]
  outcome_type_detected <- if (outcome_type == "auto") detect_outcome_type(y_raw) else outcome_type
  family_used <- if (family == "auto") infer_family(outcome_type_detected) else family
  outcome_type_mode <- if (identical(outcome_type, "auto")) "auto" else "manual"
  family_mode <- if (identical(family, "auto")) "auto" else "manual"

  if (na_action == "complete") {
    work_df <- stats::na.omit(work_df)
  }
  if (!nrow(work_df)) stop("No complete rows available after NA filtering.")

  if (outcome_type_detected == "binary") {
    y_vals <- work_df[[outcome]]
    if (is.factor(y_vals) || is.character(y_vals)) {
      levs <- unique(as.character(y_vals))
      if (length(levs) != 2) stop("Binary outcome must have exactly 2 levels.")
      work_df[[outcome]] <- factor(as.character(y_vals), levels = levs)
      y_model <- as.integer(work_df[[outcome]] == levs[2])
      positive_class <- levs[2]
      negative_class <- levs[1]
    } else {
      y_model <- as.integer(work_df[[outcome]])
      positive_class <- "1"
      negative_class <- "0"
    }
  } else {
    y_model <- work_df[[outcome]]
    positive_class <- NULL
    negative_class <- NULL
  }

  feature_names_used <- setdiff(colnames(feat_micro), ".SampleID")
  if (length(meta_vars)) {
    feature_names_used <- c(feature_names_used, meta_vars)
  }

  has_random <- !is.null(random_effect) && nzchar(random_effect)

  if (!has_random) {
    if (!requireNamespace("glmnet", quietly = TRUE)) {
      stop("Package `glmnet` is required for cv.glmnet models.")
    }

    predictor_df <- work_df[, feature_names_used, drop = FALSE]
    predictor_df <- predictor_df[, vapply(predictor_df, function(col) {
      vals <- unique(col[!is.na(col)])
      length(vals) > 1
    }, logical(1)), drop = FALSE]
    if (!ncol(predictor_df)) stop("No informative predictors remain after filtering.")

    if (identical(score_method, "weighted_sum")) {
      if (outcome_type_detected != "binary") stop("score_method = 'weighted_sum' currently supports binary outcomes only.")
      if (length(meta_vars)) stop("score_method = 'weighted_sum' does not support `meta_vars`.")
      if (!is.null(random_effect) && nzchar(random_effect)) stop("score_method = 'weighted_sum' does not support `random_effect`.")
      if (identical(validation, "oof")) stop("score_method = 'weighted_sum' requires validation = 'apparent'.")
      x_mat <- as.matrix(predictor_df)
    } else {
      x_mat <- stats::model.matrix(~ . - 1, data = predictor_df)
    }
    cv_fit <- glmnet::cv.glmnet(
      x = x_mat,
      y = y_model,
      family = family_used,
      alpha = alpha,
      nfolds = nfolds,
      standardize = standardize
    )

    if (identical(score_method, "weighted_sum")) {
      coef_mat <- as.matrix(stats::coef(cv_fit, s = "lambda.min"))
      coef_df <- data.frame(
        feature = rownames(coef_mat),
        estimate = as.numeric(coef_mat[, 1]),
        stringsAsFactors = FALSE,
        row.names = NULL
      )
      coef_df <- coef_df[coef_df$feature != "(Intercept)" & coef_df$estimate != 0, , drop = FALSE]
      coef_df <- coef_df[order(abs(coef_df$estimate), decreasing = TRUE), , drop = FALSE]
      if (!nrow(coef_df)) stop("No non-zero coefficients remain for legacy MRS scoring.")
      feat_mat <- as.matrix(predictor_df[, coef_df$feature, drop = FALSE])
      weight_mat <- matrix(coef_df$estimate, ncol = 1, dimnames = list(coef_df$feature, "estimate"))
      score <- as.numeric(feat_mat %*% weight_mat)
    } else if (identical(validation, "oof")) {
      nfolds_oof <- max(2, min(nfolds, nrow(x_mat)))
      foldid <- sample(rep(seq_len(nfolds_oof), length.out = nrow(x_mat)))
      score <- rep(NA_real_, nrow(x_mat))

      for (k in seq_len(nfolds_oof)) {
        train_idx <- foldid != k
        test_idx <- foldid == k
        if (!any(test_idx)) next
        if (family_used == "binomial" && length(unique(y_model[train_idx])) < 2) next

        fit_k <- glmnet::cv.glmnet(
          x = x_mat[train_idx, , drop = FALSE],
          y = y_model[train_idx],
          family = family_used,
          alpha = alpha,
          nfolds = max(2, min(nfolds, sum(train_idx))),
          standardize = standardize
        )

        score[test_idx] <- as.numeric(
          stats::predict(
            fit_k,
            newx = x_mat[test_idx, , drop = FALSE],
            s = "lambda.min",
            type = score_scale
          )
        )
      }

      if (anyNA(score)) {
        warning("Some OOF predictions were unavailable; missing values were filled with apparent predictions.")
        apparent_score <- as.numeric(
          stats::predict(cv_fit, newx = x_mat, s = "lambda.min", type = score_scale)
        )
        score[is.na(score)] <- apparent_score[is.na(score)]
      }
    } else {
      score <- as.numeric(
        stats::predict(cv_fit, newx = x_mat, s = "lambda.min", type = score_scale)
      )
    }

    if (!identical(score_method, "weighted_sum")) {
      coef_mat <- as.matrix(stats::coef(cv_fit, s = "lambda.min"))
      coef_df <- data.frame(
        feature = rownames(coef_mat),
        estimate = as.numeric(coef_mat[, 1]),
        stringsAsFactors = FALSE,
        row.names = NULL
      )
      coef_df <- coef_df[coef_df$feature != "(Intercept)" & coef_df$estimate != 0, , drop = FALSE]
      coef_df <- coef_df[order(abs(coef_df$estimate), decreasing = TRUE), , drop = FALSE]
    }

    engine <- if (identical(score_method, "weighted_sum")) "weighted_sum_glmnet" else "cv.glmnet"
    formula_used <- stats::as.formula(
      paste0("~ ", paste(sprintf("`%s`", colnames(predictor_df)), collapse = " + "))
    )
    feature_names_used <- colnames(predictor_df)
    model_object <- cv_fit
  } else {
    if (identical(validation, "oof")) {
      stop("validation = 'oof' is not yet implemented for mixed models. Use validation = 'apparent' for lmer/glmer.")
    }
    if (!requireNamespace("lme4", quietly = TRUE)) {
      stop("Package `lme4` is required for mixed-effect models.")
    }

    predictor_df <- work_df[, feature_names_used, drop = FALSE]
    predictor_df <- predictor_df[, vapply(predictor_df, function(col) {
      vals <- unique(col[!is.na(col)])
      length(vals) > 1
    }, logical(1)), drop = FALSE]
    feature_names_used <- colnames(predictor_df)
    if (!length(feature_names_used)) stop("No informative predictors remain after filtering.")

    fixed_terms <- paste(sprintf("`%s`", feature_names_used), collapse = " + ")
    random_term <- paste0("(1 | `", random_effect, "`)")
    formula_text <- paste(sprintf("`%s`", outcome), "~", fixed_terms, "+", random_term)
    formula_used <- stats::as.formula(formula_text)

    if (outcome_type_detected == "continuous") {
      fit <- lme4::lmer(formula_used, data = work_df, REML = FALSE)
      engine <- "lmer"
      score <- as.numeric(stats::predict(fit, type = "response"))
    } else if (outcome_type_detected == "binary") {
      fit <- lme4::glmer(formula_used, data = work_df, family = stats::binomial())
      engine <- "glmer"
      score <- as.numeric(stats::predict(fit, type = "response"))
    } else if (outcome_type_detected == "count") {
      fit <- lme4::glmer(formula_used, data = work_df, family = stats::poisson())
      engine <- "glmer"
      score <- as.numeric(stats::predict(fit, type = "response"))
    } else {
      stop("Unsupported mixed-model outcome type.")
    }

    coef_mat <- as.matrix(summary(fit)$coefficients)
    coef_df <- data.frame(
      feature = rownames(coef_mat),
      estimate = coef_mat[, "Estimate"],
      stringsAsFactors = FALSE,
      row.names = NULL
    )
    coef_df <- coef_df[coef_df$feature != "(Intercept)", , drop = FALSE]
    coef_df <- coef_df[order(abs(coef_df$estimate), decreasing = TRUE), , drop = FALSE]
    model_object <- fit
  }

  metrics <- compute_metrics(
    y = y_model,
    score = score,
    outcome_type_detected = outcome_type_detected
  )

  pred_df <- data.frame(
    SampleID = rownames(work_df),
    observed = y_model,
    score = score,
    stringsAsFactors = FALSE
  )
  if (!is.null(random_effect)) pred_df[[random_effect]] <- work_df[[random_effect]]

  subtitle <- format_subtitle(
    engine = engine,
    family_used = family_used,
    family_mode = family_mode,
    outcome_type_detected = outcome_type_detected,
    outcome_type_mode = outcome_type_mode,
    random_effect_used = random_effect,
    validation_used = validation,
    n_obs = nrow(work_df),
    n_features = length(feature_names_used),
    taxrank_used = taxrank,
    transform_used = transform,
    metrics = metrics
  )

  out <- list(
    model = model_object,
    model_info = list(
      engine = engine,
      family = family_used,
      family_mode = family_mode,
      outcome = outcome,
      outcome_type = outcome_type_detected,
      outcome_type_mode = outcome_type_mode,
      random_effect = random_effect,
      validation = validation,
      taxrank = taxrank,
      transform = transform,
      score_method = score_method,
      features = feature_names_used,
      formula = formula_used,
      positive_class = positive_class,
      negative_class = negative_class,
      n = nrow(work_df)
    ),
    coefficients = coef_df,
    predictions = pred_df,
    metrics = metrics,
    subtitle = subtitle,
    data_used = work_df
  )

  class(out) <- c("Go_MRS_fit", class(out))
  out
}
