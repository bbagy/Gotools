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
#'   to the predictor matrix. In \code{score_method = "deep_mrs"}, these
#'   covariates are expanded with \code{model.matrix()} when needed and are
#'   included alongside microbial features.
#' @param random_effect Optional character scalar; grouping variable in
#'   \code{sample_data(psIN)} for mixed models.
#' @param outcome_type One of \code{"auto"}, \code{"continuous"},
#'   \code{"binary"}, or \code{"count"}.
#' @param family One of \code{"auto"}, \code{"gaussian"}, \code{"binomial"},
#'   or \code{"poisson"}.
#' @param transform Feature transform applied before model fitting. One of
#'   \code{"clr"} (default; centered log-ratio, compositionally correct),
#'   \code{"log_rel"} (\code{log(relab * 100 + 1)}), \code{"relative"}
#'   (raw relative abundance), or \code{"count"} (raw OTU counts).
#' @param top_n Optional integer; retain at most the top N taxa by variance
#'   (\code{"clr"}) or mean abundance (other transforms) before model fitting.
#'   \code{NULL} (default) passes all prevalence/abundance-filtered taxa to
#'   glmnet, which then applies L1 sparsity internally.
#' @param prevalence_min Optional prevalence filter (proportion of samples with
#'   non-zero relative abundance). Defaults to \code{0}.
#' @param abundance_min Optional mean relative abundance filter. Defaults to
#'   \code{0}.
#' @param alpha Elastic-net mixing parameter passed to \code{cv.glmnet()}.
#'   \code{0.5} (default) balances L1 sparsity and L2 ridge shrinkage.
#' @param nfolds Number of CV folds for \code{cv.glmnet()}. Defaults to
#'   \code{5}.
#' @param seed Integer seed for reproducible fold assignment when
#'   \code{foldid} is not supplied. Defaults to \code{1234}.
#' @param foldid Optional integer vector passed to \code{cv.glmnet()} to fix
#'   fold membership. Length must match the analysis sample count after
#'   filtering/NA removal.
#' @param standardize Logical; passed to \code{cv.glmnet()}. Defaults to
#'   \code{TRUE}.
#' @param validation One of \code{"oof"} (default; out-of-fold cross-validation
#'   for unbiased score estimation) or \code{"apparent"} (train-set scoring,
#'   inflates AUC).
#' @param score_scale One of \code{"response"} or \code{"link"} for returned
#'   scores from glmnet predict (used only when
#'   \code{score_method = "predict"}).
#' @param score_method One of \code{"deep_mrs"} (default; elastic net
#'   coefficients as weights, weighted-sum MRS), \code{"weighted_sum"} or
#'   \code{"legacy_mrs"} (backward-compatible aliases for \code{"deep_mrs"}),
#'   or \code{"predict"} (uses \code{stats::predict.cv.glmnet()} directly).
#' @param project Optional project name used to save the fitted object.
#' @param name Optional file tag for saving the fitted object. When both
#'   \code{project} and \code{name} are provided, the fitted object is saved to
#'   \code{<project>_YYMMDD/table/MRS_tab/}.
#' @param na_action Currently only \code{"complete"} is supported.
#'
#' @return A list of class \code{"Go_MRS_fit"}.
#'
#' @export
Go_MRS_fit <- function(psIN,
                       outcome,
                       taxrank       = "ASV",
                       meta_vars     = NULL,
                       random_effect = NULL,
                       outcome_type  = c("auto", "continuous", "binary", "count"),
                       family        = c("auto", "gaussian", "binomial", "poisson"),
                       transform     = c("clr", "log_rel", "relative", "count"),
                       top_n         = NULL,
                       prevalence_min = 0,
                       abundance_min  = 0,
                       alpha         = 0.5,
                       nfolds        = 5,
                       seed          = 1234,
                       foldid        = NULL,
                       standardize   = TRUE,
                       validation    = c("oof", "apparent"),
                       score_scale   = c("response", "link"),
                       score_method  = "deep_mrs",
                       project       = NULL,
                       name          = NULL,
                       na_action     = c("complete"),
                       verbose       = TRUE) {

  `%||%` <- function(x, y) if (is.null(x) || length(x) == 0) y else x
  log_msg <- function(...) if (isTRUE(verbose)) message(...)
  clean_tag <- function(x) {
    x <- gsub("[^A-Za-z0-9._-]+", "_", x)
    x <- gsub("_+", "_", x)
    x <- gsub("^_|_$", "", x)
    x
  }

  outcome_type <- match.arg(outcome_type)
  family <- match.arg(family)
  transform <- match.arg(transform)
  validation <- match.arg(validation)
  score_scale <- match.arg(score_scale)
  score_method <- match.arg(score_method, c("deep_mrs", "weighted_sum", "predict", "legacy_mrs"))
  if (score_method %in% c("weighted_sum", "legacy_mrs")) score_method <- "deep_mrs"
  na_action <- match.arg(na_action)
  if (!is.null(seed)) seed <- as.integer(seed)[1]

  if (!inherits(psIN, "phyloseq")) stop("`psIN` must be a phyloseq object.")
  if (!requireNamespace("phyloseq", quietly = TRUE)) stop("Package `phyloseq` is required.")

  meta <- data.frame(phyloseq::sample_data(psIN), check.names = FALSE, stringsAsFactors = FALSE)
  log_msg(sprintf("[Go_MRS_fit] start: outcome='%s' taxrank='%s' n_samples=%d", outcome, taxrank, nrow(meta)))
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

  safe_coef_df <- function(cv_fit_obj) {
    for (s_val in c("lambda.min", "lambda.1se")) {
      cm <- as.matrix(stats::coef(cv_fit_obj, s = s_val))
      df <- data.frame(feature = rownames(cm), estimate = as.numeric(cm[, 1]),
                       stringsAsFactors = FALSE, row.names = NULL)
      df <- df[df$feature != "(Intercept)" & df$estimate != 0, , drop = FALSE]
      if (nrow(df) > 0) {
        if (!identical(s_val, "lambda.min"))
          warning(sprintf("[Go_MRS_fit] lambda.min gave 0 features; using %s.", s_val))
        return(df[order(abs(df$estimate), decreasing = TRUE), , drop = FALSE])
      }
    }
    nz_idx <- which(cv_fit_obj$nzero > 0)
    if (length(nz_idx) > 0) {
      lam_nz <- cv_fit_obj$lambda[nz_idx[length(nz_idx)]]
      cm <- as.matrix(stats::coef(cv_fit_obj, s = lam_nz))
      df <- data.frame(feature = rownames(cm), estimate = as.numeric(cm[, 1]),
                       stringsAsFactors = FALSE, row.names = NULL)
      df <- df[df$feature != "(Intercept)" & df$estimate != 0, , drop = FALSE]
      if (nrow(df) > 0) {
        warning("[Go_MRS_fit] lambda.min/1se gave 0 features; using most-regularized non-zero lambda.")
        return(df[order(abs(df$estimate), decreasing = TRUE), , drop = FALSE])
      }
    }
    NULL
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
        roc_obj <- pROC::roc(response = y_bin, predictor = score,
                             direction = "<", quiet = TRUE)
        out$roc <- roc_obj
        out$auc <- as.numeric(pROC::auc(roc_obj))
        auc_ci <- tryCatch(as.numeric(pROC::ci.auc(roc_obj)), error = function(e) NULL)
        out$auc_ci <- auc_ci
        youden_idx <- which.max(roc_obj$sensitivities + roc_obj$specificities - 1)
        out$threshold <- roc_obj$thresholds[youden_idx]
      } else {
        out$roc <- NULL
        out$auc <- NA_real_
        out$auc_ci <- NULL
        out$threshold <- 0.5
      }
      thresh <- if (!is.null(out$threshold) && is.finite(out$threshold)) out$threshold else 0.5
      pred_class <- as.integer(score >= thresh)
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
                              taxrank_used, transform_used, metrics,
                              covars_used) {
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
      " | covars = ", if (length(covars_used)) paste(covars_used, collapse = ",") else "none",
      " | validation = ", validation_used,
      " | features = ", n_features,
      " | n = ", n_obs,
      " | ", metric_txt
    )
  }

  make_tax_labels <- function(ps_obj, taxrank_used) {
    taxa_ids <- phyloseq::taxa_names(ps_obj)
    tt <- phyloseq::tax_table(ps_obj, errorIfNULL = FALSE)
    if (is.null(tt)) {
      return(stats::setNames(taxa_ids, taxa_ids))
    }
    tt <- as.data.frame(tt, stringsAsFactors = FALSE)
    if (identical(taxrank_used, "ASV")) {
      genus_vals <- if ("Genus" %in% colnames(tt)) tt$Genus else rep(NA_character_, nrow(tt))
      species_vals <- if ("Species" %in% colnames(tt)) tt$Species else rep(NA_character_, nrow(tt))
      family_vals <- if ("Family" %in% colnames(tt)) tt$Family else rep(NA_character_, nrow(tt))
      lab <- ifelse(
        !is.na(species_vals) & species_vals != "",
        paste(genus_vals, species_vals),
        ifelse(!is.na(genus_vals) & genus_vals != "",
               paste(genus_vals, "sp."),
               ifelse(!is.na(family_vals) & family_vals != "",
                      paste(family_vals, "taxon"),
                      taxa_ids))
      )
      lab[is.na(lab) | lab == ""] <- taxa_ids[is.na(lab) | lab == ""]
      lab <- make.unique(as.character(lab))
      return(stats::setNames(lab, taxa_ids))
    }
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
      NArm = score_method != "deep_mrs"
    )
  }

  otu_mat <- as(phyloseq::otu_table(ps_work), "matrix")
  if (phyloseq::taxa_are_rows(ps_work)) otu_mat <- t(otu_mat)

  if (!nrow(otu_mat) || !ncol(otu_mat)) stop("No OTU data found in `psIN`.")
  log_msg(sprintf("[Go_MRS_fit] OTU matrix: n_samples=%d n_features=%d", nrow(otu_mat), ncol(otu_mat)))

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
  log_msg(sprintf("[Go_MRS_fit] after prevalence/abundance filter: n_features=%d", ncol(relab)))

  if (!is.null(top_n) && score_method != "deep_mrs") {
    top_n <- min(as.integer(top_n), ncol(relab))
    keep_top <- order(colMeans(relab), decreasing = TRUE)[seq_len(top_n)]
    relab   <- relab[, keep_top, drop = FALSE]
    otu_mat <- otu_mat[, keep_top, drop = FALSE]
  }

  feat_micro <- switch(
    transform,
    log_rel  = log(relab * 100 + 1),
    relative = relab,
    count    = otu_mat,
    clr      = {
      p_clr <- relab + 1e-6
      log(p_clr) - rowMeans(log(p_clr))
    }
  )

  if (!is.null(top_n) && score_method == "deep_mrs") {
    top_n <- min(as.integer(top_n), ncol(feat_micro))
    rank_vec <- if (identical(transform, "clr")) {
      apply(feat_micro, 2, var)
    } else {
      colMeans(feat_micro)
    }
    keep_top <- order(rank_vec, decreasing = TRUE)[seq_len(top_n)]
    feat_micro <- feat_micro[, keep_top, drop = FALSE]
    relab <- relab[, keep_top, drop = FALSE]
    otu_mat <- otu_mat[, keep_top, drop = FALSE]
  }
  log_msg(sprintf("[Go_MRS_fit] after top_n: n_features=%d", ncol(feat_micro)))

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
  log_msg(sprintf("[Go_MRS_fit] after NA filter: n_samples=%d", nrow(work_df)))
  if (!is.null(foldid)) {
    if (length(foldid) != nrow(work_df)) {
      stop("`foldid` length must match the number of analysis samples after filtering/NA removal.")
    }
    foldid <- as.integer(foldid)
    if (any(is.na(foldid))) stop("`foldid` must not contain NA.")
    if (length(unique(foldid)) < 2) stop("`foldid` must contain at least 2 unique folds.")
  } else if (!is.null(seed)) {
    set.seed(seed)
    log_msg(sprintf("[Go_MRS_fit] using seed=%d", seed))
  }

  if (outcome_type_detected == "binary") {
    y_vals <- work_df[[outcome]]
    if (is.factor(y_vals) || is.character(y_vals)) {
      levs <- if (is.factor(y_vals)) levels(y_vals) else unique(as.character(y_vals))
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

  feature_names_used <- colnames(feat_micro)
  if (length(meta_vars)) {
    feature_names_used <- c(feature_names_used, meta_vars)
  }

  has_random <- !is.null(random_effect) && nzchar(random_effect)

  if (!has_random) {
    log_msg(sprintf("[Go_MRS_fit] fitting glmnet: family=%s validation=%s score_method=%s", family_used, validation, score_method))
    if (!requireNamespace("glmnet", quietly = TRUE)) {
      stop("Package `glmnet` is required for cv.glmnet models.")
    }

    predictor_df <- work_df[, feature_names_used, drop = FALSE]
    predictor_df <- predictor_df[, vapply(predictor_df, function(col) {
      vals <- unique(col[!is.na(col)])
      length(vals) > 1
    }, logical(1)), drop = FALSE]
    if (!ncol(predictor_df)) stop("No informative predictors remain after filtering.")

    if (score_method == "deep_mrs") {
      if (outcome_type_detected != "binary") stop("score_method = 'deep_mrs' currently supports binary outcomes only.")
      if (!is.null(random_effect) && nzchar(random_effect)) stop("score_method = 'deep_mrs' does not support `random_effect`.")
      micro_cols <- intersect(setdiff(colnames(feat_micro), ".SampleID"), colnames(predictor_df))
      meta_cols <- setdiff(colnames(predictor_df), micro_cols)
      x_micro <- if (length(micro_cols)) as.matrix(predictor_df[, micro_cols, drop = FALSE]) else NULL
      x_meta <- if (length(meta_cols)) stats::model.matrix(~ . - 1, data = predictor_df[, meta_cols, drop = FALSE]) else NULL
      if (is.null(x_micro) && is.null(x_meta)) stop("No predictors remain for deep_mrs.")
      if (is.null(x_micro)) {
        x_mat <- x_meta
      } else if (is.null(x_meta)) {
        x_mat <- x_micro
      } else {
        x_mat <- cbind(x_micro, x_meta)
      }
    } else {
      x_mat <- stats::model.matrix(~ . - 1, data = predictor_df)
    }
    cv_fit <- glmnet::cv.glmnet(
      x = x_mat,
      y = y_model,
      family = family_used,
      alpha = alpha,
      nfolds = nfolds,
      foldid = foldid,
      standardize = standardize
    )

    if (score_method == "deep_mrs" && identical(validation, "oof")) {
      nfolds_oof <- max(2, min(nfolds, nrow(x_mat)))
      foldid_oof <- if (!is.null(foldid)) {
        foldid
      } else {
        sample(rep(seq_len(nfolds_oof), length.out = nrow(x_mat)))
      }
      folds_oof <- sort(unique(foldid_oof))
      score <- rep(NA_real_, nrow(x_mat))

      for (k in folds_oof) {
        train_idx <- foldid_oof != k
        test_idx  <- foldid_oof == k
        if (!any(test_idx)) next
        if (length(unique(y_model[train_idx])) < 2) next

        foldid_train <- NULL
        if (!is.null(foldid)) {
          foldid_train <- foldid_oof[train_idx]
          foldid_train <- match(foldid_train, unique(foldid_train))
        }

        fit_k <- glmnet::cv.glmnet(
          x = x_mat[train_idx, , drop = FALSE],
          y = y_model[train_idx],
          family = family_used,
          alpha = alpha,
          nfolds = max(2, min(nfolds, sum(train_idx))),
          foldid = foldid_train,
          standardize = standardize
        )
        coef_k <- as.matrix(stats::coef(fit_k, s = "lambda.min"))
        coef_k_df <- data.frame(
          feature  = rownames(coef_k),
          estimate = as.numeric(coef_k[, 1]),
          stringsAsFactors = FALSE
        )
        coef_k_df <- coef_k_df[coef_k_df$feature != "(Intercept)" & coef_k_df$estimate != 0, , drop = FALSE]
        if (!nrow(coef_k_df)) next
        avail_feats <- intersect(coef_k_df$feature, colnames(x_mat))
        if (!length(avail_feats)) next
        coef_k_df <- coef_k_df[coef_k_df$feature %in% avail_feats, , drop = FALSE]
        wt_k <- stats::setNames(coef_k_df$estimate, coef_k_df$feature)
        score[test_idx] <- as.numeric(x_mat[test_idx, names(wt_k), drop = FALSE] %*% wt_k)
      }

      if (anyNA(score)) {
        warning("[Go_MRS_fit] Some OOF deep_mrs folds produced no features; filling NAs with apparent weighted-sum.")
        coef_app_df <- safe_coef_df(cv_fit)
        if (!is.null(coef_app_df) && nrow(coef_app_df) > 0) {
          wt_app    <- stats::setNames(coef_app_df$estimate, coef_app_df$feature)
          avail_app <- intersect(names(wt_app), colnames(x_mat))
          app_score <- as.numeric(x_mat[, avail_app, drop = FALSE] %*% wt_app[avail_app])
          score[is.na(score)] <- app_score[is.na(score)]
        }
        if (anyNA(score)) {
          warning("[Go_MRS_fit] Apparent fallback also unavailable; remaining NAs set to 0.")
          score[is.na(score)] <- 0
        }
      }

      coef_df <- safe_coef_df(cv_fit)
      if (is.null(coef_df) || !nrow(coef_df)) {
        warning("[Go_MRS_fit] OOF: no features in final model coefficient report; coef table will be empty.")
        coef_df <- data.frame(feature = character(0), estimate = numeric(0),
                              stringsAsFactors = FALSE)
      }
    } else if (score_method == "deep_mrs") {
      coef_df <- safe_coef_df(cv_fit)
      if (is.null(coef_df) || !nrow(coef_df)) {
        stop("[Go_MRS_fit] No features selected by glmnet at any lambda. Check data quality, sample size, or increase top_n.")
      }
      avail_feats <- intersect(coef_df$feature, colnames(x_mat))
      coef_df     <- coef_df[coef_df$feature %in% avail_feats, , drop = FALSE]
      feat_mat    <- x_mat[, avail_feats, drop = FALSE]
      weight_mat  <- matrix(coef_df$estimate, ncol = 1,
                            dimnames = list(coef_df$feature, "estimate"))
      score <- as.numeric(feat_mat %*% weight_mat)
    } else if (identical(validation, "oof")) {
      nfolds_oof <- max(2, min(nfolds, nrow(x_mat)))
      foldid_oof <- if (!is.null(foldid)) {
        foldid
      } else {
        sample(rep(seq_len(nfolds_oof), length.out = nrow(x_mat)))
      }
      folds_oof <- sort(unique(foldid_oof))
      score <- rep(NA_real_, nrow(x_mat))

      for (k in folds_oof) {
        train_idx <- foldid_oof != k
        test_idx <- foldid_oof == k
        if (!any(test_idx)) next
        if (family_used == "binomial" && length(unique(y_model[train_idx])) < 2) next

        foldid_train <- NULL
        if (!is.null(foldid)) {
          foldid_train <- foldid_oof[train_idx]
          foldid_train <- match(foldid_train, unique(foldid_train))
        }

        fit_k <- glmnet::cv.glmnet(
          x = x_mat[train_idx, , drop = FALSE],
          y = y_model[train_idx],
          family = family_used,
          alpha = alpha,
          nfolds = max(2, min(nfolds, sum(train_idx))),
          foldid = foldid_train,
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

    if (score_method != "deep_mrs") {
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

    engine <- if (score_method == "deep_mrs") "mrs_glmnet" else "cv.glmnet"
    formula_used <- stats::as.formula(
      paste0("~ ", paste(sprintf("`%s`", colnames(predictor_df)), collapse = " + "))
    )
    feature_names_used <- colnames(x_mat)
    model_object <- cv_fit
  } else {
    log_msg(sprintf("[Go_MRS_fit] fitting mixed model: family=%s validation=%s", family_used, validation))
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
    metrics = metrics,
    covars_used = meta_vars
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
      seed = seed,
      foldid = foldid,
      taxrank = taxrank,
      transform = transform,
      meta_vars = meta_vars,
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
  if (!is.null(project) && !is.null(name)) {
    out_dirs <- Go_path(project = project, pdf = "no", table = "yes", path = NULL)
    out_dir <- file.path(out_dirs$tab, "MRS_tab")
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    save_path <- file.path(out_dir, paste0(clean_tag(name), ".rds"))
    saveRDS(out, save_path)
    out$saved_rds <- save_path
    log_msg(sprintf("[Go_MRS_fit] saved RDS: %s", save_path))
  }
  out
}
