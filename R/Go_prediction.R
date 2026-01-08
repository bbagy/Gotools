#' Go_prediction
#'
#' Unified predictor for microbiome (+ optional clinical) using Random Forest
#' or XGBoost with StudyID-aware cross-validation/holdout, random search tuning,
#' and standardized outputs.
#'
#' @description
#' Builds a binary classifier from a \code{phyloseq} object and (optionally)
#' clinical covariates. You can toggle inclusion of microbiome features via
#' \code{bacteriaSet}, choose the algorithm (\code{"randomforest"} or \code{"xgboost"}),
#' perform K-fold CV with subject-aware folds, run a random hyperparameter search,
#' and (optionally) evaluate a 75/25 holdout split. The function writes dated
#' output folders mirroring \code{Go_randomforest}/\code{Go_xgboost}.
#'
#' @param psIN A \code{phyloseq} object.
#' @param project Character; project prefix used to create an output directory
#'   \verb{<project_YYMMDD>/Randomforest} or \verb{<project_YYMMDD>/XGBoost}.
#' @param outcome Character; name of the binary outcome column in sample metadata.
#'   Must contain levels specified in \code{orders}.
#' @param method Character; one of \code{c("randomforest","xgboost")}.
#' @param bacteriaSet Logical; include microbiome features (relative abundances
#'   after filtering). If \code{FALSE}, only clinical features are used.
#' @param testSet Logical; if \code{FALSE} (default) runs CV-only. If \code{TRUE},
#'   performs a 75/25 StudyID-aware holdout split with inner CV on the train set.
#' @param clinical_vari Character vector of metadata columns to include as
#'   covariates (optional). Character/factor covariates are one-hot encoded,
#'   NAs imputed by column median, then scaled.
#' @param StudyID_col Character; subject/group id column. If missing, rownames
#'   are used. Default \code{"StudyID"}.
#' @param taxrank Either \code{"ASV"} to keep features as is, or a taxonomic rank
#'   present in \code{tax_table(psIN)} (e.g., \code{"Genus"}). Default \code{"ASV"}.
#' @param prev_min Numeric in \([0,1]\); minimum prevalence threshold (presence on
#'   relative abundance). Default \code{0.01}.
#' @param relab_min Numeric; minimum mean relative abundance threshold. Default \code{1e-4}.
#' @param n_folds Integer; number of CV folds. Default \code{5}.
#' @param seed Integer; random seed. Default \code{123}.
#' @param n_candidates Integer; random search candidates. Default \code{40}.
#' @param num.trees Integer; RF trees / XGB max \code{nrounds} (early stopping used
#'   for XGB). Default \code{1000}.
#' @param orders Character vector of length 2 giving outcome level order
#'   (e.g., \code{c("Control","Case")}); positive class is \code{orders[2]}.
#'
#' @details
#' \strong{StudyID-aware CV:} When \code{StudyID_col} varies across rows,
#' folds keep samples from the same subject together and preserve class balance
#' at the group level when possible; otherwise, standard stratified CV is used.
#'
#' \strong{Filtering:} If \code{bacteriaSet=TRUE}, taxa are kept if they pass
#' \code{prev_min} OR \code{relab_min}.
#'
#' \strong{Tuning:} RF tunes \code{mtry}, \code{min.node.size}, \code{sample.fraction}.
#' XGB tunes \code{eta}, \code{max_depth}, \code{min_child_weight}, \code{subsample},
#' \code{colsample_bytree}, \code{gamma}, \code{lambda}, \code{alpha}, and
#' \code{early_stopping_rounds}.
#'
#' \strong{Outputs (side effects):}
#' \itemize{
#'   \item \code{random_search_results.csv}
#'   \item \code{predictions.csv} (OOF for CV-only; Train/Test for holdout)
#'   \item model + metadata: \code{rf_final_model.rds}/\code{xgb_final_model.rds},
#'         \code{rf_meta.rds}/\code{xgb_meta.rds}
#'   \item feature importance with direction (\code{importance_feature.csv}) and
#'         aggregated by taxon (\code{importance_taxon.csv})
#'   \item ROC/PR PNGs (OOF or Train/Test variants)
#' }
#'
#' @return (Invisibly) a list containing mode (\code{"CV-only"} or \code{"Holdout+CV"}),
#'   metrics (AUC/AUPRC for OOF or Train/Test), \code{best_param}, and \code{outdir}.
#'
#' @examples
#' \dontrun{
#' res <- Go_prediction(
#'   psIN = ps, project = "IBD", outcome = "Status",
#'   method = "randomforest", bacteriaSet = TRUE,
#'   clinical_vari = c("Age","BMI","Sex"),
#'   StudyID_col = "SubjectID", taxrank = "Genus",
#'   orders = c("Control","Case")
#' )
#' res$AUC
#' }
#'
#' @seealso \code{\link{Go_randomforest}}, \code{\link{Go_xgboost}}
#'
#' @name Go_prediction
#' @export

suppressPackageStartupMessages({
  library(phyloseq)
  library(ranger)
  library(xgboost)
  library(pROC)
  library(PRROC)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(stringr)
})

Go_prediction <- function(
    psIN,
    project,
    outcome,
    method = c("randomforest","xgboost"),
    name = NULL,
    bacteriaSet = TRUE,             # 마이크로바이옴 피처 포함 여부
    testSet = FALSE,                # TRUE면 75/25 holdout + inner-CV
    clinical_vari = character(0),
    StudyID_col   = "StudyID",
    taxrank       = "ASV",
    prev_min      = 0.01,
    relab_min     = 1e-4,
    n_folds       = 5,
    seed          = 123,
    n_candidates  = 40,
    num.trees     = 1000,
    orders        = c("Control","Case")
){
  method <- match.arg(method)
  set.seed(seed)
  stopifnot(inherits(psIN, "phyloseq"))

  ## --- 내부 유틸 --------------------------------------------------------------
  .map_direction_safe <- function(feature_names, dir_lookup) {
    # 1차: 원본 이름 매칭
    res <- unname(dir_lookup[feature_names])
    miss <- is.na(res)
    if (any(miss)) {
      # 2차: make.names()로 정규화 이름 매칭
      dn <- dir_lookup
      names(dn) <- make.names(names(dn))
      res[miss] <- unname(dn[make.names(feature_names[miss])])
    }
    res[is.na(res)] <- "Neutral"
    res
  }
  safe_spearman <- function(a, b){
    if (is.null(a) || is.null(b)) return(NA_real_)
    if (all(is.na(a)) || all(is.na(b))) return(NA_real_)
    if (length(unique(a[!is.na(a)])) < 2) return(NA_real_)
    suppressWarnings(cor(a, b, method="spearman", use="pairwise.complete.obs"))
  }
  get_auc <- function(obs_fac, prob, orders){
    prob <- as.numeric(prob)
    if (all(is.na(prob)) || length(unique(stats::na.omit(prob))) < 2) return(0.5)
    as.numeric(pROC::auc(pROC::roc(obs_fac, prob, levels=orders, direction="<", quiet=TRUE)))
  }
  make_group_strat_folds <- function(groups, ybin, K=5, seed=123){
    set.seed(seed)
    df <- data.frame(g=groups, y=ybin)
    rep_lab <- aggregate(y ~ g, df, function(z) round(mean(z)))
    u0 <- rep_lab$g[rep_lab$y==0]; u1 <- rep_lab$g[rep_lab$y==1]
    if (length(u0) == 0 || length(u1) == 0) {
      labs <- unique(groups)
      foldg <- sample(rep(1:K, length.out = length(labs)))
      return(lapply(1:K, function(k) which(foldg[match(groups, labs)] == k)))
    }
    k0 <- if (length(u0) < K) max(2, length(u0)) else K
    k1 <- if (length(u1) < K) max(2, length(u1)) else K
    f0 <- sample(rep(1:k0, length.out=length(u0)))
    f1 <- sample(rep(1:k1, length.out=length(u1)))
    fmap <- rbind(data.frame(g=u0, f=f0), data.frame(g=u1, f=f1))
    ff <- fmap$f[match(groups, fmap$g)]
    ff[is.na(ff)] <- sample(1:max(fmap$f), sum(is.na(ff)), replace=TRUE)
    lapply(1:max(ff), function(k) which(ff == k))
  }
  make_strat_folds <- function(ybin, K=5, seed=123){
    set.seed(seed)
    idx0 <- which(ybin==0); idx1 <- which(ybin==1)
    k0 <- if (length(idx0) < K) max(2, length(idx0)) else K
    k1 <- if (length(idx1) < K) max(2, length(idx1)) else K
    f0 <- sample(rep(1:k0, length.out=length(idx0)))
    f1 <- sample(rep(1:k1, length.out=length(idx1)))
    ff <- integer(length(ybin)); ff[idx0] <- f0; ff[idx1] <- f1
    lapply(1:max(ff), function(k) which(ff == k))
  }



  # ---- CV vs OOF 요약 저장 ----
  log_cv_oof_summary <- function(outdir, model_name,
                                 CV_best_AUC, CV_best_iter = NA_integer_,
                                 oof_prob, yfac, positive_class, orders) {
    if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
    # OOF 성능
    AUC_oof <- get_auc(yfac, oof_prob, orders)
    pr <- PRROC::pr.curve(scores.class0 = oof_prob[yfac==positive_class],
                          scores.class1 = oof_prob[yfac!=positive_class], curve=TRUE)
    AUPRC_oof <- pr$auc.integral

    # 저장
    sum_df <- data.frame(
      Model        = model_name,
      CV_best_AUC  = CV_best_AUC,
      CV_best_iter = CV_best_iter,
      OOF_AUC      = AUC_oof,
      OOF_AUPRC    = AUPRC_oof
    )

    write.csv(sum_df, file.path(outdir, paste0(model_name, "_cv_oof_summary.csv")), row.names = FALSE)

    # 그림(선택)
    png(file.path(outdir, paste0("ROC_OOF_", model_name, ".png")), width=900, height=900, res=130)
    plot.roc(pROC::roc(yfac, oof_prob, levels=orders, direction="<", quiet=TRUE),
             main=sprintf("%s OOF ROC (AUC=%.3f)", model_name, AUC_oof), lwd=2); abline(0,1,lty=2); dev.off()
    png(file.path(outdir, paste0("PR_OOF_", model_name, ".png")), width=900, height=900, res=130)
    plot(pr, main=sprintf("%s OOF PR (AUPRC=%.3f)", model_name, AUPRC_oof)); dev.off()

    invisible(list(CV_best_AUC=CV_best_AUC, CV_best_iter=CV_best_iter,
                   OOF_AUC=AUC_oof, OOF_AUPRC=AUPRC_oof))
  }
  ## --- 0) 메타 준비: 레벨/StudyID 고정 --------------------------------------
  meta0 <- data.frame(sample_data(psIN), check.names = FALSE, stringsAsFactors = FALSE)
  stopifnot(outcome %in% names(meta0))
  if (!(StudyID_col %in% names(meta0)) || all(is.na(meta0[[StudyID_col]]))) {
    message(sprintf("[Info] '%s'가 없어 rownames를 사용합니다.", StudyID_col))
    meta0[[StudyID_col]] <- rownames(meta0)
  }
  yy0 <- meta0[[outcome]]; if (!is.factor(yy0)) yy0 <- factor(yy0)
  if (!all(orders %in% levels(yy0))) {
    stop(sprintf("[Error] outcome(levels=%s)에 orders(%s) 중 일부가 없습니다.",
                 paste(levels(yy0), collapse=", "), paste(orders, collapse=", ")))
  }
  yy0 <- factor(yy0, levels = orders)
  meta0[[outcome]] <- yy0
  sample_data(psIN)[[outcome]] <- yy0
  if (!(StudyID_col %in% colnames(sample_data(psIN))))
    sample_data(psIN)[[StudyID_col]] <- meta0[[StudyID_col]]
  positive_class <- orders[2]

  ## --- 1) 출력 폴더 ----------------------------------------------------------
  today <- format(Sys.Date(), "%y%m%d")
  method_label <- if (method=="randomforest") "randomforest" else "xgboost"

  if (!is.null(name) && nzchar(name)) {
    method_label <- paste0(method_label, "_", name)
  }

  root <- file.path(
    sprintf("%s_%s", project, today),
    method_label
  )
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  message("[OutDir] ", root)

  ## --- 2) phyloseq → feature table ------------------------------------------
  ps <- psIN
  if (!identical(taxrank, "ASV")) {
    if (is.null(tax_table(ps)) || !(taxrank %in% colnames(tax_table(ps)))) {
      stop(sprintf("[Error] taxrank '%s' not found in tax_table.", taxrank))
    }
    ps <- tax_glom(ps, taxrank = taxrank, NArm = TRUE)
  }
  relab <- NULL
  if (bacteriaSet) {
    OTU <- as(otu_table(ps), "matrix")
    if (taxa_are_rows(ps)) OTU <- t(OTU)
    relab0 <- sweep(OTU, 1, pmax(1e-12, rowSums(OTU)), "/")
    keep  <- (colMeans(relab0 > 0) >= prev_min) | (colMeans(relab0) >= relab_min)
    relab <- relab0[, keep, drop = FALSE]
  }

  ## --- 3) Taxon 라벨 맵 ------------------------------------------------------
  tax_map <- NULL
  if (!is.null(tax_table(ps)) && bacteriaSet) {
    tt <- as.data.frame(tax_table(ps), stringsAsFactors = FALSE)
    prio <- intersect(c("Species","Genus","Family","Order","Class","Phylum","Kingdom"),
                      colnames(tt))
    make_tax_label <- function(row) {
      g <- row[["Genus"]]; s <- row[["Species"]]
      if (!is.null(g) && nzchar(g) && !is.na(g) && !is.null(s) && nzchar(s) && !is.na(s)) return(paste(g, s))
      if (!is.null(g) && nzchar(g) && !is.na(g)) return(g)
      if (!is.null(s) && nzchar(s) && !is.na(s)) return(s)
      for (nm in prio) {
        v <- row[[nm]]
        if (!is.null(v) && nzchar(v) && !is.na(v)) return(paste0(v, " sp."))
      }
      return(rownames(row))
    }
    tt$TaxLabel <- apply(tt, 1, make_tax_label)
    tax_map <- setNames(tt$TaxLabel, rownames(tt))
  }

  ## --- 4) 디자인 매트릭스 ----------------------------------------------------
  meta <- data.frame(sample_data(ps), check.names = FALSE, stringsAsFactors = FALSE)
  X_micro <- if (bacteriaSet) relab else NULL

  X_cli <- NULL
  if (length(clinical_vari)) {
    miss_cli <- setdiff(clinical_vari, colnames(meta))
    if (length(miss_cli)) warning(sprintf("[Warn] clinical vars not found: %s", paste(miss_cli, collapse=", ")))
    clinical_vari <- setdiff(clinical_vari, miss_cli)
    if (length(clinical_vari)) {
      X_cli <- meta[, clinical_vari, drop = FALSE]
      for (j in colnames(X_cli)) if (is.factor(X_cli[[j]])) X_cli[[j]] <- as.character(X_cli[[j]])
      char_cols <- names(Filter(is.character, X_cli))
      if (length(char_cols)) {
        dummies <- do.call(cbind, lapply(char_cols, function(cn){
          mm <- model.matrix(~.-1, data = data.frame(val = factor(X_cli[[cn]])))
          colnames(mm) <- paste(cn, colnames(mm), sep="__"); mm
        }))
        X_cli <- cbind(X_cli[, setdiff(colnames(X_cli), char_cols), drop=FALSE], dummies)
      }
      X_cli <- data.matrix(X_cli)
      for (j in seq_len(ncol(X_cli))) { v <- X_cli[, j]; if (anyNA(v)) v[is.na(v)] <- median(v, na.rm=TRUE); X_cli[, j] <- v }
      X_cli <- scale(X_cli)
    }
  }
  if (!bacteriaSet && is.null(X_cli)) stop("[Error] bacteriaSet=FALSE 인데 clinical_vari가 비어 있습니다.")
  X <- if (!is.null(X_cli) && !is.null(X_micro)) cbind(X_micro, X_cli) else (X_cli %||% X_micro)
  X <- as.matrix(X); storage.mode(X) <- "double"

  yfac <- meta[[outcome]]                 # factor(levels=orders)
  gid  <- meta[[StudyID_col]]
  stopifnot(nrow(X) == length(yfac), length(gid) == length(yfac))
  ybin <- as.integer(yfac == positive_class)

  ## --- 5) 폴드 ---------------------------------------------------------------
  group_mode <- !all(meta[[StudyID_col]] == rownames(meta))
  folds_main <- if (group_mode) make_group_strat_folds(gid, ybin, n_folds, seed) else make_strat_folds(ybin, n_folds, seed)

  ## --- 6) 탐색공간 -----------------------------------------------------------
  sample_grid_rf <- function(n, p){
    data.frame(
      mtry            = pmax(1, round(runif(n, 0.1, 0.5) * p)),
      min.node.size   = sample(c(1,2,3,5,10), n, replace=TRUE),
      sample.fraction = runif(n, 0.6, 0.95),
      stringsAsFactors = FALSE
    )
  }
  sample_grid_xgb <- function(n){
    data.frame(
      eta               = runif(n, 0.03, 0.3),
      max_depth         = sample(2:6, n, replace=TRUE),
      min_child_weight  = sample(c(1,2,3,5), n, replace=TRUE),
      subsample         = runif(n, 0.6, 0.95),
      colsample_bytree  = runif(n, 0.5, 0.95),
      gamma             = runif(n, 0, 3),
      lambda            = 10^runif(n, -2, 1),
      alpha             = 10^runif(n, -3, 1),
      early_stopping_rounds = sample(c(30,40,50,75), n, replace=TRUE),
      stringsAsFactors = FALSE
    )
  }

  ## ====================== RandomForest =======================================
  if (method == "randomforest") {
    taby <- table(yfac)
    class.weights <- as.numeric(max(taby) / taby); names(class.weights) <- names(taby)

    eval_cv_auc_rf <- function(X, yfac, groups=NULL, params, n_folds=5, seed=123){
      folds <- if (!is.null(groups)) make_group_strat_folds(groups, as.integer(yfac==positive_class), n_folds, seed)
      else make_strat_folds(as.integer(yfac==positive_class), n_folds, seed)
      oof <- rep(NA_real_, length(yfac))
      for (k in seq_along(folds)) {
        te <- folds[[k]]; tr <- setdiff(seq_along(yfac), te)
        df_tr <- data.frame(X[tr,,drop=FALSE]); df_tr[[outcome]] <- yfac[tr]
        rf <- ranger::ranger(
          formula = as.formula(paste(outcome, "~ .")),
          data = df_tr,
          num.trees = num.trees,
          mtry = params$mtry,
          min.node.size = params$min.node.size,
          sample.fraction = params$sample.fraction,
          probability = TRUE,
          classification = TRUE,
          importance = "permutation",
          class.weights = class.weights,
          seed = seed
        )
        oof[te] <- predict(rf, data.frame(X[te,,drop=FALSE]))$predictions[, positive_class]
      }
      auc <- get_auc(yfac, oof, orders)
      list(auc=auc, pred=oof, folds=folds)
    }

    grid <- sample_grid_rf(n_candidates, ncol(X))
    cat(sprintf("[TUNE] random search (RF): %d candidates\n", nrow(grid)))
    best_auc <- -Inf; best_ix <- 1; evals <- vector("list", nrow(grid))
    for (i in seq_len(nrow(grid))) {
      cv <- eval_cv_auc_rf(X, yfac, if (group_mode) gid else NULL, grid[i,], n_folds, seed)
      evals[[i]] <- cbind(grid[i,], AUC=cv$auc)
      if (!is.na(cv$auc) && cv$auc > best_auc) { best_auc <- cv$auc; best_ix <- i }
      cat(sprintf("  - cand %2d/%2d: AUC=%.3f (mtry=%d, min.node=%d, frac=%.2f)\n",
                  i, nrow(grid), cv$auc, grid$mtry[i], grid$min.node.size[i], grid$sample.fraction[i]))
    }
    eval_df <- do.call(rbind, evals); eval_df <- eval_df[order(-eval_df$AUC), ]
    write.csv(eval_df, file.path(root, "random_search_results.csv"), row.names = FALSE)
    best_param <- grid[best_ix, ]
    cat(sprintf("[Best CV] AUC=%.3f (RF)\n", best_auc))

    # OOF 예측 + PR
    cv <- eval_cv_auc_rf(X, yfac, if (group_mode) gid else NULL, best_param, n_folds, seed)
    oof <- cv$pred
    AUC_oof <- get_auc(yfac, oof, orders)
    pr <- PRROC::pr.curve(scores.class0 = oof[yfac==positive_class],
                          scores.class1 = oof[yfac!=positive_class], curve=TRUE)
    AUPRC_oof <- pr$auc.integral

    # 최종 모델 (전체)
    df_all <- data.frame(X); df_all[[outcome]] <- yfac
    final_rf <- ranger::ranger(
      formula = as.formula(paste(outcome, "~ .")),
      data = df_all,
      num.trees = num.trees,
      mtry = best_param$mtry,
      min.node.size = best_param$min.node.size,
      sample.fraction = best_param$sample.fraction,
      probability = TRUE,
      classification = TRUE,
      importance = "permutation",
      class.weights = class.weights,
      oob.error = TRUE,
      seed = seed
    )
    saveRDS(final_rf, file.path(root, "rf_final_model.rds"))
    saveRDS(list(outcome=outcome, StudyID_col=StudyID_col, positive_class=positive_class,
                 best_param=best_param, num.trees=num.trees, folds=cv$folds,
                 seed=seed, levels=orders),
            file.path(root, "rf_meta.rds"))

    write.csv(data.frame(SampleID=rownames(X), StudyID=gid, outcome=yfac, set="OOF", pred=oof),
              file.path(root, "predictions.csv"), row.names = FALSE)

    # 중요도 + 방향성 (OOF 기반) — 안전 매핑
    dir_val <- sapply(seq_len(ncol(X)), function(j) safe_spearman(X[,j], oof))
    dir_lab <- ifelse(dir_val > 0, "Positive", ifelse(dir_val < 0, "Negative", "Neutral"))
    imp <- as.data.frame(ranger::importance(final_rf))
    colnames(imp) <- "PermImportance"; imp$Feature <- rownames(imp)
    imp$Direction <- .map_direction_safe(imp$Feature, setNames(dir_lab, colnames(X)))

    base_feat <- colnames(relab %||% matrix(nrow=nrow(X), ncol=0))
    lab_vec <- if (!is.null(tax_map) && length(base_feat)) {
      ifelse(base_feat %in% names(tax_map), tax_map[base_feat], base_feat)
    } else base_feat
    feat2lab <- setNames(as.character(lab_vec), base_feat)
    imp$Taxon <- ifelse(imp$Feature %in% names(feat2lab), feat2lab[imp$Feature], imp$Feature)
    imp <- imp[order(-imp$PermImportance), ]
    write.csv(imp, file.path(root, "importance_feature.csv"), row.names = FALSE)

    importance_taxon <- imp %>%
      group_by(Taxon) %>%
      summarise(
        SumImportance  = sum(PermImportance, na.rm=TRUE),
        MeanImportance = mean(PermImportance, na.rm=TRUE),
        Direction = {
          s <- sum(ifelse(Direction=="Positive", 1, ifelse(Direction=="Negative",-1,0)) * PermImportance, na.rm=TRUE)
          ifelse(s > 0, "Positive", ifelse(s < 0, "Negative", "Neutral"))
        },
        .groups="drop"
      ) %>% arrange(desc(SumImportance))
    write.csv(importance_taxon, file.path(root, "importance_taxon.csv"), row.names = FALSE)

    png(file.path(root, "ROC_OOF.png"), width=900, height=900, res=130)
    plot.roc(pROC::roc(yfac, oof, levels=orders, direction="<", quiet=TRUE),
             main=sprintf("OOF ROC (AUC=%.3f)", AUC_oof), col="#1f77b4", lwd=2); abline(0,1,lty=2,col="grey60"); dev.off()
    png(file.path(root, "PR_OOF.png"), width=900, height=900, res=130)
    plot(pr, main=sprintf("OOF PR (AUPRC=%.3f)", AUPRC_oof)); dev.off()

    if (!testSet) {
      return(invisible(list(mode="CV-only", AUC=AUC_oof, AUPRC=AUPRC_oof,
                            best_param=best_param, outdir=root)))
    }

    ## --- Holdout 75/25 (Train=OOB, Test=Holdout) -----------------------------
    set.seed(seed)
    if (group_mode) {
      df_g <- data.frame(gid=gid, y=ybin)
      rep_lab <- aggregate(y ~ gid, df_g, function(z) round(mean(z)))
      labs <- rep_lab$gid
      g0 <- rep_lab$gid[rep_lab$y==0]; g1 <- rep_lab$gid[rep_lab$y==1]
      n0_tr <- floor(length(g0)*0.75); n1_tr <- floor(length(g1)*0.75)
      tr_g <- c(sample(g0, n0_tr), sample(g1, n1_tr))
      te_g <- setdiff(labs, tr_g)
      tr_idx <- which(gid %in% tr_g); te_idx <- which(gid %in% te_g)
    } else {
      pos_idx <- which(ybin==1); neg_idx <- which(ybin==0)
      tr_pos <- sample(pos_idx, floor(0.75*length(pos_idx)))
      tr_neg <- sample(neg_idx, floor(0.75*length(neg_idx)))
      tr_idx <- sort(c(tr_pos, tr_neg)); te_idx <- setdiff(seq_along(ybin), tr_idx)
    }

    Xtr <- X[tr_idx,,drop=FALSE]; ytr <- yfac[tr_idx]; Xte <- X[te_idx,,drop=FALSE]; yte <- yfac[te_idx]

    df_tr <- data.frame(Xtr); df_tr[[outcome]] <- ytr
    bst <- ranger::ranger(
      formula = as.formula(paste(outcome, "~ .")),
      data = df_tr,
      num.trees = num.trees,
      mtry = best_param$mtry,
      min.node.size = best_param$min.node.size,
      sample.fraction = best_param$sample.fraction,
      probability = TRUE,
      classification = TRUE,
      importance = "permutation",
      class.weights = class.weights,
      oob.error = TRUE,     # Train 성능은 OOB로
      seed = seed
    )
    ph_tr <- bst$predictions[, positive_class]                                   # OOB
    ph_te <- predict(bst, data.frame(Xte))$predictions[, positive_class]         # Holdout
    auc_tr <- get_auc(ytr, ph_tr, orders); auc_te <- get_auc(yte, ph_te, orders)
    pr_tr  <- PRROC::pr.curve(scores.class0 = ph_tr[ytr==positive_class],
                              scores.class1 = ph_tr[ytr!=positive_class], curve=TRUE)
    pr_te  <- PRROC::pr.curve(scores.class0 = ph_te[yte==positive_class],
                              scores.class1 = ph_te[yte!=positive_class], curve=TRUE)

    write.csv(rbind(
      data.frame(SampleID=rownames(Xtr), StudyID=gid[tr_idx], outcome=ytr, set="Train", pred=ph_tr),
      data.frame(SampleID=rownames(Xte), StudyID=gid[te_idx], outcome=yte, set="Test",  pred=ph_te)
    ), file.path(root, "predictions.csv"), row.names = FALSE)

    # 중요도 + 방향성 (테스트 기준 가능하면 test, 아니면 train OOB)
    dir_base <- if (length(ph_te) >= 5) ph_te else ph_tr
    X_base   <- if (length(ph_te) >= 5) Xte   else Xtr
    dir_val  <- sapply(seq_len(ncol(X_base)), function(j) safe_spearman(X_base[,j], dir_base))
    dir_lab  <- ifelse(dir_val > 0, "Positive", ifelse(dir_val < 0, "Negative", "Neutral"))
    imp2 <- as.data.frame(ranger::importance(bst))
    colnames(imp2) <- "PermImportance"; imp2$Feature <- rownames(imp2)
    imp2$Direction <- .map_direction_safe(imp2$Feature, setNames(dir_lab, colnames(X_base)))

    base_feat <- colnames(relab %||% matrix(nrow=nrow(X), ncol=0))
    lab_vec <- if (!is.null(tax_map) && length(base_feat)) {
      ifelse(base_feat %in% names(tax_map), tax_map[base_feat], base_feat)
    } else base_feat
    feat2lab <- setNames(as.character(lab_vec), base_feat)
    imp2$Taxon <- ifelse(imp2$Feature %in% names(feat2lab), feat2lab[imp2$Feature], imp2$Feature)
    imp2 <- imp2[order(-imp2$PermImportance), ]
    write.csv(imp2, file.path(root, "importance_feature.csv"), row.names = FALSE)

    importance_taxon2 <- imp2 %>%
      group_by(Taxon) %>%
      summarise(
        SumImportance  = sum(PermImportance, na.rm=TRUE),
        MeanImportance = mean(PermImportance, na.rm=TRUE),
        Direction = {
          s <- sum(ifelse(Direction=="Positive", 1, ifelse(Direction=="Negative",-1,0)) * PermImportance, na.rm=TRUE)
          ifelse(s > 0, "Positive", ifelse(s < 0, "Negative", "Neutral"))
        },
        .groups="drop"
      ) %>% arrange(desc(SumImportance))
    write.csv(importance_taxon2, file.path(root, "importance_taxon.csv"), row.names = FALSE)

    png(file.path(root, "ROC_Train.png"), width=900, height=900, res=130)
    plot.roc(pROC::roc(ytr, ph_tr, levels=orders, direction="<", quiet=TRUE),
             main=sprintf("Train ROC (AUC=%.3f)", auc_tr), col="#1f77b4", lwd=2); abline(0,1,lty=2,col="grey60"); dev.off()
    png(file.path(root, "ROC_Test.png"), width=900, height=900, res=130)
    plot.roc(pROC::roc(yte, ph_te, levels=orders, direction="<", quiet=TRUE),
             main=sprintf("Test ROC (AUC=%.3f)", auc_te), col="#d62728", lwd=2); abline(0,1,lty=2,col="grey60"); dev.off()
    png(file.path(root, "PR_Train.png"), width=900, height=900, res=130)
    plot(pr_tr, main=sprintf("Train PR (AUPRC=%.3f)", pr_tr$auc.integral)); dev.off()
    png(file.path(root, "PR_Test.png"), width=900, height=900, res=130)
    plot(pr_te, main=sprintf("Test PR (AUPRC=%.3f)", pr_te$auc.integral)); dev.off()

    return(invisible(list(mode="Holdout+CV", AUC_train=auc_tr, AUC_test=auc_te,
                          AUPRC_train=pr_tr$auc.integral, AUPRC_test=pr_te$auc.integral,
                          best_param=best_param, outdir=root)))
  }

  ## ======================== XGBoost ==========================================
  dmat_all <- xgb.DMatrix(data = X, label = ybin, missing = NA)
  scale_pos_weight <- if (sum(ybin==1)>0) sum(ybin==0)/sum(ybin==1) else 1

  eval_candidate <- function(par, folds){
    param <- list(
      objective = "binary:logistic",
      eval_metric = "auc",
      tree_method = "hist",
      scale_pos_weight = scale_pos_weight,
      eta = par$eta,
      max_depth = par$max_depth,
      min_child_weight = par$min_child_weight,
      subsample = par$subsample,
      colsample_bytree = par$colsample_bytree,
      gamma = par$gamma,
      lambda = par$lambda,
      alpha  = par$alpha
    )
    cv <- xgb.cv(params = param, data = dmat_all, nrounds = num.trees, folds = folds,
                 early_stopping_rounds = par$early_stopping_rounds, verbose = 0, maximize = TRUE)
    list(auc = cv$evaluation_log$test_auc_mean[cv$best_iteration],
         best_iter = cv$best_iteration,
         param = param,
         esr = par$early_stopping_rounds)
  }

  set.seed(seed)
  grid <- sample_grid_xgb(n_candidates)
  cat(sprintf("[TUNE] random search (XGB): %d candidates\n", nrow(grid)))
  best <- list(auc=-Inf)
  res_list <- vector("list", nrow(grid))
  for (i in seq_len(nrow(grid))) {
    r <- eval_candidate(grid[i,], folds_main)
    res_list[[i]] <- data.frame(idx=i, AUC=r$auc, best_iter=r$best_iter, grid[i,])
    cat(sprintf("  - cand %2d/%2d: AUC=%.3f | iter=%d (md=%d, eta=%.2f, mcw=%d, sub=%.2f, col=%.2f)\n",
                i, nrow(grid), r$auc, r$best_iter, grid$max_depth[i], grid$eta[i],
                grid$min_child_weight[i], grid$subsample[i], grid$colsample_bytree[i]))
    if (!is.na(r$auc) && r$auc > best$auc) best <- r
  }
  res_df <- do.call(rbind, res_list); res_df <- res_df[order(-res_df$AUC), ]
  write.csv(res_df, file.path(root, "random_search_results.csv"), row.names = FALSE)
  cat(sprintf("[Best CV] AUC=%.3f | best_iter=%d (XGB)\n", best$auc, best$best_iter))



  ## OOF 예측
  oof <- rep(NA_real_, length(ybin))
  best_param <- best$param; best_iter <- best$best_iter
  for (k in seq_along(folds_main)) {
    te <- folds_main[[k]]; tr <- setdiff(seq_along(ybin), te)
    dtr <- xgb.DMatrix(X[tr,,drop=FALSE], label=ybin[tr], missing=NA)
    dte <- xgb.DMatrix(X[te,,drop=FALSE], label=ybin[te], missing=NA)
    bst <- xgb.train(params = best_param, data = dtr, nrounds = best_iter, verbose = 0)
    oof[te] <- predict(bst, dte)
  }
  AUC_oof <- get_auc(yfac, oof, orders)
  pr <- PRROC::pr.curve(scores.class0 = oof[yfac==positive_class],
                        scores.class1 = oof[yfac!=positive_class], curve=TRUE)
  AUPRC_oof <- pr$auc.integral

  full_bst <- xgb.train(params = best_param, data = dmat_all, nrounds = best_iter, verbose = 0)
  saveRDS(full_bst, file.path(root, "xgb_final_model.rds"))
  saveRDS(list(outcome=outcome, StudyID_col=StudyID_col, positive_class=positive_class,
               best_param=best_param, best_iter=best_iter, seed=seed, folds=folds_main,
               levels=orders),
          file.path(root, "xgb_meta.rds"))

  write.csv(data.frame(SampleID=rownames(X), StudyID=gid, outcome=yfac, set="OOF", pred=oof),
            file.path(root, "predictions.csv"), row.names = FALSE)




  CV_best_AUC  <- best$auc
  CV_best_iter <- best$best_iter
  log_cv_oof_summary(outdir = root, model_name = "xgb",
                     CV_best_AUC = CV_best_AUC, CV_best_iter = CV_best_iter,
                     oof_prob = oof, yfac = yfac,
                     positive_class = positive_class, orders = orders)



  ## 중요도 + 방향성 (SHAP 우선, 실패시 Spearman) — 안전 매핑
  shap_ok <- TRUE
  shap_mat <- try(predict(full_bst, dmat_all, predcontrib=TRUE), silent=TRUE)
  if (inherits(shap_mat, "try-error")) shap_ok <- FALSE
  feat_names <- colnames(X)
  if (shap_ok) {
    shap_mat <- matrix(shap_mat, ncol = ncol(X)+1, byrow = TRUE)
    colnames(shap_mat) <- c(feat_names, "BIAS")
    MeanAbsShap <- colMeans(abs(shap_mat[, feat_names, drop=FALSE]))
    MeanShap    <- colMeans(shap_mat[, feat_names, drop=FALSE])
    Direction <- ifelse(MeanShap > 0, "Positive", ifelse(MeanShap < 0, "Negative", "Neutral"))
    imp_tbl <- data.frame(Feature=feat_names, PermImportance=MeanAbsShap, Direction=Direction, check.names=FALSE)
  } else {
    dir_val <- sapply(seq_len(ncol(X)), function(j) safe_spearman(X[,j], oof))
    Direction <- ifelse(dir_val > 0, "Positive", ifelse(dir_val < 0, "Negative", "Neutral"))
    names(Direction) <- colnames(X)
    imp_xgb <- try(xgb.importance(model = full_bst), silent=TRUE)
    if (inherits(imp_xgb, "try-error") || is.null(imp_xgb)) {
      imp_tbl <- data.frame(Feature=feat_names, PermImportance=NA_real_, Direction=Direction, check.names=FALSE)
    } else {
      imp_tbl <- imp_xgb[, c("Feature","Gain")]; colnames(imp_tbl) <- c("Feature","PermImportance")
      imp_tbl$Direction <- .map_direction_safe(imp_tbl$Feature, Direction)
    }
  }

  base_feat <- colnames(relab %||% matrix(nrow=nrow(X), ncol=0))
  lab_vec <- if (!is.null(tax_map) && length(base_feat)) {
    ifelse(base_feat %in% names(tax_map), tax_map[base_feat], base_feat)
  } else base_feat
  feat2lab <- setNames(as.character(lab_vec), base_feat)
  imp_tbl$Taxon <- ifelse(imp_tbl$Feature %in% names(feat2lab), feat2lab[imp_tbl$Feature], imp_tbl$Feature)
  imp_tbl <- imp_tbl[order(-imp_tbl$PermImportance), ]
  write.csv(imp_tbl, file.path(root, "importance_feature.csv"), row.names = FALSE)

  importance_taxon <- imp_tbl %>%
    group_by(Taxon) %>%
    summarise(
      SumImportance  = sum(PermImportance, na.rm=TRUE),
      MeanImportance = mean(PermImportance, na.rm=TRUE),
      Direction = {
        s <- sum(ifelse(Direction=="Positive", 1, ifelse(Direction=="Negative",-1, 0)) * PermImportance, na.rm=TRUE)
        ifelse(s > 0, "Positive", ifelse(s < 0, "Negative", "Neutral"))
      },
      .groups="drop"
    ) %>% arrange(desc(SumImportance))
  write.csv(importance_taxon, file.path(root, "importance_taxon.csv"), row.names = FALSE)

  png(file.path(root, "ROC_OOF.png"), width=900, height=900, res=130)
  plot.roc(pROC::roc(yfac, oof, levels=orders, direction="<", quiet=TRUE),
           main=sprintf("OOF ROC (AUC=%.3f)", AUC_oof), col="#1f77b4", lwd=2); abline(0,1,lty=2,col="grey60"); dev.off()
  png(file.path(root, "PR_OOF.png"), width=900, height=900, res=130)
  plot(pr, main=sprintf("OOF PR (AUPRC=%.3f)", AUPRC_oof)); dev.off()

  if (!testSet) {
    return(invisible(list(mode="CV-only", AUC=AUC_oof, AUPRC=AUPRC_oof,
                          best_param=best_param, outdir=root)))
  }

  ## --- Holdout (Train OOF, Test Holdout) ------------------------------------
  set.seed(seed)
  if (group_mode) {
    df_g <- data.frame(gid=gid, y=ybin)
    rep_lab <- aggregate(y ~ gid, df_g, function(z) round(mean(z)))
    labs <- rep_lab$gid
    g0 <- rep_lab$gid[rep_lab$y==0]; g1 <- rep_lab$gid[rep_lab$y==1]
    n0_tr <- floor(length(g0)*0.75); n1_tr <- floor(length(g1)*0.75)
    tr_g <- c(sample(g0, n0_tr), sample(g1, n1_tr))
    te_g <- setdiff(labs, tr_g)
    tr_idx <- which(gid %in% tr_g); te_idx <- which(gid %in% te_g)
  } else {
    pos_idx <- which(ybin==1); neg_idx <- which(ybin==0)
    tr_pos <- sample(pos_idx, floor(0.75*length(pos_idx)))
    tr_neg <- sample(neg_idx, floor(0.75*length(neg_idx)))
    tr_idx <- sort(c(tr_pos, tr_neg)); te_idx <- setdiff(seq_along(ybin), tr_idx)
  }

  Xtr <- X[tr_idx,,drop=FALSE]; ytr <- ybin[tr_idx]; ytr_fac <- yfac[tr_idx]
  Xte <- X[te_idx,,drop=FALSE]; yte <- ybin[te_idx]; yte_fac <- yfac[te_idx]

  dtr_all <- xgb.DMatrix(Xtr, label=ytr, missing=NA)
  folds_tr <- if (group_mode) make_group_strat_folds(gid[tr_idx], ytr, n_folds, seed) else make_strat_folds(ytr, n_folds, seed)
  grid_in <- sample_grid_xgb(n_candidates)
  best_in <- list(auc=-Inf)
  for (i in seq_len(nrow(grid_in))) {
    par <- grid_in[i, ]
    param <- list(
      objective="binary:logistic", eval_metric="auc", tree_method="hist",
      scale_pos_weight = if (sum(ytr==1)>0) sum(ytr==0)/sum(ytr==1) else 1,
      eta=par$eta, max_depth=par$max_depth, min_child_weight=par$min_child_weight,
      subsample=par$subsample, colsample_bytree=par$colsample_bytree,
      gamma=par$gamma, lambda=par$lambda, alpha=par$alpha
    )
    cv <- xgb.cv(params=param, data=dtr_all, nrounds=num.trees, folds=folds_tr,
                 early_stopping_rounds=par$early_stopping_rounds, verbose=0, maximize=TRUE)
    auc_i <- cv$evaluation_log$test_auc_mean[cv$best_iteration]
    if (!is.na(auc_i) && auc_i > best_in$auc)
      best_in <- list(auc=auc_i, best_iter=cv$best_iteration, param=param, esr=par$early_stopping_rounds)
  }
  cat(sprintf("[InnerCV] best AUC=%.3f | best_iter=%d (XGB)\n", best_in$auc, best_in$best_iter))

  dtr <- xgb.DMatrix(Xtr, label=ytr, missing=NA)
  dte <- xgb.DMatrix(Xte, label=yte, missing=NA)

  # Train 성능: 내부 폴드 OOF로 계산 (ROC=1 방지)
  ph_tr <- rep(NA_real_, length(ytr))
  for (k in seq_along(folds_tr)) {
    te_k <- folds_tr[[k]]
    tr_k <- setdiff(seq_along(ytr), te_k)
    dtr_k <- xgb.DMatrix(Xtr[tr_k,,drop=FALSE], label=ytr[tr_k], missing=NA)
    dte_k <- xgb.DMatrix(Xtr[te_k,,drop=FALSE], label=ytr[te_k], missing=NA)
    bst_k <- xgb.train(params=best_in$param, data=dtr_k, nrounds=best_in$best_iter, verbose=0)
    ph_tr[te_k] <- predict(bst_k, dte_k)
  }
  # Test: train 전체로 학습 후 holdout 예측
  bst_tr <- xgb.train(params=best_in$param, data=dtr, nrounds=best_in$best_iter, verbose=0)
  ph_te  <- predict(bst_tr, dte)

  auc_tr <- get_auc(ytr_fac, ph_tr, orders); auc_te <- get_auc(yte_fac, ph_te, orders)
  pr_tr  <- PRROC::pr.curve(scores.class0 = ph_tr[ytr_fac==positive_class],
                            scores.class1 = ph_tr[ytr_fac!=positive_class], curve=TRUE)
  pr_te  <- PRROC::pr.curve(scores.class0 = ph_te[yte_fac==positive_class],
                            scores.class1 = ph_te[yte_fac!=positive_class], curve=TRUE)

  write.csv(rbind(
    data.frame(SampleID=rownames(Xtr), StudyID=gid[tr_idx], outcome=ytr_fac, set="Train", pred=ph_tr),
    data.frame(SampleID=rownames(Xte), StudyID=gid[te_idx], outcome=yte_fac, set="Test",  pred=ph_te)
  ), file.path(root, "predictions.csv"), row.names = FALSE)

  # 중요도 + 방향성 (Train OOF 기반 방향, 안전 매핑)
  dir_val <- sapply(seq_len(ncol(Xtr)), function(j) safe_spearman(Xtr[,j], ph_tr))
  Direction <- ifelse(dir_val > 0, "Positive", ifelse(dir_val < 0, "Negative", "Neutral"))
  names(Direction) <- colnames(Xtr)
  imp_x <- try(xgb.importance(model = bst_tr), silent=TRUE)
  if (inherits(imp_x, "try-error") || is.null(imp_x)) {
    imp_hold <- data.frame(Feature=colnames(Xtr), PermImportance=NA_real_, Direction=Direction, check.names=FALSE)
  } else {
    imp_hold <- imp_x[, c("Feature","Gain")]; colnames(imp_hold) <- c("Feature","PermImportance")
    imp_hold$Direction <- .map_direction_safe(imp_hold$Feature, Direction)
  }
  base_feat <- colnames(relab %||% matrix(nrow=nrow(X), ncol=0))
  lab_vec <- if (!is.null(tax_map) && length(base_feat)) {
    ifelse(base_feat %in% names(tax_map), tax_map[base_feat], base_feat)
  } else base_feat
  feat2lab <- setNames(as.character(lab_vec), base_feat)
  imp_hold$Taxon <- ifelse(imp_hold$Feature %in% names(feat2lab), feat2lab[imp_hold$Feature], imp_hold$Feature)
  imp_hold <- imp_hold[order(-imp_hold$PermImportance), ]
  write.csv(imp_hold, file.path(root, "importance_feature.csv"), row.names = FALSE)

  importance_taxon_h <- imp_hold %>%
    group_by(Taxon) %>%
    summarise(
      SumImportance  = sum(PermImportance, na.rm=TRUE),
      MeanImportance = mean(PermImportance, na.rm=TRUE),
      Direction = {
        s <- sum(ifelse(Direction=="Positive", 1, ifelse(Direction=="Negative",-1, 0)) * PermImportance, na.rm=TRUE)
        ifelse(s > 0, "Positive", ifelse(s < 0, "Negative", "Neutral"))
      },
      .groups="drop"
    ) %>% arrange(desc(SumImportance))
  write.csv(importance_taxon_h, file.path(root, "importance_taxon.csv"), row.names = FALSE)

  png(file.path(root, "ROC_Train.png"), width=900, height=900, res=130)
  plot.roc(pROC::roc(ytr_fac, ph_tr, levels=orders, direction="<", quiet=TRUE),
           main=sprintf("Train ROC (AUC=%.3f)", auc_tr), lwd=2); abline(0,1,lty=2); dev.off()
  png(file.path(root, "ROC_Test.png"), width=900, height=900, res=130)
  plot.roc(pROC::roc(yte_fac, ph_te, levels=orders, direction="<", quiet=TRUE),
           main=sprintf("Test ROC (AUC=%.3f)", auc_te), lwd=2); abline(0,1,lty=2); dev.off()
  png(file.path(root, "PR_Train.png"), width=900, height=900, res=130)
  plot(pr_tr, main=sprintf("Train PR (AUPRC=%.3f)", pr_tr$auc.integral)); dev.off()
  png(file.path(root, "PR_Test.png"), width=900, height=900, res=130)
  plot(pr_te, main=sprintf("Test PR (AUPRC=%.3f)", pr_te$auc.integral)); dev.off()

  return(invisible(list(mode="Holdout+CV", AUC_train=auc_tr, AUC_test=auc_te,
                        AUPRC_train=pr_tr$auc.integral, AUPRC_test=pr_te$auc.integral,
                        best_param=best_in$param, outdir=root)))
}

`%||%` <- function(a, b) if (!is.null(a)) a else b
