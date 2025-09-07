#' Go_xgboost
#'
#' StudyID-aware XGBoost classification for microbiome (+ optional clinical covariates)
#' with stratified/grouped cross-validation, random hyperparameter search, SHAP-based
#' (or gain+Spearman) feature direction, and rich outputs.
#'
#' @description
#' Builds a binary XGBoost model from a \code{phyloseq} object’s relative abundances
#' (optionally aggregated to a taxonomic rank) and optional clinical variables.
#' The function:
#' \itemize{
#'   \item optionally aggregates features by \code{taxrank} (e.g., \code{"Genus"}),
#'   \item converts to relative abundances and filters by prevalence/mean abundance,
#'   \item one-hot encodes character clinical covariates, imputes NA with medians, and scales,
#'   \item fixes outcome level order (\code{orders}; positive class is \code{orders[2]}),
#'   \item does random hyperparameter search via \code{xgb.cv} with K folds that preserve
#'         subject groups when \code{StudyID_col} varies across rows,
#'   \item trains the final model and saves OOF (or Train/Test) predictions, feature
#'         importance (with direction), aggregated importance by taxon, and ROC/PR plots.
#' }
#'
#' @param psIN A \code{phyloseq} object.
#' @param project Character; project prefix used to create an output directory
#'   \verb{<project_YYMMDD>/XGBoost}.
#' @param outcome Character; name of the binary outcome column in sample metadata.
#'   Levels must include \code{orders}.
#' @param testSet Logical; if \code{FALSE} (default) runs CV-only. If \code{TRUE},
#'   creates a 75/25 group-aware holdout split plus inner CV on the train set.
#' @param clinical_vari Character vector of metadata columns to add as covariates
#'   (optional). Character/factor covariates are one-hot encoded.
#' @param StudyID_col Character; subject/group id column. If missing, rownames are used. Default \code{"StudyID"}.
#' @param taxrank Either \code{"ASV"} to keep features as is, or a taxonomic rank
#'   present in \code{tax_table(psIN)} (e.g., \code{"Genus"}). Default \code{"ASV"}.
#' @param prev_min Numeric in \([0,1]\); minimum prevalence threshold on relative
#'   abundance presence. Default \code{0.01}.
#' @param relab_min Numeric; minimum mean relative abundance threshold. Default \code{1e-4}.
#' @param n_folds Integer; number of CV folds. Default \code{5}.
#' @param seed Integer; random seed. Default \code{123}.
#' @param n_candidates Integer; number of random hyperparameter candidates. Default \code{40}.
#' @param num.trees Integer; maximum boosting rounds used for CV/early stopping
#'   (best iteration is selected). Default \code{1000}.
#' @param orders Character vector of length 2 giving the outcome level order,
#'   e.g., \code{c("Control","Case")}. The positive class is \code{orders[2]}.
#'
#' @details
#' \strong{Grouping and stratification:} If \code{StudyID_col} does not equal rownames
#' for all samples, folds are created to keep samples from the same StudyID together,
#' while preserving class balance at the group level when possible. Otherwise,
#' standard stratified folds are used. Class imbalance is addressed via
#' \code{scale_pos_weight = neg/pos}.
#'
#' \strong{Feature direction:} If SHAP values are available (using \code{predcontrib=TRUE}),
#' feature “Direction” is assigned by the mean SHAP value sign; otherwise, direction
#' falls back to Spearman correlation between feature values and OOF (or Train) predictions.
#'
#' \strong{Outputs:} The function writes a dated output folder with the tuned results,
#' predictions, feature/taxon importances (including Direction), and ROC/PR plots.
#'
#' @return (Invisibly) a list with run \code{mode}, metrics (AUC/AUPRC), \code{best_param},
#'   and \code{outdir}. Side effects: writes to \verb{<project_YYMMDD>/XGBoost/}:
#' \itemize{
#'   \item \code{random_search_results.csv}
#'   \item \code{predictions.csv} (OOF for CV-only; Train/Test for holdout mode)
#'   \item \code{xgb_final_model.rds} (or train model in holdout mode), \code{xgb_meta.rds}
#'   \item \code{importance_feature.csv}, \code{importance_taxon.csv}
#'   \item ROC/PR PNGs (\code{ROC_OOF.png}, \code{PR_OOF.png}, or Train/Test variants)
#' }
#'
#' @examples
#' \dontrun{
#' # Binary outcome with levels c("Control","Case") in metadata column "Status"
#' res <- Go_xgboost(
#'   psIN = ps,
#'   project = "IBD",
#'   outcome = "Status",
#'   clinical_vari = c("Age","BMI","Sex"),
#'   StudyID_col = "SubjectID",
#'   taxrank = "Genus",
#'   prev_min = 0.05,
#'   relab_min = 1e-4,
#'   n_folds = 5,
#'   n_candidates = 40,
#'   num.trees = 1000,
#'   orders = c("Control","Case")
#' )
#' res$AUC
#' }
#'
#' @importFrom phyloseq sample_data tax_table taxa_are_rows otu_table tax_glom
#' @importFrom xgboost xgb.cv xgb.train xgb.importance xgb.DMatrix
#' @importFrom pROC roc auc plot.roc
#' @importFrom PRROC pr.curve
#' @importFrom dplyr group_by summarise arrange
#' @importFrom stats model.matrix median cor
#' @importFrom utils write.csv
#' @importFrom grDevices png dev.off
#' @export



`%||%` <- function(a, b) if (!is.null(a)) a else b

Go_xgboost <- function(
    psIN,
    project,
    outcome,
    testSet = FALSE,
    clinical_vari = character(0),
    StudyID_col   = "StudyID",
    taxrank       = "ASV",
    prev_min      = 0.01,
    relab_min     = 1e-4,
    n_folds       = 5,
    seed          = 123,
    n_candidates  = 40,     # 랜덤 탐색 후보 수
    num.trees     = 1000,   # 최대 nrounds (조기중단으로 best_iter 선택)
    orders        = c("Control","Case")
){
  set.seed(seed)

  ## --- 0) 메타 준비: 레벨/StudyID 고정 (psIN에도 반영) -----------------------
  meta0 <- data.frame(sample_data(psIN), check.names = FALSE, stringsAsFactors = FALSE)
  if (!(StudyID_col %in% names(meta0)) || all(is.na(meta0[[StudyID_col]]))) {
    message(sprintf("[Info] '%s'가 없어 rownames를 사용합니다.", StudyID_col))
    meta0[[StudyID_col]] <- rownames(meta0)
  }
  stopifnot(outcome %in% names(meta0))
  yy0 <- meta0[[outcome]]; if (!is.factor(yy0)) yy0 <- factor(yy0)
  if (!all(orders %in% levels(yy0))) {
    stop(sprintf("[Error] outcome(levels=%s)에 orders(%s) 중 일부가 없습니다.",
                 paste(levels(yy0), collapse=", "), paste(orders, collapse=", ")))
  }
  yy0 <- factor(yy0, levels = orders)
  meta0[[outcome]] <- yy0
  sample_data(psIN)[[outcome]] <- yy0
  if (!(StudyID_col %in% colnames(sample_data(psIN)))) {
    sample_data(psIN)[[StudyID_col]] <- meta0[[StudyID_col]]
  }
  positive_class <- orders[2]

  ## --- 1) 출력 폴더 ----------------------------------------------------------
  today <- format(Sys.Date(), "%y%m%d")
  root  <- file.path(sprintf("%s_%s", project, today), "XGBoost")
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  message("[OutDir] ", root)

  ## --- 2) phyloseq → feature table ------------------------------------------
  stopifnot(inherits(psIN, "phyloseq"))
  ps <- psIN
  if (!identical(taxrank, "ASV")) {
    if (is.null(tax_table(ps)) || !(taxrank %in% colnames(tax_table(ps)))) {
      stop(sprintf("[Error] taxrank '%s' not found in tax_table.", taxrank))
    }
    ps <- tax_glom(ps, taxrank = taxrank, NArm = TRUE)
  }
  OTU <- as(otu_table(ps), "matrix"); if (taxa_are_rows(ps)) OTU <- t(OTU)
  relab <- sweep(OTU, 1, pmax(1e-12, rowSums(OTU)), "/")
  prev <- colMeans(relab > 0)
  keep  <- (prev >= prev_min) | (colMeans(relab) >= relab_min)
  relab <- relab[, keep, drop = FALSE]

  ## Taxon label map
  tax_map <- NULL
  if (!is.null(tax_table(ps))) {
    tt <- as.data.frame(tax_table(ps), stringsAsFactors = FALSE)
    rank_cols <- colnames(tt)
    prio <- intersect(c("Species","Genus","Family","Order","Class","Phylum","Kingdom"),
                      rank_cols)
    make_tax_label <- function(row) {
      g <- row[["Genus"]]; s <- row[["Species"]]
      if (!is.null(g) && !is.na(g) && nzchar(g) && !is.null(s) && !is.na(s) && nzchar(s)) return(paste(g, s))
      if (!is.null(g) && !is.na(g) && nzchar(g)) return(g)
      if (!is.null(s) && !is.na(s) && nzchar(s)) return(s)
      for (nm in prio) {
        v <- row[[nm]]
        if (!is.null(v) && !is.na(v) && nzchar(v)) return(paste0(v, " sp."))
      }
      return(rownames(row))
    }
    tt$TaxLabel <- apply(tt, 1, make_tax_label)
    tax_map <- setNames(tt$TaxLabel, rownames(tt))
  }

  ## X: 마이크로바이옴 + 임상
  meta <- data.frame(sample_data(ps), check.names = FALSE, stringsAsFactors = FALSE)
  X_micro <- relab
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
  X <- if (!is.null(X_cli)) cbind(X_micro, X_cli) else X_micro
  X <- as.matrix(X); storage.mode(X) <- "double"
  yfac <- meta[[outcome]]                # factor with levels==orders
  ybin <- as.integer(yfac == positive_class)  # 0/1
  gid  <- meta[[StudyID_col]]
  stopifnot(length(ybin) == nrow(X), length(gid) == nrow(X))

  ## class weights → scale_pos_weight
  pos <- sum(ybin==1); neg <- sum(ybin==0)
  scale_pos_weight <- if (pos>0) neg/pos else 1

  ## --- 3) 폴드 유틸 ----------------------------------------------------------
  make_group_strat_folds <- function(groups, y, K=5, seed=123){
    set.seed(seed)
    df <- data.frame(g=groups, y=y)
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
  make_strat_folds <- function(y, K=5, seed=123){
    set.seed(seed)
    idx0 <- which(y==0); idx1 <- which(y==1)
    k0 <- if (length(idx0) < K) max(2, length(idx0)) else K
    k1 <- if (length(idx1) < K) max(2, length(idx1)) else K
    f0 <- sample(rep(1:k0, length.out=length(idx0)))
    f1 <- sample(rep(1:k1, length.out=length(idx1)))
    ff <- integer(length(y)); ff[idx0] <- f0; ff[idx1] <- f1
    lapply(1:max(ff), function(k) which(ff == k))
  }
  group_mode <- !all(meta[[StudyID_col]] == rownames(meta))
  folds_main <- if (group_mode) make_group_strat_folds(gid, ybin, n_folds, seed) else make_strat_folds(ybin, n_folds, seed)

  dmat_all <- xgb.DMatrix(data = X, label = ybin, missing = NA)

  ## --- 4) 랜덤 탐색 공간 (XGB) ----------------------------------------------
  sample_param_grid <- function(n){
    data.frame(
      eta               = runif(n, 0.03, 0.3),
      max_depth         = sample(2:6, n, replace=TRUE),
      min_child_weight  = sample(c(1,2,3,5), n, replace=TRUE),
      subsample         = runif(n, 0.6, 0.95),
      colsample_bytree  = runif(n, 0.5, 0.95),
      gamma             = runif(n, 0, 3),
      lambda            = 10^runif(n, -2, 1),  # L2
      alpha             = 10^runif(n, -3, 1),  # L1
      early_stopping_rounds = sample(c(30, 40, 50, 75), n, replace=TRUE),
      stringsAsFactors = FALSE
    )
  }

  ## --- 5) xgb.cv 평가자 ------------------------------------------------------
  eval_candidate <- function(par){
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
    cv <- xgb.cv(
      params = param,
      data = dmat_all,
      nrounds = num.trees,
      folds = folds_main,
      early_stopping_rounds = par$early_stopping_rounds,
      verbose = 0,
      maximize = TRUE
    )
    list(auc = cv$evaluation_log$test_auc_mean[cv$best_iteration],
         best_iter = cv$best_iteration,
         param = param,
         esr = par$early_stopping_rounds)
  }

  ## --- 6) 랜덤 탐색 ----------------------------------------------------------
  grid <- sample_param_grid(n_candidates)
  best <- list(auc = -Inf)
  cat(sprintf("[TUNE] random search: %d candidates\n", nrow(grid)))
  res_list <- vector("list", nrow(grid))
  for (i in seq_len(nrow(grid))) {
    r <- eval_candidate(grid[i, ])
    res_list[[i]] <- data.frame(
      idx=i,
      AUC = r$auc,
      best_iter = r$best_iter,
      eta = grid$eta[i],
      max_depth = grid$max_depth[i],
      min_child_weight = grid$min_child_weight[i],
      subsample = grid$subsample[i],
      colsample_bytree = grid$colsample_bytree[i],
      gamma = grid$gamma[i],
      lambda = grid$lambda[i],
      alpha  = grid$alpha[i],
      early_stopping_rounds = grid$early_stopping_rounds[i]
    )
    cat(sprintf("  - cand %2d/%2d: AUC=%.3f | iter=%d (md=%d, eta=%.2f, mcw=%d, sub=%.2f, col=%.2f)\n",
                i, nrow(grid), r$auc, r$best_iter, grid$max_depth[i], grid$eta[i],
                grid$min_child_weight[i], grid$subsample[i], grid$colsample_bytree[i]))
    if (!is.na(r$auc) && r$auc > best$auc) best <- r
  }
  res_df <- do.call(rbind, res_list)
  res_df <- res_df[order(-res_df$AUC), ]
  write.csv(res_df, file.path(root, "random_search_results.csv"), row.names = FALSE)
  cat(sprintf("[Best CV] AUC=%.3f | best_iter=%d\n", best$auc, best$best_iter))

  ## --- 7) OOF 예측(베스트 파라미터로) ----------------------------------------
  oof <- rep(NA_real_, length(ybin))
  best_param <- best$param; best_iter <- best$best_iter
  for (k in seq_along(folds_main)) {
    te <- folds_main[[k]]; tr <- setdiff(seq_along(ybin), te)
    dtr <- xgb.DMatrix(X[tr,,drop=FALSE], label=ybin[tr], missing=NA)
    dte <- xgb.DMatrix(X[te,,drop=FALSE], label=ybin[te], missing=NA)
    bst <- xgb.train(params = best_param, data = dtr, nrounds = best_iter, verbose = 0)
    oof[te] <- predict(bst, dte)
  }

  # ROC (그대로)
  roc_obj <- pROC::roc(response = yfac, predictor = oof,
                       levels = orders, direction = "<", quiet = TRUE)
  AUC_oof <- as.numeric(pROC::auc(roc_obj))

  # PR (scores.class0=negative, scores.class1=positive)
  pr <- PRROC::pr.curve(
    scores.class0 = oof[yfac != positive_class],
    scores.class1 = oof[yfac == positive_class],
    curve = TRUE
  )
  AUPRC_oof <- pr$auc.integral

  ## --- 8) 최종 모델(전체 학습) + 중요도/방향성 -------------------------------
  full_bst <- xgb.train(params = best_param, data = dmat_all, nrounds = best_iter, verbose = 0)
  saveRDS(full_bst, file.path(root, "xgb_final_model.rds"))

  meta_obj <- list(
    outcome = outcome,
    StudyID_col = StudyID_col,
    positive_class = positive_class,
    best_param = best_param,
    best_iter = best_iter,
    seed = seed,
    folds = folds_main,
    levels = orders
  )
  saveRDS(meta_obj, file.path(root, "xgb_meta.rds"))

  # 예측 저장 (OOF)
  pred_df <- data.frame(
    SampleID = rownames(X),
    StudyID  = gid,
    outcome  = yfac,
    set      = "OOF",
    pred     = oof,
    stringsAsFactors = FALSE
  )
  write.csv(pred_df, file.path(root, "predictions.csv"), row.names = FALSE)

  # 중요도 + 방향성 (SHAP 우선; 실패 시 gain+Spearman) — NA 방지
  shap_ok <- TRUE
  shap_mat <- try(predict(full_bst, dmat_all, predcontrib=TRUE), silent=TRUE)
  if (inherits(shap_mat, "try-error")) shap_ok <- FALSE
  feat_names <- colnames(X)
  Direction <- rep("Neutral", ncol(X)); names(Direction) <- feat_names

  if (shap_ok) {
    shap_mat <- matrix(shap_mat, ncol = ncol(X)+1, byrow = TRUE)
    colnames(shap_mat) <- c(feat_names, "BIAS")
    MeanAbsShap <- colMeans(abs(shap_mat[, feat_names, drop=FALSE]))
    MeanShap    <- colMeans(shap_mat[, feat_names, drop=FALSE])
    Direction[feat_names] <- ifelse(MeanShap > 0, "Positive",
                                    ifelse(MeanShap < 0, "Negative", "Neutral"))
    importance_feature <- data.frame(
      Feature=feat_names,
      PermImportance=MeanAbsShap,
      Direction=Direction,
      check.names=FALSE
    )
    importance_feature <- importance_feature[order(-importance_feature$PermImportance), ]
  } else {
    imp <- try(xgb.importance(model = full_bst), silent=TRUE)
    if (inherits(imp, "try-error") || is.null(imp)) {
      importance_feature <- data.frame(Feature = feat_names, PermImportance = NA_real_, stringsAsFactors = FALSE)
    } else {
      importance_feature <- imp[, c("Feature","Gain")]
      colnames(importance_feature) <- c("Feature","PermImportance")
    }
    safe_spearman <- function(a, b){
      if (is.null(a) || is.null(b)) return(NA_real_)
      if (all(is.na(a)) || all(is.na(b))) return(NA_real_)
      if (length(unique(a[!is.na(a)])) < 2) return(NA_real_)
      suppressWarnings(cor(a, b, method="spearman", use="pairwise.complete.obs"))
    }
    dir_val <- sapply(seq_len(ncol(X)), function(j) safe_spearman(X[,j], oof))
    dir_val[is.na(dir_val)] <- 0
    Direction <- ifelse(dir_val > 0, "Positive", ifelse(dir_val < 0, "Negative", "Neutral"))
    importance_feature$Direction <- Direction[importance_feature$Feature] %||% "Neutral"
  }

  base_feat <- colnames(relab)
  lab_vec <- if (!is.null(tax_map)) {
    ifelse(base_feat %in% names(tax_map), tax_map[base_feat], base_feat)
  } else base_feat
  feat2lab <- setNames(as.character(lab_vec), base_feat)

  importance_feature$Taxon <- ifelse(importance_feature$Feature %in% names(feat2lab),
                                     feat2lab[importance_feature$Feature],
                                     importance_feature$Feature)
  write.csv(importance_feature, file.path(root, "importance_feature.csv"), row.names = FALSE)

  importance_taxon <- importance_feature %>%
    group_by(Taxon) %>%
    summarise(
      SumImportance  = sum(PermImportance, na.rm=TRUE),
      MeanImportance = mean(PermImportance, na.rm=TRUE),
      Direction = {
        s <- sum(ifelse(Direction=="Positive",  1,
                        ifelse(Direction=="Negative", -1, 0)) * PermImportance, na.rm=TRUE)
        ifelse(s > 0, "Positive", ifelse(s < 0, "Negative", "Neutral"))
      },
      .groups="drop"
    ) %>% arrange(desc(SumImportance))
  write.csv(importance_taxon, file.path(root, "importance_taxon.csv"), row.names = FALSE)

  ## --- 9) 그림 저장 ----------------------------------------------------------
  png(file.path(root, "ROC_OOF.png"), width=900, height=700)
  plot.roc(roc_obj, main = sprintf("OOF ROC (AUC=%.3f)", AUC_oof), col="#1f77b4", lwd=2)
  abline(0,1,lty=2,col="grey60"); dev.off()

  png(file.path(root, "PR_OOF.png"), width=900, height=700)
  plot(pr, main = sprintf("OOF PR (AUPRC=%.3f)", pr$auc.integral)); dev.off()

  ## --- 10) Holdout 모드 ------------------------------------------------------
  results <- list()
  if (!testSet) {
    results$mode  <- "CV-only"
    results$AUC   <- AUC_oof
    results$AUPRC <- AUPRC_oof
    results$best_param <- best_param
    results$outdir <- root
    cat(sprintf("[CV-only] AUC=%.3f | AUPRC=%.3f\n", AUC_oof, AUPRC_oof))
    return(invisible(results))
  }

  ## Holdout 75/25 (그룹 보존 + 층화)
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
    tr_idx <- sort(c(tr_pos, tr_neg))
    te_idx <- setdiff(seq_along(ybin), tr_idx)
  }

  Xtr <- X[tr_idx,,drop=FALSE]; ytr <- ybin[tr_idx]; ytr_fac <- yfac[tr_idx]
  Xte <- X[te_idx,,drop=FALSE]; yte <- ybin[te_idx]; yte_fac <- yfac[te_idx]

  # inner-CV (train only)
  dtr_all <- xgb.DMatrix(Xtr, label=ytr, missing=NA)
  folds_tr <- if (group_mode) make_group_strat_folds(gid[tr_idx], ytr, n_folds, seed) else make_strat_folds(ytr, n_folds, seed)

  grid_in <- sample_param_grid(n_candidates)
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
  cat(sprintf("[InnerCV] best AUC=%.3f | best_iter=%d\n", best_in$auc, best_in$best_iter))

  # ★ Train OOF 예측으로 Train 성능 계산 (AUC=1 방지 핵심)
  oof_tr <- rep(NA_real_, length(ytr))
  for (k in seq_along(folds_tr)) {
    te_k <- folds_tr[[k]]; tr_k <- setdiff(seq_along(ytr), te_k)
    dtr_k <- xgb.DMatrix(Xtr[tr_k, , drop = FALSE], label = ytr[tr_k], missing = NA)
    dte_k <- xgb.DMatrix(Xtr[te_k, , drop = FALSE], label = ytr[te_k], missing = NA)
    bst_k <- xgb.train(params = best_in$param, data = dtr_k,
                       nrounds = best_in$best_iter, verbose = 0)
    oof_tr[te_k] <- predict(bst_k, dte_k)
  }
  roc_tr <- pROC::roc(ytr_fac, oof_tr, levels = orders, direction = "<", quiet = TRUE)
  auc_tr <- as.numeric(pROC::auc(roc_tr))
  pr_tr  <- PRROC::pr.curve(
    scores.class0 = oof_tr[ytr_fac != positive_class],
    scores.class1 = oof_tr[ytr_fac == positive_class],
    curve = TRUE
  )

  # 최종 모델로 TEST 예측
  dtr <- xgb.DMatrix(Xtr, label = ytr, missing = NA)
  dte <- xgb.DMatrix(Xte, label = yte, missing = NA)
  bst_tr <- xgb.train(params=best_in$param, data=dtr, nrounds=best_in$best_iter, verbose=0)
  ph_te <- predict(bst_tr, dte)

  roc_te <- pROC::roc(yte_fac, ph_te, levels = orders, direction = "<", quiet = TRUE)
  auc_te <- as.numeric(pROC::auc(roc_te))
  pr_te  <- PRROC::pr.curve(
    scores.class0 = ph_te[yte_fac != positive_class],
    scores.class1 = ph_te[yte_fac == positive_class],
    curve = TRUE
  )

  # 저장: 예측 (Train=OOF, Test=최종)
  pred_hold <- rbind(
    data.frame(SampleID=rownames(Xtr), StudyID=gid[tr_idx], outcome=ytr_fac, set="Train", pred=oof_tr),
    data.frame(SampleID=rownames(Xte), StudyID=gid[te_idx], outcome=yte_fac, set="Test",  pred=ph_te)
  )
  write.csv(pred_hold, file.path(root, "predictions.csv"), row.names = FALSE)

  # 중요도/방향성 (SHAP 우선; 실패 시 gain+Spearman) — Train 기준
  shap_ok2 <- TRUE
  shap_tr <- try(predict(bst_tr, dtr, predcontrib=TRUE), silent=TRUE)
  if (inherits(shap_tr, "try-error")) shap_ok2 <- FALSE

  if (shap_ok2) {
    shap_tr <- matrix(shap_tr, ncol = ncol(Xtr)+1, byrow = TRUE)
    colnames(shap_tr) <- c(colnames(Xtr), "BIAS")
    MeanAbsShap <- colMeans(abs(shap_tr[, colnames(Xtr), drop=FALSE]))
    MeanShap    <- colMeans(shap_tr[, colnames(Xtr), drop=FALSE])
    Direction_tr <- ifelse(MeanShap > 0, "Positive", ifelse(MeanShap < 0, "Negative", "Neutral"))
    importance_feature <- data.frame(
      Feature = colnames(Xtr),
      PermImportance = MeanAbsShap,
      Direction = Direction_tr,
      stringsAsFactors = FALSE
    )
  } else {
    safe_spearman <- function(a, b){
      if (is.null(a) || is.null(b)) return(NA_real_)
      if (all(is.na(a)) || all(is.na(b))) return(NA_real_)
      if (length(unique(a[!is.na(a)])) < 2) return(NA_real_)
      suppressWarnings(cor(a, b, method="spearman", use="pairwise.complete.obs"))
    }
    dir_val <- sapply(seq_len(ncol(Xtr)), function(j) safe_spearman(Xtr[,j], oof_tr))
    dir_val[is.na(dir_val)] <- 0
    Direction_tr <- ifelse(dir_val > 0, "Positive", ifelse(dir_val < 0, "Negative", "Neutral"))
    imp <- try(xgb.importance(model=bst_tr), silent=TRUE)
    if (inherits(imp, "try-error") || is.null(imp)) {
      importance_feature <- data.frame(Feature = colnames(Xtr), PermImportance = NA_real_, Direction = Direction_tr, stringsAsFactors = FALSE)
    } else {
      importance_feature <- imp[, c("Feature","Gain")]; colnames(importance_feature) <- c("Feature","PermImportance")
      importance_feature$Direction <- Direction_tr[importance_feature$Feature] %||% "Neutral"
    }
  }
  base_feat <- colnames(relab)
  lab_vec <- if (!is.null(tax_map)) { ifelse(base_feat %in% names(tax_map), tax_map[base_feat], base_feat) } else base_feat
  feat2lab <- setNames(as.character(lab_vec), base_feat)
  importance_feature$Taxon <- ifelse(importance_feature$Feature %in% names(feat2lab),
                                     feat2lab[importance_feature$Feature],
                                     importance_feature$Feature)
  importance_feature <- importance_feature[order(-importance_feature$PermImportance), ]
  write.csv(importance_feature, file.path(root, "importance_feature.csv"), row.names = FALSE)

  importance_taxon <- importance_feature %>%
    group_by(Taxon) %>%
    summarise(
      SumImportance  = sum(PermImportance, na.rm=TRUE),
      MeanImportance = mean(PermImportance, na.rm=TRUE),
      Direction = {
        s <- sum(ifelse(Direction=="Positive",  1,
                        ifelse(Direction=="Negative", -1, 0)) * PermImportance, na.rm=TRUE)
        ifelse(s > 0, "Positive", ifelse(s < 0, "Negative", "Neutral"))
      },
      .groups="drop"
    ) %>% arrange(desc(SumImportance))
  write.csv(importance_taxon, file.path(root, "importance_taxon.csv"), row.names = FALSE)

  # 그림
  png(file.path(root, "ROC_Train.png"), width=900, height=700)
  plot.roc(roc_tr,
           main = sprintf("Train ROC (AUC=%.3f)", auc_tr), col="#1f77b4", lwd=2)
  abline(0,1,lty=2,col="grey60"); dev.off()
  png(file.path(root, "ROC_Test.png"), width=900, height=700)
  plot.roc(roc_te,
           main = sprintf("Test ROC (AUC=%.3f)", auc_te), col="#d62728", lwd=2)
  abline(0,1,lty=2,col="grey60"); dev.off()
  png(file.path(root, "PR_Train.png"), width=900, height=700)
  plot(pr_tr, main = sprintf("Train PR (AUPRC=%.3f)", pr_tr$auc.integral)); dev.off()
  png(file.path(root, "PR_Test.png"), width=900, height=700)
  plot(pr_te, main = sprintf("Test PR (AUPRC=%.3f)", pr_te$auc.integral)); dev.off()

  results$mode <- "Holdout+CV"
  results$AUC_train <- auc_tr
  results$AUC_test  <- auc_te
  results$AUPRC_train <- pr_tr$auc.integral
  results$AUPRC_test  <- pr_te$auc.integral
  results$best_param <- best_in$param
  results$outdir <- root
  cat(sprintf("[Holdout] AUC(train)=%.3f | AUC(test)=%.3f | AUPRC(test)=%.3f\n", auc_tr, auc_te, pr_te$auc.integral))
  invisible(results)
}
