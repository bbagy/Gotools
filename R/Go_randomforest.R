#' Go_randomforest
#'
#' StudyID-aware Random Forest classification for microbiome (+ optional clinical covariates)
#' with stratified/grouped cross-validation, random hyperparameter search, and rich outputs.
#'
#' @description
#' Builds a random forest (via \pkg{ranger}) to predict a binary outcome from a \code{phyloseq}
#' object’s relative abundances (optionally aggregated to a taxonomic rank) and optional
#' clinical variables. The function:
#' \itemize{
#'   \item optionally aggregates features by \code{taxrank} (e.g., Genus),
#'   \item converts to relative abundances and filters by prevalence/mean abundance,
#'   \item one-hot encodes character clinical covariates, imputes NA with medians, and scales,
#'   \item fixes outcome level order (\code{orders}; positive class is \code{orders[2]}),
#'   \item performs random hyperparameter search using K-fold CV; folds preserve subject groups
#'         when \code{StudyID_col} varies across rows,
#'   \item trains the final model (CV-only or Holdout+CV) and saves predictions, importance,
#'         aggregated importance by taxon (with direction), and ROC/PR plots.
#' }
#'
#' @param psIN A \code{phyloseq} object.
#' @param project Character; project prefix used to create an output directory
#'   \verb{<project_YYMMDD>/Randomforest}.
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
#' @param num.trees Integer; number of trees for \pkg{ranger}. Default \code{1000}.
#' @param orders Character vector of length 2 giving the outcome level order,
#'   e.g., \code{c("Control","Case")}. The positive class is \code{orders[2]}.
#'
#' @details
#' \strong{Folds and grouping:} If \code{StudyID_col} does not equal rownames for all samples,
#' folds are created to keep samples from the same StudyID together, while preserving class balance
#' at the group level when possible. Otherwise, standard stratified folds are used.
#'
#' \strong{Filtering:} Features are kept if they pass \code{prev_min} OR \code{relab_min}.
#'
#' \strong{Direction of effect:} Feature “Direction” (Positive/Negative/Neutral) is assigned by
#' Spearman correlation between feature values and out-of-fold (or test) predicted probabilities.
#' \code{NA} correlations are treated as Neutral.
#'
#' @return (Invisibly) a list with run \code{mode}, metrics (AUC/AUPRC), \code{best_param},
#'   and \code{outdir}. Side effects: writes to \verb{<project_YYMMDD>/Randomforest/}:
#' \itemize{
#'   \item \code{random_search_results.csv}
#'   \item \code{predictions.csv} (OOF or Train/Test)
#'   \item \code{rf_final_model.rds}, \code{rf_meta.rds}
#'   \item \code{importance_feature.csv}, \code{importance_taxon.csv}
#'   \item ROC/PR PNGs (\code{ROC_OOF.png}, \code{PR_OOF.png}, or Train/Test variants)
#' }
#'
#' @examples
#' \dontrun{
#' # Binary outcome with levels c("Control","Case") in metadata column "Status"
#' res <- Go_randomforest(
#'   psIN = ps,
#'   project = "IBD",
#'   outcome = "Status",
#'   clinical_vari = c("Age","BMI","Sex"),
#'   StudyID_col = "SubjectID",
#'   taxrank = "Genus",
#'   prev_min = 0.05,
#'   relab_min = 1e-4,
#'   n_folds = 5,
#'   n_candidates = 30,
#'   num.trees = 1200,
#'   orders = c("Control","Case")
#' )
#' res$AUC
#' }
#'
#' @importFrom phyloseq sample_data tax_table taxa_are_rows otu_table tax_glom sample_sums
#' @importFrom ranger ranger importance
#' @importFrom pROC roc auc plot.roc
#' @importFrom PRROC pr.curve
#' @importFrom dplyr group_by summarise arrange
#' @importFrom stats cor model.matrix median as.formula
#' @importFrom utils write.csv
#' @importFrom grDevices png dev.off
#' @export


Go_randomforest <- function(
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
    n_candidates  = 40,
    num.trees     = 1000,
    orders        = c("Control","Case")   # 레벨 순서 (control, case)
){
  set.seed(seed)

  ## --- 0) 메타 준비: 레벨/StudyID 고정 (psIN에도 반영) -----------------------
  meta0 <- data.frame(sample_data(psIN), check.names = FALSE, stringsAsFactors = FALSE)

  # StudyID 없으면 샘플ID로 대체
  if (!(StudyID_col %in% names(meta0)) || all(is.na(meta0[[StudyID_col]]))) {
    message(sprintf("[Info] '%s'가 없어 rownames를 사용합니다.", StudyID_col))
    meta0[[StudyID_col]] <- rownames(meta0)
  }

  # outcome 레벨 고정 (orders 적용)
  stopifnot(outcome %in% names(meta0))
  yy0 <- meta0[[outcome]]
  if (!is.factor(yy0)) yy0 <- factor(yy0)

  if (!all(orders %in% levels(yy0))) {
    stop(sprintf("[Error] outcome(levels=%s)에 orders(%s) 중 일부가 없습니다.",
                 paste(levels(yy0), collapse=", "),
                 paste(orders, collapse=", ")))
  }
  yy0 <- factor(yy0, levels = orders)         # 레벨 고정
  meta0[[outcome]] <- yy0
  sample_data(psIN)[[outcome]] <- yy0         # psIN에도 반영 (이후 재로딩해도 유지)
  if (!(StudyID_col %in% colnames(sample_data(psIN)))) {
    sample_data(psIN)[[StudyID_col]] <- meta0[[StudyID_col]]
  }
  positive_class <- orders[2]

  ## --- 1) 출력 폴더 ----------------------------------------------------------
  today <- format(Sys.Date(), "%y%m%d")
  root  <- file.path(sprintf("%s_%s", project, today), "Randomforest")
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  message("[OutDir] ", root)

  ## --- 2) phyloseq → feature table ------------------------------------------
  stopifnot(inherits(psIN, "phyloseq"))
  ps <- psIN

  # taxrank 처리
  if (!identical(taxrank, "ASV")) {
    if (is.null(tax_table(ps)) || !(taxrank %in% colnames(tax_table(ps)))) {
      stop(sprintf("[Error] taxrank '%s' not found in tax_table.", taxrank))
    }
    ps <- tax_glom(ps, taxrank = taxrank, NArm = TRUE)
  }

  # feature table (samples x taxa)
  OTU <- as(otu_table(ps), "matrix")
  if (taxa_are_rows(ps)) OTU <- t(OTU)

  # 상대풍부도
  relab <- sweep(OTU, 1, pmax(1e-12, rowSums(OTU)), "/")

  # 필터링: prevalence & abundance
  prev <- colMeans(relab > 0)
  keep1 <- prev >= prev_min
  keep2 <- colMeans(relab) >= relab_min
  keep  <- keep1 | keep2
  relab <- relab[, keep, drop = FALSE]

  # Taxon 라벨 매핑
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

  # X: 마이크로바이옴 + 임상
  meta <- data.frame(sample_data(ps), check.names = FALSE, stringsAsFactors = FALSE)
  # 여기서 outcome/StudyID는 이미 psIN에 반영되어 있으므로 레벨이 유지됨
  X_micro <- relab

  # 임상 변수
  X_cli <- NULL
  if (length(clinical_vari)) {
    miss_cli <- setdiff(clinical_vari, colnames(meta))
    if (length(miss_cli)) {
      warning(sprintf("[Warn] clinical vars not found: %s", paste(miss_cli, collapse=", ")))
    }
    clinical_vari <- setdiff(clinical_vari, miss_cli)
    if (length(clinical_vari)) {
      X_cli <- meta[, clinical_vari, drop = FALSE]
      for (j in colnames(X_cli)) {
        if (is.factor(X_cli[[j]])) X_cli[[j]] <- as.character(X_cli[[j]])
      }
      char_cols <- names(Filter(is.character, X_cli))
      if (length(char_cols)) {
        dummies <- do.call(cbind, lapply(char_cols, function(cn){
          mm <- model.matrix(~.-1, data = data.frame(val = factor(X_cli[[cn]])))
          colnames(mm) <- paste(cn, colnames(mm), sep="__")
          mm
        }))
        X_cli <- cbind(
          X_cli[, setdiff(colnames(X_cli), char_cols), drop = FALSE],
          dummies
        )
      }
      X_cli <- data.matrix(X_cli)
      for (j in seq_len(ncol(X_cli))) {
        v <- X_cli[, j]
        if (anyNA(v)) v[is.na(v)] <- median(v, na.rm = TRUE)
        X_cli[, j] <- v
      }
      X_cli <- scale(X_cli)
    }
  }

  X <- if (!is.null(X_cli)) cbind(X_micro, X_cli) else X_micro
  X <- as.matrix(X); storage.mode(X) <- "double"

  # y, 그룹ID
  y   <- meta[[outcome]]             # factor with levels == orders
  gid <- meta[[StudyID_col]]
  stopifnot(length(y) == nrow(X), length(gid) == nrow(X))

  # class weights (불균형 시)
  taby <- table(y)
  class.weights <- as.numeric(max(taby) / taby)
  names(class.weights) <- names(taby)

  ## --- 3) 폴드 유틸 ----------------------------------------------------------
  make_group_strat_folds <- function(groups, y, K=5, seed=123){
    set.seed(seed)
    df <- data.frame(g=groups, y=as.integer(y==levels(y)[2]))
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
    ybin <- as.integer(y==levels(y)[2])
    idx0 <- which(ybin==0); idx1 <- which(ybin==1)
    k0 <- if (length(idx0) < K) max(2, length(idx0)) else K
    k1 <- if (length(idx1) < K) max(2, length(idx1)) else K
    f0 <- sample(rep(1:k0, length.out=length(idx0)))
    f1 <- sample(rep(1:k1, length.out=length(idx1)))
    ff <- integer(length(y))
    ff[idx0] <- f0; ff[idx1] <- f1
    lapply(1:max(ff), function(k) which(ff == k))
  }

  group_mode <- !all(meta[[StudyID_col]] == rownames(meta))

  ## --- 4) 랜덤 탐색 공간 -----------------------------------------------------
  sample_param_grid <- function(n){
    data.frame(
      mtry            = pmax(1, round(runif(n, 0.05, 0.5) * ncol(X))),  # 5~50% of features
      min.node.size   = sample(c(1,2,3,5,10), n, replace=TRUE),
      sample.fraction = runif(n, 0.6, 0.95),
      stringsAsFactors = FALSE
    )
  }

  ## --- 5) AUC 계산자 (levels 고정) ------------------------------------------
  get_auc <- function(obs, prob){
    ro <- pROC::roc(response = obs, predictor = prob,
                    levels = orders, direction = "<", quiet = TRUE)
    as.numeric(pROC::auc(ro))
  }

  ## --- 6) CV 평가(OOT/OOF) --------------------------------------------------
  eval_cv_auc <- function(X, y, groups=NULL, params, n_folds=5, seed=123){
    folds <- if (!is.null(groups)) make_group_strat_folds(groups, y, n_folds, seed)
    else make_strat_folds(y, n_folds, seed)
    oof <- rep(NA_real_, length(y))
    for (k in seq_along(folds)) {
      te <- folds[[k]]; tr <- setdiff(seq_along(y), te)
      df_tr <- data.frame(X[tr,,drop=FALSE])
      df_tr[[outcome]] <- y[tr]

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
      phat <- predict(rf, data.frame(X[te,,drop=FALSE]))$predictions[, positive_class]
      oof[te] <- phat
    }
    auc <- get_auc(y, oof)
    list(auc=auc, pred=oof, folds=folds)
  }

  ## --- 7) 랜덤 탐색 ----------------------------------------------------------
  grid <- sample_param_grid(n_candidates)
  evals <- vector("list", nrow(grid))
  best_auc <- -Inf; best_ix <- 1

  cat(sprintf("[TUNE] random search: %d candidates\n", nrow(grid)))
  for (i in seq_len(nrow(grid))) {
    params_i <- grid[i, ]
    cv <- eval_cv_auc(X, y, if (group_mode) gid else NULL, params_i, n_folds, seed)
    evals[[i]] <- cbind(params_i, AUC=cv$auc)
    if (!is.na(cv$auc) && cv$auc > best_auc) {
      best_auc <- cv$auc; best_ix <- i
    }
    cat(sprintf("  - cand %2d/%2d: AUC=%.3f (mtry=%d, min.node=%d, frac=%.2f)\n",
                i, nrow(grid), cv$auc, params_i$mtry, params_i$min.node.size, params_i$sample.fraction))
  }
  eval_df <- do.call(rbind, evals)
  eval_df <- eval_df[order(-eval_df$AUC), ]
  write.csv(eval_df, file.path(root, "random_search_results.csv"), row.names = FALSE)
  best_param <- grid[best_ix, ]
  cat(sprintf("[Best CV] AUC=%.3f | mtry=%d, min.node=%d, frac=%.2f\n",
              best_auc, best_param$mtry, best_param$min.node.size, best_param$sample.fraction))

  ## --- 8) 최종 학습/평가 + 저장 ----------------------------------------------
  results <- list()
  if (!testSet) {
    # CV-only: best_param으로 OOF 예측 생성
    cv <- eval_cv_auc(X, y, if (group_mode) gid else NULL, best_param, n_folds, seed)
    oof <- cv$pred
    auc_roc <- get_auc(y, oof)

    # PR
    pr <- PRROC::pr.curve(scores.class0 = oof[y==positive_class],
                          scores.class1 = oof[y!=positive_class], curve = TRUE)
    auprc <- pr$auc.integral

    # 최종 모델 (전체로)
    df_all <- data.frame(X); df_all[[outcome]] <- y
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
      keep.inbag = TRUE,
      seed = seed
    )

    # 저장: 모델/메타/예측
    saveRDS(final_rf, file.path(root, "rf_final_model.rds"))
    meta_obj <- list(
      outcome = outcome,
      StudyID_col = StudyID_col,
      positive_class = positive_class,
      best_param = best_param,
      num.trees = num.trees,
      folds = cv$folds,
      seed = seed,
      levels = orders
    )
    saveRDS(meta_obj, file.path(root, "rf_meta.rds"))

    pred_df <- data.frame(
      SampleID = rownames(X),
      StudyID  = gid,
      outcome  = y,
      set      = "OOF",
      pred     = oof,
      stringsAsFactors = FALSE
    )
    write.csv(pred_df, file.path(root, "predictions.csv"), row.names = FALSE)

    # 중요도 + 방향성 (OOF 기반) — NA 방지(Neutral)
    safe_spearman <- function(a, b){
      if (is.null(a) || is.null(b)) return(NA_real_)
      if (all(is.na(a)) || all(is.na(b))) return(NA_real_)
      if (length(unique(a[!is.na(a)])) < 2) return(NA_real_)  # 상수열
      suppressWarnings(cor(a, b, method="spearman", use="pairwise.complete.obs"))
    }
    dir_val <- sapply(seq_len(ncol(X)), function(j) safe_spearman(X[,j], oof))
    names(dir_val) <- colnames(X)
    dir_val[is.na(dir_val)] <- 0  # NA → 0
    dir_lab <- ifelse(dir_val > 0, "Positive", ifelse(dir_val < 0, "Negative", "Neutral"))

    imp <- as.data.frame(ranger::importance(final_rf))
    colnames(imp) <- "PermImportance"
    imp$Feature <- rownames(imp)
    imp$Direction <- dir_lab[imp$Feature] %||% "Neutral"

    base_feat <- colnames(relab)
    lab_vec <- if (!is.null(tax_map)) {
      ifelse(base_feat %in% names(tax_map), tax_map[base_feat], base_feat)
    } else base_feat
    feat2lab <- setNames(as.character(lab_vec), base_feat)

    imp$Taxon <- ifelse(imp$Feature %in% names(feat2lab), feat2lab[imp$Feature], imp$Feature)
    imp <- imp[order(-imp$PermImportance), ]
    write.csv(imp, file.path(root, "importance_feature.csv"), row.names = FALSE)

    agg <- imp %>%
      group_by(Taxon) %>%
      summarise(
        SumImportance  = sum(PermImportance, na.rm=TRUE),
        MeanImportance = mean(PermImportance, na.rm=TRUE),
        # Direction 가중 평균 부호로 결정, 0이면 Neutral
        Direction = {
          s <- sum(ifelse(Direction=="Positive",  1,
                          ifelse(Direction=="Negative", -1, 0)) * PermImportance,
                   na.rm=TRUE)
          ifelse(s > 0, "Positive", ifelse(s < 0, "Negative", "Neutral"))
        },
        .groups="drop"
      ) %>% arrange(desc(SumImportance))
    write.csv(agg, file.path(root, "importance_taxon.csv"), row.names = FALSE)

    # 간단 ROC/PR PNG (levels=orders)
    roc_obj <- pROC::roc(y, oof, levels = orders, direction="<", quiet=TRUE)
    png(file.path(root, "ROC_OOF.png"), width=900, height=700)
    plot.roc(roc_obj, main = sprintf("OOF ROC (AUC=%.3f)", as.numeric(pROC::auc(roc_obj))), col="#1f77b4", lwd=2)
    abline(0,1,lty=2,col="grey60"); dev.off()

    png(file.path(root, "PR_OOF.png"), width=900, height=700)
    plot(pr, main = sprintf("OOF PR (AUPRC=%.3f)", pr$auc.integral))
    dev.off()

    results$mode  <- "CV-only"
    results$AUC   <- auc_roc
    results$AUPRC <- auprc
    results$best_param <- best_param
    results$outdir <- root

    cat(sprintf("[CV-only] AUC=%.3f | AUPRC=%.3f\n", auc_roc, auprc))

  } else {
    ## Holdout 75/25 (그룹 보존 + 층화)
    set.seed(seed)
    if (group_mode) {
      df_g <- data.frame(gid=gid, y=as.integer(y==positive_class))
      rep_lab <- aggregate(y ~ gid, df_g, function(z) round(mean(z)))
      labs <- rep_lab$gid
      g0 <- rep_lab$gid[rep_lab$y==0]; g1 <- rep_lab$gid[rep_lab$y==1]
      n0_tr <- floor(length(g0)*0.75); n1_tr <- floor(length(g1)*0.75)
      tr_g <- c(sample(g0, n0_tr), sample(g1, n1_tr))
      te_g <- setdiff(labs, tr_g)
      tr_idx <- which(gid %in% tr_g); te_idx <- which(gid %in% te_g)
    } else {
      pos_idx <- which(y==positive_class); neg_idx <- which(y!=positive_class)
      tr_pos <- sample(pos_idx, floor(0.75*length(pos_idx)))
      tr_neg <- sample(neg_idx, floor(0.75*length(neg_idx)))
      tr_idx <- sort(c(tr_pos, tr_neg))
      te_idx <- setdiff(seq_along(y), tr_idx)
    }

    Xtr <- X[tr_idx,,drop=FALSE]; ytr <- y[tr_idx]; gid_tr <- gid[tr_idx]
    Xte <- X[te_idx,,drop=FALSE]; yte <- y[te_idx]; gid_te <- gid[te_idx]

    # inner-CV (훈련세트만)
    grid_in <- sample_param_grid(n_candidates)
    best_auc_in <- -Inf; best_ix_in <- 1
    for (i in seq_len(nrow(grid_in))) {
      params_i <- grid_in[i, ]
      cv <- eval_cv_auc(Xtr, ytr, if (group_mode) gid_tr else NULL, params_i, n_folds, seed)
      if (!is.na(cv$auc) && cv$auc > best_auc_in) { best_auc_in <- cv$auc; best_ix_in <- i }
    }
    best_param <- grid_in[best_ix_in, ]
    cat(sprintf("[InnerCV] best AUC=%.3f | mtry=%d, min.node=%d, frac=%.2f\n",
                best_auc_in, best_param$mtry, best_param$min.node.size, best_param$sample.fraction))

    # 최종 훈련 (훈련 세트 전체)
    df_tr <- data.frame(Xtr); df_tr[[outcome]] <- ytr
    final_rf <- ranger::ranger(
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
      seed = seed
    )
    saveRDS(final_rf, file.path(root, "rf_final_model.rds"))
    meta_obj <- list(
      outcome = outcome,
      StudyID_col = StudyID_col,
      positive_class = positive_class,
      best_param = best_param,
      num.trees = num.trees,
      seed = seed,
      split = list(train_idx=tr_idx, test_idx=te_idx),
      levels = orders
    )
    saveRDS(meta_obj, file.path(root, "rf_meta.rds"))

    # 예측
    ph_tr <- final_rf$predictions[, positive_class]
    ph_te <- predict(final_rf, data.frame(Xte))$predictions[, positive_class]

    # ROC/PR
    auc_tr <- get_auc(ytr, ph_tr)
    auc_te <- get_auc(yte, ph_te)
    pr_tr <- PRROC::pr.curve(scores.class0 = ph_tr[ytr==positive_class],
                             scores.class1 = ph_tr[ytr!=positive_class], curve=TRUE)
    pr_te <- PRROC::pr.curve(scores.class0 = ph_te[yte==positive_class],
                             scores.class1 = ph_te[yte!=positive_class], curve=TRUE)

    # 저장: 예측
    pred_df <- rbind(
      data.frame(SampleID=rownames(Xtr), StudyID=gid_tr, outcome=ytr, set="Train", pred=ph_tr),
      data.frame(SampleID=rownames(Xte), StudyID=gid_te, outcome=yte, set="Test",  pred=ph_te)
    )
    write.csv(pred_df, file.path(root, "predictions.csv"), row.names = FALSE)

    # 중요도 + 방향성 (테스트 기준; NA→Neutral)
    safe_spearman <- function(a, b){
      if (is.null(a) || is.null(b)) return(NA_real_)
      if (all(is.na(a)) || all(is.na(b))) return(NA_real_)
      if (length(unique(a[!is.na(a)])) < 2) return(NA_real_)
      suppressWarnings(cor(a, b, method="spearman", use="pairwise.complete.obs"))
    }
    dir_base <- if (length(ph_te) >= 5) ph_te else ph_tr
    X_base   <- if (length(ph_te) >= 5) Xte   else Xtr
    dir_val  <- sapply(seq_len(ncol(X_base)), function(j) safe_spearman(X_base[,j], dir_base))
    names(dir_val) <- colnames(X_base)

    # X_base에 없던 피처 보정
    miss_dir <- setdiff(colnames(X), names(dir_val))
    if (length(miss_dir)) {
      add_val <- sapply(miss_dir, function(feat) safe_spearman(Xtr[,feat], ph_tr))
      dir_val <- c(dir_val, add_val)
    }
    dir_val[is.na(dir_val)] <- 0
    dir_lab <- ifelse(dir_val > 0, "Positive", ifelse(dir_val < 0, "Negative", "Neutral"))

    imp <- as.data.frame(ranger::importance(final_rf))
    colnames(imp) <- "PermImportance"
    imp$Feature <- rownames(imp)
    imp$Direction <- dir_lab[imp$Feature] %||% "Neutral"

    base_feat <- colnames(relab)
    lab_vec <- if (!is.null(tax_map)) {
      ifelse(base_feat %in% names(tax_map), tax_map[base_feat], base_feat)
    } else base_feat
    feat2lab <- setNames(as.character(lab_vec), base_feat)

    imp$Taxon <- ifelse(imp$Feature %in% names(feat2lab), feat2lab[imp$Feature], imp$Feature)
    imp <- imp[order(-imp$PermImportance), ]
    write.csv(imp, file.path(root, "importance_feature.csv"), row.names = FALSE)

    agg <- imp %>%
      group_by(Taxon) %>%
      summarise(
        SumImportance  = sum(PermImportance, na.rm=TRUE),
        MeanImportance = mean(PermImportance, na.rm=TRUE),
        Direction = {
          s <- sum(ifelse(Direction=="Positive",  1,
                          ifelse(Direction=="Negative", -1, 0)) * PermImportance,
                   na.rm=TRUE)
          ifelse(s > 0, "Positive", ifelse(s < 0, "Negative", "Neutral"))
        },
        .groups="drop"
      ) %>% arrange(desc(SumImportance))
    write.csv(agg, file.path(root, "importance_taxon.csv"), row.names = FALSE)

    # ROC/PR 그림 (levels=orders)
    roc_tr <- pROC::roc(ytr, ph_tr, levels=orders, direction="<", quiet=TRUE)
    roc_te <- pROC::roc(yte, ph_te, levels=orders, direction="<", quiet=TRUE)
    png(file.path(root, "ROC_Train.png"), width=900, height=700)
    plot.roc(roc_tr, main = sprintf("Train ROC (AUC=%.3f)", as.numeric(pROC::auc(roc_tr))), col="#1f77b4", lwd=2)
    abline(0,1,lty=2,col="grey60"); dev.off()
    png(file.path(root, "ROC_Test.png"), width=900, height=700)
    plot.roc(roc_te, main = sprintf("Test ROC (AUC=%.3f)", as.numeric(pROC::auc(roc_te))), col="#d62728", lwd=2)
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
    results$best_param <- best_param
    results$outdir <- root

    cat(sprintf("[Holdout] AUC(train)=%.3f | AUC(test)=%.3f | AUPRC(test)=%.3f\n",
                auc_tr, auc_te, pr_te$auc.integral))
  }

  invisible(results)
}

# 작은 유틸: NULL coalescing
`%||%` <- function(a, b) if (!is.null(a)) a else b
