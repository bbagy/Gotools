## =========================================================
## Go_randomforest(): phyloseq → RandomForest full pipeline
##  - UID-wise 누수방지 분할
##  - ASV 필터 + CLR
##  - (옵션) 임상변수 머지(표준화/원핫/결측대치)
##  - 랜덤탐색 튜닝 + OOF/TEST ROC/PR
##  - SHAP 기반 중요도(방향성 포함), ASV→TaxLabel 라벨링
##  - 모든 결과 저장
## =========================================================

suppressPackageStartupMessages({
  library(phyloseq)
  library(ranger)
  library(pROC)
  library(PRROC)
  library(ggplot2)
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(fastshap)
})

Go_randomforest <- function(psIN,
                            project,
                            outcome,                 # sample_data() 컬럼명(이진)
                            testSet = FALSE,
                            clinical_vari = c(),     # sample_data()의 임상변수명들(옵션)
                            uid_col = "UID",
                            K = 5,
                            seed = 123,
                            prev_thr = 0.05,         # prevalence >= 5%
                            mean_abund_thr = 1e-4,   # mean rel.abund >= 0.01%
                            pseudo = 1e-6,           # CLR pseudo-count
                            n_random = 60,           # 랜덤 탐색 횟수
                            num.trees = 1000,
                            num.threads = NULL) {
  
  set.seed(seed)
  
  ## ---------------- namespace-like helpers ----------------
  RFx <- new.env(parent = emptyenv())
  
  RFx$safe_dir <- function(...) {
    p <- file.path(...)
    if (!dir.exists(p)) dir.create(p, recursive = TRUE, showWarnings = FALSE)
    p
  }
  
  RFx$today_tag <- function() format(Sys.Date(), "%y%m%d")
  
  RFx$levels_01 <- function(yfac) {
    # control=0, case=1 로 고정 (알파벳순으로 0/1 배정)
    y <- factor(yfac)
    if (nlevels(y) != 2) stop("Outcome must be binary.")
    y <- factor(y, levels = sort(levels(y)))
    list(y = as.integer(y == levels(y)[2]), levels = levels(y))
  }
  
  RFx$get_uid <- function(sdt, uid_col) {
    sdt <- as.data.frame(sdt)
    if (!uid_col %in% colnames(sdt)) return(NULL)
    as.character(sdt[[uid_col]])
  }
  
  RFx$asv_filter_clr <- function(ps, prev_thr, mean_abund_thr, pseudo) {
    # 1) prevalence filter
    otu <- as(otu_table(ps), "matrix")
    if (taxa_are_rows(ps)) otu <- t(otu)
    # prevalence
    prev <- colMeans(otu > 0, na.rm = TRUE)
    keep1 <- prev >= prev_thr
    # rel.abund & mean_abund
    rel <- sweep(otu, 1, pmax(rowSums(otu), 1), "/")
    mean_ab <- colMeans(rel, na.rm = TRUE)
    keep2 <- mean_ab >= mean_abund_thr
    keep <- keep1 & keep2
    if (!any(keep)) stop("No taxa left after filtering. Loosen thresholds.")
    otu_f <- rel[, keep, drop = FALSE]
    # CLR
    otu_f <- otu_f + pseudo
    gm <- exp(rowMeans(log(otu_f)))
    clr <- log(otu_f) - log(gm)
    # 행 = sample, 열 = taxa(ASV)
    list(X = clr, taxa = colnames(clr))
  }
  
  RFx$build_tax_label <- function(ps, asv_vec) {
    tt <- as.data.frame(tax_table(ps), stringsAsFactors = FALSE)
    tt$ASV <- rownames(tt)
    tt <- tt[tt$ASV %in% asv_vec, , drop = FALSE]
    prio <- intersect(c("Species","Genus","Family","Order","Class","Phylum","Kingdom"),
                      colnames(tt))
    make_lab <- function(row) {
      g <- row[["Genus"]]; s <- row[["Species"]]
      if (!is.na(g) && nzchar(g) && !is.na(s) && nzchar(s)) return(paste(g, s))
      if (!is.na(g) && nzchar(g)) return(g)
      if (!is.na(s) && nzchar(s)) return(s)
      for (nm in prio) {
        v <- row[[nm]]
        if (!is.na(v) && nzchar(v)) return(paste0(v, " sp."))
      }
      row[["ASV"]]
    }
    tt$TaxLabel <- apply(tt, 1, make_lab)
    setNames(tt$TaxLabel, tt$ASV)
  }
  
  RFx$prep_clinical <- function(sdt, clinical_vari, train_ids) {
    # 반환: list( Xc_train, Xc_all, scaler, levels_map, cols )
    if (length(clinical_vari) == 0) {
      return(list(Xc_all = NULL, Xc_train = NULL, scaler = NULL, levels_map = NULL, cols = NULL))
    }
    df <- as.data.frame(sdt)[, clinical_vari, drop = FALSE]
    # 구분
    num_cols <- names(df)[sapply(df, is.numeric)]
    cat_cols <- setdiff(names(df), num_cols)
    
    # 수치: 훈련셋 기준 평균/표준편차
    scaler <- NULL; levels_map <- NULL
    if (length(num_cols)) {
      mu <- sapply(df[train_ids, num_cols, drop=FALSE], function(z) mean(z, na.rm=TRUE))
      sdv <- sapply(df[train_ids, num_cols, drop=FALSE], function(z) sd(z, na.rm=TRUE))
      sdv[sdv == 0 | is.na(sdv)] <- 1
      # 결측 대치(평균)
      for (nm in num_cols) df[[nm]][is.na(df[[nm]])] <- mu[[nm]]
      # 표준화
      for (nm in num_cols) df[[nm]] <- (df[[nm]] - mu[[nm]]) / sdv[[nm]]
      scaler <- list(mu = mu, sd = sdv, cols = num_cols)
    }
    
    # 범주: 훈련셋 레벨 고정 + One-hot (결측은 "Missing")
    Xc <- NULL
    if (length(cat_cols)) {
      df[cat_cols] <- lapply(df[cat_cols], function(z) {
        z <- as.character(z); z[is.na(z) | z==""] <- "(Missing)"; factor(z)
      })
      levels_map <- lapply(df[train_ids, cat_cols, drop=FALSE], levels)
      # test 쪽 unseen level은 "(Missing)"로
      for (nm in cat_cols) {
        lv <- levels_map[[nm]]
        df[[nm]][!(df[[nm]] %in% lv)] <- "(Missing)"
        df[[nm]] <- factor(df[[nm]], levels = unique(c(lv, "(Missing)")))
      }
      Xc <- model.matrix(~ . - 1, data = df[, cat_cols, drop = FALSE])
    }
    
    # 합치기
    if (length(num_cols) && !is.null(Xc)) {
      Xc_all <- cbind(as.matrix(df[, num_cols, drop=FALSE]), Xc)
      cols <- c(num_cols, colnames(Xc))
    } else if (length(num_cols)) {
      Xc_all <- as.matrix(df[, num_cols, drop=FALSE]); cols <- num_cols
    } else if (!is.null(Xc)) {
      Xc_all <- Xc; cols <- colnames(Xc)
    } else {
      Xc_all <- NULL; cols <- NULL
    }
    list(Xc_all = Xc_all, Xc_train = if (!is.null(Xc_all)) Xc_all[train_ids, , drop=FALSE] else NULL,
         scaler = scaler, levels_map = levels_map, cols = cols)
  }
  
  RFx$uid_stratify <- function(uid_vec, y01, K, seed) {
    set.seed(seed)
    uid <- as.character(uid_vec)
    dfu <- aggregate(y01 ~ uid, FUN = function(z) round(mean(z)))
    colnames(dfu) <- c("uid", "label")
    u0 <- dfu$uid[dfu$label == 0]
    u1 <- dfu$uid[dfu$label == 1]
    k0 <- if (length(u0) < K) max(2, length(u0)) else K
    k1 <- if (length(u1) < K) max(2, length(u1)) else K
    f0 <- sample(rep(1:k0, length.out = length(u0)))
    f1 <- sample(rep(1:k1, length.out = length(u1)))
    fmap <- rbind(data.frame(uid = u0, fold = f0),
                  data.frame(uid = u1, fold = f1))
    ff <- fmap$fold[match(uid, fmap$uid)]
    ff[is.na(ff)] <- sample(1:max(fmap$fold), sum(is.na(ff)), replace = TRUE)
    lapply(1:max(ff), function(k) which(ff == k))
  }
  
  RFx$auc_pr <- function(y, p) {
    roc_obj <- pROC::roc(response = y, predictor = p, quiet = TRUE, direction = "<")
    auc_roc <- as.numeric(pROC::auc(roc_obj))
    pr_obj  <- PRROC::pr.curve(scores.class0 = p[y == 1], scores.class1 = p[y == 0], curve = FALSE)
    list(auc = auc_roc, auprc = as.numeric(pr_obj$auc.integral))
  }
  
  RFx$rand_params <- function(p, n_random) {
    # p = feature 수 (mtry 샘플링 범위에 사용)
    data.frame(
      mtry = pmax(1L, round(runif(n_random, 0.5*sqrt(p), 2*sqrt(p)))),
      min.node.size = sample(c(1,3,5,10), n_random, replace = TRUE),
      sample.fraction = runif(n_random, 0.6, 1.0),
      splitrule = sample(c("gini","extratrees"), n_random, replace = TRUE),
      stringsAsFactors = FALSE
    )
  }
  
  RFx$oof_pred_ranger <- function(X, y, folds, class.weights, num.trees, parms, num.threads) {
    # 각 fold에서 OOF 확률 생성
    pred <- rep(NA_real_, length(y))
    for (k in seq_along(folds)) {
      te <- folds[[k]]; tr <- setdiff(seq_along(y), te)
      df_tr <- data.frame(y = factor(y, levels = c(0,1))[tr], X[tr, , drop=FALSE])
      fit <- ranger::ranger(
        y ~ .,
        data = df_tr,
        probability = TRUE,
        num.trees = num.trees,
        mtry = parms$mtry,
        min.node.size = parms$min.node.size,
        sample.fraction = parms$sample.fraction,
        splitrule = parms$splitrule,
        class.weights = class.weights,
        num.threads = num.threads,
        seed = 1
      )
      pr <- predict(fit, data = data.frame(X[te, , drop=FALSE]))$predictions
      # 두 레벨 prob 행렬 → case(=level 2) 확률
      p1 <- if (is.matrix(pr)) pr[,2] else pr
      pred[te] <- p1
    }
    pred
  }
  
  RFx$fit_final_ranger <- function(X, y, class.weights, num.trees, parms, num.threads) {
    df_tr <- data.frame(y = factor(y, levels = c(0,1)), X)
    ranger::ranger(
      y ~ .,
      data = df_tr,
      probability = TRUE,
      num.trees = num.trees,
      mtry = parms$mtry,
      min.node.size = parms$min.node.size,
      sample.fraction = parms$sample.fraction,
      splitrule = parms$splitrule,
      class.weights = class.weights,
      num.threads = num.threads,
      seed = 1,
      importance = "impurity"
    )
  }
  
  RFx$shap_table <- function(model, X, ps=NULL, topN=20) {
    # fastshap: 적용 데이터 일부 샘플링 가능(여기선 전부 사용)
    set.seed(1)
    ffun <- function(object, newdata) {
      p <- predict(object, data = data.frame(newdata))$predictions
      if (is.matrix(p)) p[,2] else p
    }
    S <- fastshap::explain(
      object = model,
      X = as.data.frame(X),
      pred_wrapper = ffun,
      nsim = 100, # 빠른 근사
      adjust = TRUE
    )
    MeanAbsShap <- apply(abs(S), 2, mean)
    MeanShap    <- colMeans(S)
    df <- data.frame(
      Feature = colnames(X),
      MeanAbsShap = MeanAbsShap,
      MeanShap = MeanShap,
      stringsAsFactors = FALSE
    )
    # ASV 라벨
    base_asv <- sub("^.*__", "", df$Feature)      # tp prefix 제거가 있을 경우 대비
    base_asv <- sub("__d.*$", "", base_asv)       # Δ 꼬리표 제거
    df$BaseASV <- base_asv
    if (!is.null(ps)) {
      labmap <- RFx$build_tax_label(ps, unique(df$BaseASV))
      df$TaxLabel <- ifelse(df$BaseASV %in% names(labmap), labmap[df$BaseASV], df$BaseASV)
    } else {
      df$TaxLabel <- df$BaseASV
    }
    df$Direction <- ifelse(df$MeanShap >= 0, "Positive", "Negative")
    df <- df[order(-df$MeanAbsShap), ]
    list(shap = S, table = df, top = head(df, topN))
  }
  
  RFx$plot_bar <- function(df_top, title, out_png) {
    p <- ggplot(df_top, aes(x = reorder(TaxLabel, MeanAbsShap), y = MeanAbsShap, fill = Direction)) +
      geom_col() + coord_flip() +
      labs(title = title, x = NULL, y = "Mean |SHAP|") +
      theme_bw(base_size = 12)
    ggsave(out_png, p, width = 6, height = 5, dpi = 300)
    p
  }
  
  RFx$plot_rocpr <- function(y, p, out_png_prefix, title_prefix) {
    roc_obj <- pROC::roc(response = y, predictor = p, quiet=TRUE, direction = "<")
    auc_roc <- as.numeric(pROC::auc(roc_obj))
    png(paste0(out_png_prefix, "_ROC.png"), width=1400, height=1200, res=180)
    plot(roc_obj, col="#1f77b4", lwd=2, main=sprintf("%s ROC (AUC=%.3f)", title_prefix, auc_roc))
    abline(0,1,lty=2,col="grey60"); dev.off()
    
    pr_obj <- PRROC::pr.curve(scores.class0 = p[y==1], scores.class1 = p[y==0], curve=TRUE)
    png(paste0(out_png_prefix, "_PR.png"), width=1400, height=1200, res=180)
    plot(pr_obj, main=sprintf("%s PR (AUPRC=%.3f)", title_prefix, pr_obj$auc.integral))
    dev.off()
    list(auc=auc_roc, auprc=as.numeric(pr_obj$auc.integral))
  }
  
  ## ---------------- output dirs ----------------
  out_root <- sprintf("%s_%s", project, RFx$today_tag())
  out_dir  <- RFx$safe_dir(out_root, "Randomforest")
  out_tab  <- RFx$safe_dir(out_dir, "table")
  out_fig  <- RFx$safe_dir(out_dir, "fig")
  out_rds  <- RFx$safe_dir(out_dir, "rds")
  message("[OUT] ", normalizePath(out_dir))
  
  ## ---------------- data prep ----------------
  stopifnot(inherits(psIN, "phyloseq"))
  sdt <- sample_data(psIN)
  if (!outcome %in% colnames(sdt)) stop("Outcome column not found in sample_data().")
  
  # 라벨(0/1), UID
  yinfo <- RFx$levels_01(sdt[[outcome]])
  y01   <- yinfo$y
  ylv   <- yinfo$levels
  uid_vec <- RFx$get_uid(sdt, uid_col)
  if (is.null(uid_vec)) uid_vec <- rownames(sdt)  # UID 없으면 샘플ID로 대체(누수방지X)
  
  # ASV 필터 + CLR
  x_asv <- RFx$asv_filter_clr(psIN, prev_thr, mean_abund_thr, pseudo)
  Xmb   <- x_asv$X  # samples x features
  # 행 맞추기(안전)
  common_samples <- intersect(rownames(sdt), rownames(Xmb))
  sdt <- sdt[common_samples, , drop=FALSE]
  Xmb <- Xmb[common_samples, , drop=FALSE]
  y01 <- y01[match(common_samples, rownames(sdt))]
  uid_vec <- uid_vec[match(common_samples, rownames(sdt))]
  
  # 임상 변수
  tr_ids_all <- seq_len(nrow(sdt)) # 임시 (train id는 split 후 업데이트)
  clin <- RFx$prep_clinical(sdt, clinical_vari, train_ids = tr_ids_all)
  Xall <- if (!is.null(clin$Xc_all)) cbind(Xmb, clin$Xc_all) else Xmb
  colnames(Xall) <- make.unique(colnames(Xall))
  
  # 클래스 가중치
  pos <- sum(y01==1); neg <- sum(y01==0)
  class.weights <- c("0" = 1, "1" = if (pos>0) neg/pos else 1)
  
  # num.threads
  if (is.null(num.threads)) {
    num.threads <- max(1L, parallel::detectCores(logical = TRUE) - 1L)
  }
  
  ## ---------------- split (UID-wise) ----------------
  if (!testSet) {
    tr_idx <- seq_len(nrow(Xall)); te_idx <- integer(0)
  } else {
    set.seed(seed)
    # UID-wise stratified 75/25
    uids <- unique(uid_vec)
    # 라벨 by UID
    dfu <- aggregate(y01 ~ uid_vec, FUN = function(z) round(mean(z)))
    names(dfu) <- c("uid","lab")
    # 각 라벨 그룹에서 75%
    pick_uid <- function(U) sample(U, size = max(1, floor(0.75*length(U))), replace = FALSE)
    u_tr <- c(pick_uid(dfu$uid[dfu$lab==0]), pick_uid(dfu$uid[dfu$lab==1]))
    is_tr <- uid_vec %in% u_tr
    tr_idx <- which(is_tr); te_idx <- which(!is_tr)
    if (length(unique(y01[te_idx])) < 2) {
      # fallback: 샘플 wise로 75/25
      n <- length(y01); tr_idx <- sample(seq_len(n), size = floor(0.75*n))
      te_idx <- setdiff(seq_len(n), tr_idx)
    }
  }
  
  # CV folds (train set only, UID-wise)
  folds <- RFx$uid_stratify(uid_vec[tr_idx], y01[tr_idx], K = K, seed = seed)
  
  ## ---------------- random search & OOF ----------------
  P <- ncol(Xall)
  grid <- RFx$rand_params(P, n_random)
  oof_best <- -Inf
  best_param <- NULL
  OOF_pred_best <- NULL
  
  for (i in seq_len(nrow(grid))) {
    parms <- grid[i, ]
    pred <- RFx$oof_pred_ranger(
      X = Xall[tr_idx, , drop=FALSE],
      y = y01[tr_idx],
      folds = folds,
      class.weights = class.weights,
      num.trees = num.trees,
      parms = parms,
      num.threads = num.threads
    )
    mets <- RFx$auc_pr(y01[tr_idx], pred)
    if (mets$auc > oof_best || (abs(mets$auc - oof_best) < 1e-6 && mets$auprc > (attr(OOF_pred_best,"auprc") %||% -Inf))) {
      oof_best <- mets$auc
      best_param <- parms
      OOF_pred_best <- pred
      attr(OOF_pred_best,"auprc") <- mets$auprc
    }
    cat(sprintf("[RS %3d/%3d] AUC=%.3f, AUPRC=%.3f  (mtry=%d, min.node=%d, frac=%.2f, split=%s)\n",
                i, nrow(grid), mets$auc, mets$auprc,
                parms$mtry, parms$min.node.size, parms$sample.fraction, parms$splitrule))
  }
  cat(sprintf("\n[Best OOF] AUC=%.3f, AUPRC=%.3f\n", oof_best, attr(OOF_pred_best,"auprc")))
  
  ## ---------------- train final & (option) test ----------------
  final_model <- RFx$fit_final_ranger(
    X = Xall[tr_idx, , drop=FALSE],
    y = y01[tr_idx],
    class.weights = class.weights,
    num.trees = num.trees,
    parms = best_param,
    num.threads = num.threads
  )
  
  # OOF curves
  oof_curves <- RFx$plot_rocpr(y01[tr_idx], OOF_pred_best,
                               out_png_prefix = file.path(out_fig,"OOF"),
                               title_prefix = "OOF")
  
  # TEST curves
  test_curves <- NULL; test_pred <- NULL
  if (length(te_idx) > 0) {
    pr <- predict(final_model, data = data.frame(Xall[te_idx, , drop=FALSE]))$predictions
    p1 <- if (is.matrix(pr)) pr[,2] else pr
    test_curves <- RFx$plot_rocpr(y01[te_idx], p1,
                                  out_png_prefix = file.path(out_fig,"TEST"),
                                  title_prefix = "TEST")
    test_pred <- p1
  }
  
  ## ---------------- importance via SHAP ----------------
  shap_out <- RFx$shap_table(final_model, Xall[tr_idx, , drop=FALSE], ps = psIN, topN = 20)
  write.csv(shap_out$table, file.path(out_tab, "importance_shap_feature_table.csv"), row.names = FALSE)
  RFx$plot_bar(shap_out$top, "Top features by mean |SHAP| (OOF train)", file.path(out_fig, "importance_feature.png"))
  
  # 택사 단위 집계
  agg_tax <- shap_out$table %>%
    group_by(TaxLabel) %>%
    summarise(MeanAbsShap = sum(MeanAbsShap), MeanShap = sum(MeanShap), .groups="drop") %>%
    arrange(desc(MeanAbsShap)) %>%
    mutate(Direction = ifelse(MeanShap >= 0, "Positive", "Negative"))
  write.csv(agg_tax, file.path(out_tab, "importance_shap_taxon_table.csv"), row.names = FALSE)
  RFx$plot_bar(head(agg_tax, 20), "Top taxa by mean |SHAP| (OOF train)", file.path(out_fig, "importance_taxon.png"))
  
  ## ---------------- save artifacts ----------------
  # 예측 저장
  oof_tab <- data.frame(
    Sample = rownames(Xall)[tr_idx],
    UID    = uid_vec[tr_idx],
    y      = y01[tr_idx],
    pred   = OOF_pred_best
  )
  fwrite(oof_tab, file.path(out_tab, "oof_pred.csv"))
  
  if (!is.null(test_pred)) {
    te_tab <- data.frame(
      Sample = rownames(Xall)[te_idx],
      UID    = uid_vec[te_idx],
      y      = y01[te_idx],
      pred   = test_pred
    )
    fwrite(te_tab, file.path(out_tab, "test_pred.csv"))
  }
  
  # 메타 저장(재현성용)
  rf_meta <- list(
    seed = seed,
    date = Sys.time(),
    project = project,
    uid_col = uid_col,
    outcome = outcome,
    testSet = testSet,
    K = K,
    prev_thr = prev_thr,
    mean_abund_thr = mean_abund_thr,
    pseudo = pseudo,
    n_random = n_random,
    num.trees = num.trees,
    best_param = best_param,
    oof_metrics = oof_curves,
    test_metrics = test_curves,
    class.weights = class.weights,
    y_levels = ylv,
    clinical = list(
      used = length(clinical_vari) > 0,
      variables = clinical_vari,
      scaler = clin$scaler,
      levels_map = clin$levels_map,
      cols = clin$cols
    ),
    features = colnames(Xall),
    train_index = tr_idx,
    test_index  = te_idx
  )
  saveRDS(rf_meta, file.path(out_rds, "rf_meta.rds"))
  saveRDS(final_model, file.path(out_rds, "rf_final.rds"))
  
  cat("\n[SAVED] artifacts under:\n",
      " - ", normalizePath(out_dir), "\n",
      "   * rds/rf_final.rds, rds/rf_meta.rds\n",
      "   * fig/OOF_ROC.png, fig/OOF_PR.png, (TEST_*.png if testSet=TRUE)\n",
      "   * fig/importance_feature.png, fig/importance_taxon.png\n",
      "   * table/oof_pred.csv, table/test_pred.csv (if test)\n",
      "   * table/importance_shap_feature_table.csv, table/importance_shap_taxon_table.csv\n", sep = "")
  
  invisible(list(
    model = final_model,
    meta = rf_meta,
    oof = oof_tab,
    test = if (!is.null(test_pred)) te_tab else NULL,
    shap = shap_out,
    imp_taxon = agg_tax
  ))
}