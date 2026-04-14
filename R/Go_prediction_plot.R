#' Go_prediction_plot
#'
#' Standardized plotting helper for prediction results saved by
#' \code{Go_prediction()}.
#'
#' @description
#' Reads a prediction result directory, detects whether it contains CV-only
#' out-of-fold predictions or Train/Test holdout predictions, and writes
#' standardized PDF figures under
#' \verb{<project_YYMMDD>/pdf/ROC/},
#' \verb{<project_YYMMDD>/pdf/PRROC/}, and
#' \verb{<project_YYMMDD>/pdf/prediction_Importance/}. When multiple result
#' directories are provided, the function creates combined comparison plots for
#' OOF predictions only.
#'
#' @param result Character output directory returned by \code{Go_prediction()}.
#'   For backward compatibility, a list with an \code{outdir} element is also
#'   accepted. If multiple result directories are provided, the function creates
#'   combined OOF ROC/PRROC comparison plots.
#' @param topN Integer; number of top importance features to display.
#' @param name Optional character label appended to output filenames.
#' @param compare_cols Optional color vector for multi-result OOF comparison
#'   plots.
#' @param mycols Named character vector for importance direction colors.
#'   Defaults to \code{c(Positive = "#d62728", Negative = "#1f77b4",
#'   Neutral = "grey70")}.
#' @param roc_width Numeric PDF width for ROC.
#' @param roc_height Numeric PDF height for ROC.
#' @param pr_width Numeric PDF width for PR.
#' @param pr_height Numeric PDF height for PR.
#' @param imp_width Numeric PDF width for importance.
#'
#' @return Saves PDF files and returns \code{NULL} invisibly.
#'
#' @details
#' With a single \code{result}, the function creates ROC/PRROC PDFs for either
#' OOF predictions (\code{testSet = FALSE}) or Train/Test predictions
#' (\code{testSet = TRUE}), plus an importance PDF when
#' \code{importance_feature.csv} is present.
#'
#' With multiple \code{result} directories, the function switches to comparison
#' mode and draws combined ROC/PRROC PDFs using \code{set = "OOF"} only.
#' Train/Test comparison is intentionally not supported in this mode.
#'
#' Output filenames include the optional \code{name} tag and the current date.
#'
#' @examples
#' \dontrun{
#' outdir <- Go_prediction(...)
#' Go_prediction_plot(outdir, name = "single_run")
#'
#' Go_prediction_plot(
#'   result = c(
#'     "Project_250101/no_test_bac_LightGBM",
#'     "Project_250101/no_test_clin_LightGBM",
#'     "Project_250101/no_test_all_LightGBM"
#'   ),
#'   name = "OOF_compare"
#' )
#' }
#'
#' @param patchwork Logical. If \code{TRUE}, skip saving and return the plot object(s) for use with \code{Gg_patchwork()} or the \pkg{patchwork} package. Default \code{FALSE}.
#' @export
Go_prediction_plot <- function(result,
                               topN = 20,
                               name = NULL,
                               compare_cols = NULL,
                               mycols = c(Positive = "#d62728",
                                          Negative = "#1f77b4",
                                          Neutral = "grey70"),
                               roc_width = 4,
                               roc_height = 4,
                               pr_width = 4,
                               pr_height = 4.5,
                               imp_width = 7,
                               patchwork = FALSE) {
  needed <- c("pROC", "PRROC", "ggplot2")
  missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    stop(
      "[Go_prediction_plot] Required package(s) not installed: ",
      paste(missing, collapse = ", "),
      ".\nRun Gotool_dependency() first and retry."
    )
  }

  resolve_outdirs <- function(x) {
    if (is.character(x) && length(x) >= 1 && all(nzchar(x))) {
      return(unname(x))
    }
    if (is.list(x)) {
      if (!is.null(x$outdir) && length(x$outdir) == 1) {
        return(x$outdir)
      }
      vals <- vapply(x, function(el) {
        if (is.character(el) && length(el) == 1 && nzchar(el)) return(el)
        if (is.list(el) && !is.null(el$outdir) && length(el$outdir) == 1) return(el$outdir)
        NA_character_
      }, character(1))
      vals <- vals[!is.na(vals)]
      if (length(vals)) return(unname(vals))
    }
    stop("[Error] result must be a Go_prediction() output directory, a character vector of directories, or a list containing $outdir.")
  }

  detect_meta_file <- function(outdir) {
    hits <- c("rf_meta.rds", "xgb_meta.rds", "lgb_meta.rds")
    hit_path <- file.path(outdir, hits)
    hit_path[file.exists(hit_path)][1]
  }

  format_auc_label <- function(x) {
    ifelse(is.na(x), "NA", sprintf("%.3f", x))
  }

  clean_tag <- function(x) {
    x <- gsub("[^A-Za-z0-9._-]+", "_", x)
    x <- gsub("_+", "_", x)
    x <- gsub("^_|_$", "", x)
    x
  }

  default_compare_cols <- function(n) {
    base_cols <- c("#1f77b4", "#2ca02c", "#d62728", "#ff7f0e", "#9467bd", "#8c564b")
    rep(base_cols, length.out = n)
  }

  derive_compare_labels <- function(paths) {
    bases <- basename(paths)
    split_tokens <- strsplit(bases, "_", fixed = TRUE)
    min_len <- min(lengths(split_tokens))
    common_prefix_len <- 0

    if (min_len > 0) {
      for (i in seq_len(min_len)) {
        tok_i <- vapply(split_tokens, `[`, character(1), i)
        if (length(unique(tok_i)) == 1) {
          common_prefix_len <- i
        } else {
          break
        }
      }
    }

    labels <- vapply(split_tokens, function(tok) {
      rest <- tok[seq.int(common_prefix_len + 1, length(tok))]
      rest <- rest[nzchar(rest)]
      if (!length(rest)) {
        return(paste(tok, collapse = "_"))
      }
      paste(rest, collapse = "_")
    }, character(1))

    labels <- clean_tag(labels)
    labels[!nzchar(labels)] <- clean_tag(bases[!nzchar(labels)])
    labels
  }

  compute_pr <- function(prob, outcome, positive_class) {
    PRROC::pr.curve(
      scores.class0 = prob[outcome == positive_class],
      scores.class1 = prob[outcome != positive_class],
      curve = TRUE
    )
  }

  build_model_label <- function(meta_path, outdir) {
    base_nm <- basename(outdir)
    if (grepl("^Randomforest", base_nm, ignore.case = TRUE)) return("Random Forest")
    if (grepl("^Xgboost", base_nm, ignore.case = TRUE)) return("XGBoost")
    if (grepl("^LightGBM", base_nm, ignore.case = TRUE)) return("LightGBM")
    if (grepl("rf_meta\\.rds$", meta_path)) return("Random Forest")
    if (grepl("xgb_meta\\.rds$", meta_path)) return("XGBoost")
    if (grepl("lgb_meta\\.rds$", meta_path)) return("LightGBM")
    base_nm
  }

  if (isTRUE(patchwork)) {
    message("[Go_prediction_plot] patchwork = TRUE: ROC/PR plots are base-R graphics (skipped). Only the feature importance ggplot is returned.")
  }

  outdirs <- resolve_outdirs(result)
  if (!length(outdirs)) {
    stop("[Error] no valid result directory found.")
  }

  if (length(outdirs) > 1) {
    if (any(!dir.exists(outdirs))) {
      stop(sprintf("[Error] result directory not found: %s", paste(outdirs[!dir.exists(outdirs)], collapse = ", ")))
    }

    out_roots <- unique(dirname(outdirs))
    if (length(out_roots) != 1) {
      stop("[Error] multi-result comparison requires all result directories to be under the same project directory.")
    }

    comp_labels <- derive_compare_labels(outdirs)
    if (anyDuplicated(comp_labels)) {
      comp_labels <- make.unique(comp_labels, sep = "_")
    }

    comp_cols <- compare_cols %||% default_compare_cols(length(outdirs))
    if (length(comp_cols) < length(outdirs)) {
      comp_cols <- rep(comp_cols, length.out = length(outdirs))
    }
    names(comp_cols) <- comp_labels

    roc_dir <- file.path(out_roots, "pdf", "ROC")
    pr_dir <- file.path(out_roots, "pdf", "PRROC")
    dir.create(roc_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(pr_dir, recursive = TRUE, showWarnings = FALSE)

    roc_list <- vector("list", length(outdirs))
    pr_list <- vector("list", length(outdirs))
    auc_list <- numeric(length(outdirs))
    auprc_list <- numeric(length(outdirs))
    model_labels <- character(length(outdirs))
    levels_ref <- NULL
    pos_ref <- NULL

    for (i in seq_along(outdirs)) {
      outdir_i <- outdirs[i]
      meta_path_i <- detect_meta_file(outdir_i)
      if (length(meta_path_i) == 0 || !file.exists(meta_path_i)) {
        stop(sprintf("[Error] meta file not found in: %s", outdir_i))
      }
      pred_path_i <- file.path(outdir_i, "predictions.csv")
      if (!file.exists(pred_path_i)) {
        stop(sprintf("[Error] predictions.csv not found in: %s", outdir_i))
      }

      meta_i <- readRDS(meta_path_i)
      preds_i <- utils::read.csv(pred_path_i, stringsAsFactors = FALSE)
      if (!all(c("outcome", "pred", "set") %in% colnames(preds_i))) {
        stop(sprintf("[Error] predictions.csv in %s must contain outcome, pred, and set.", outdir_i))
      }
      if (!any(preds_i$set == "OOF")) {
        avail_sets <- paste(sort(unique(preds_i$set)), collapse = ", ")
        stop(sprintf(
          "[Error] Go_prediction_plot comparison mode supports OOF only. '%s' does not contain set='OOF' (found: %s). Use a single result directory for Train/Test plots.",
          outdir_i, avail_sets
        ))
      }

      levels_i <- meta_i$levels %||% unique(preds_i$outcome)
      pos_i <- meta_i$positive_class %||% levels_i[2]
      if (is.null(levels_ref)) {
        levels_ref <- levels_i
        pos_ref <- pos_i
      } else {
        if (!identical(as.character(levels_i), as.character(levels_ref))) {
          stop("[Error] all comparison results must share the same outcome level order.")
        }
        if (!identical(as.character(pos_i), as.character(pos_ref))) {
          stop("[Error] all comparison results must share the same positive class.")
        }
      }

      oo_i <- preds_i[preds_i$set == "OOF", , drop = FALSE]
      oo_i$outcome <- factor(oo_i$outcome, levels = levels_ref)
      roc_i <- pROC::roc(oo_i$outcome, oo_i$pred, levels = levels_ref, direction = "<", quiet = TRUE)
      pr_i <- compute_pr(oo_i$pred, oo_i$outcome, pos_ref)

      roc_list[[i]] <- roc_i
      pr_list[[i]] <- pr_i
      auc_list[i] <- as.numeric(roc_i$auc)
      auprc_list[i] <- pr_i$auc.integral
      model_labels[i] <- build_model_label(meta_path_i, outdir_i)
    }

    names(roc_list) <- comp_labels
    names(pr_list) <- comp_labels
    names(auc_list) <- comp_labels
    names(auprc_list) <- comp_labels

    date_tag <- format(Sys.Date(), "%y%m%d")
    name_tag <- if (!is.null(name) && nzchar(name)) clean_tag(name) else "comparison"
    file_tag <- paste(c("prediction", name_tag, date_tag), collapse = "_")
    title_model <- paste(unique(model_labels), collapse = " / ")

    first_name <- comp_labels[1]
    if (!isTRUE(patchwork)) {
      grDevices::pdf(file.path(roc_dir, sprintf("%s_ROC.pdf", file_tag)),
                     width = roc_width + 1, height = roc_height + 1)
      pROC::plot.roc(roc_list[[first_name]],
                     col = comp_cols[first_name],
                     lwd = 2,
                     asp = 1,
                     main = sprintf("OOF ROC comparison (%s)", title_model))
      if (length(roc_list) > 1) {
        for (nm in comp_labels[-1]) {
          pROC::plot.roc(roc_list[[nm]], col = comp_cols[nm], lwd = 2, add = TRUE)
        }
      }
      graphics::abline(0, 1, lty = 2, col = "grey70")
      graphics::legend("bottomright",
                       legend = sprintf("%s (AUC=%s)", comp_labels, format_auc_label(auc_list)),
                       col = comp_cols[comp_labels], lwd = 2, bty = "n")
      grDevices::dev.off()

      grDevices::pdf(file.path(pr_dir, sprintf("%s_PRROC.pdf", file_tag)),
                     width = pr_width + 1, height = pr_height + 0.5)
      graphics::plot(pr_list[[first_name]]$curve[, 1], pr_list[[first_name]]$curve[, 2],
                     type = "l", lwd = 2, col = comp_cols[first_name],
                     xlab = "Recall", ylab = "Precision",
                     main = sprintf("OOF PR comparison (%s)", title_model))
      if (length(pr_list) > 1) {
        for (nm in comp_labels[-1]) {
          graphics::lines(pr_list[[nm]]$curve[, 1], pr_list[[nm]]$curve[, 2],
                          lwd = 2, col = comp_cols[nm])
        }
      }
      graphics::legend("bottomleft",
                       legend = sprintf("%s (AUPRC=%s)", comp_labels, format_auc_label(auprc_list)),
                       col = comp_cols[comp_labels], lwd = 2, bty = "n")
      grDevices::dev.off()
    }

    return(invisible(NULL))  # comparison mode: no ggplot to return
  }

  outdir <- outdirs[1]
  if (!dir.exists(outdir)) {
    stop(sprintf("[Error] result directory not found: %s", outdir))
  }

  meta_path <- detect_meta_file(outdir)
  if (length(meta_path) == 0 || !file.exists(meta_path)) {
    stop(sprintf("[Error] meta file not found in: %s", outdir))
  }
  pred_path <- file.path(outdir, "predictions.csv")
  imp_path <- file.path(outdir, "importance_feature.csv")
  if (!file.exists(pred_path)) {
    stop(sprintf("[Error] predictions.csv not found in: %s", outdir))
  }

  meta <- readRDS(meta_path)
  preds <- utils::read.csv(pred_path, stringsAsFactors = FALSE)
  if (!all(c("outcome", "pred") %in% colnames(preds))) {
    stop("[Error] predictions.csv must contain outcome and pred columns.")
  }

  out_root <- dirname(outdir)
  pdf_root <- file.path(out_root, "pdf")
  roc_dir <- file.path(pdf_root, "ROC")
  pr_dir <- file.path(pdf_root, "PRROC")
  imp_dir <- file.path(pdf_root, "prediction_Importance")
  dir.create(roc_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(pr_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(imp_dir, recursive = TRUE, showWarnings = FALSE)

  levels_out <- meta$levels %||% unique(preds$outcome)
  if (length(levels_out) != 2) {
    stop("[Error] outcome levels must have length 2.")
  }
  positive_class <- meta$positive_class %||% levels_out[2]
  preds$outcome <- factor(preds$outcome, levels = levels_out)
  model_label <- build_model_label(meta_path, outdir)
  file_stub <- basename(outdir)
  date_tag <- format(Sys.Date(), "%y%m%d")
  name_tag <- if (!is.null(name) && nzchar(name)) clean_tag(name) else NULL
  file_tag <- paste(c("prediction", name_tag, file_stub, date_tag), collapse = "_")

  has_set <- "set" %in% colnames(preds)
  is_holdout <- has_set && any(preds$set %in% c("Train", "Test"))
  is_oof <- has_set && any(preds$set == "OOF")

  if (!is_holdout && !is_oof) {
    stop("[Error] predictions.csv must contain set='OOF' or set in c('Train','Test').")
  }

  if (is_holdout) {
    tr <- preds[preds$set == "Train", , drop = FALSE]
    te <- preds[preds$set == "Test", , drop = FALSE]

    roc_tr <- pROC::roc(tr$outcome, tr$pred, levels = levels_out, direction = "<", quiet = TRUE)
    roc_te <- pROC::roc(te$outcome, te$pred, levels = levels_out, direction = "<", quiet = TRUE)

    grDevices::pdf(file.path(roc_dir, sprintf("%s_ROC.pdf", file_tag)),
                   width = roc_width, height = roc_height)
    pROC::plot.roc(roc_tr, col = "#1f77b4", lwd = 2, asp = 1,
                   main = sprintf("%s ROC", model_label))
    pROC::plot.roc(roc_te, col = "#d62728", lwd = 2, add = TRUE, asp = 1)
    graphics::abline(0, 1, lty = 2, col = "grey70")
    graphics::legend("bottomright",
                     legend = c(sprintf("Train (AUC=%s)", format_auc_label(as.numeric(roc_tr$auc))),
                                sprintf("Test (AUC=%s)", format_auc_label(as.numeric(roc_te$auc)))),
                     col = c("#1f77b4", "#d62728"), lwd = 2, bty = "n")
    grDevices::dev.off()

    pr_tr <- compute_pr(tr$pred, tr$outcome, positive_class)
    pr_te <- compute_pr(te$pred, te$outcome, positive_class)
    grDevices::pdf(file.path(pr_dir, sprintf("%s_PRROC.pdf", file_tag)),
                   width = pr_width, height = pr_height)
    graphics::plot(pr_tr$curve[, 1], pr_tr$curve[, 2], type = "l", lwd = 2, col = "#1f77b4",
                   xlab = "Recall", ylab = "Precision",
                   main = sprintf("%s PR", model_label))
    graphics::lines(pr_te$curve[, 1], pr_te$curve[, 2], lwd = 2, col = "#d62728")
    graphics::legend("bottomleft",
                     legend = c(sprintf("Train (AUPRC=%s)", format_auc_label(pr_tr$auc.integral)),
                                sprintf("Test (AUPRC=%s)", format_auc_label(pr_te$auc.integral))),
                     col = c("#1f77b4", "#d62728"), lwd = 2, bty = "n")
    grDevices::dev.off()
  } else {
    oo <- preds[preds$set == "OOF", , drop = FALSE]
    roc_oo <- pROC::roc(oo$outcome, oo$pred, levels = levels_out, direction = "<", quiet = TRUE)
    grDevices::pdf(file.path(roc_dir, sprintf("%s_ROC.pdf", file_tag)),
                   width = roc_width, height = roc_height)
    pROC::plot.roc(roc_oo, col = "#1f77b4", lwd = 2, asp = 1,
                   main = sprintf("%s OOF ROC", model_label))
    graphics::abline(0, 1, lty = 2, col = "grey70")
    graphics::legend("bottomright",
                     legend = sprintf("OOF (AUC=%s)", format_auc_label(as.numeric(roc_oo$auc))),
                     col = "#1f77b4", lwd = 2, bty = "n")
    grDevices::dev.off()

    pr_oo <- compute_pr(oo$pred, oo$outcome, positive_class)
    grDevices::pdf(file.path(pr_dir, sprintf("%s_PRROC.pdf", file_tag)),
                   width = pr_width, height = pr_height)
    graphics::plot(pr_oo$curve[, 1], pr_oo$curve[, 2], type = "l", lwd = 2, col = "#1f77b4",
                   xlab = "Recall", ylab = "Precision",
                   main = sprintf("%s OOF PR", model_label))
    graphics::legend("bottomleft",
                     legend = sprintf("OOF (AUPRC=%s)", format_auc_label(pr_oo$auc.integral)),
                     col = "#1f77b4", lwd = 2, bty = "n")
    grDevices::dev.off()
  }

  if (file.exists(imp_path)) {
    imp_feat <- utils::read.csv(imp_path, stringsAsFactors = FALSE)
    score_col <- if ("PermImportance" %in% colnames(imp_feat)) "PermImportance" else {
      num_cols <- names(Filter(is.numeric, imp_feat))
      if (!length(num_cols)) NA_character_ else num_cols[1]
    }
    label_col <- if ("Taxon" %in% colnames(imp_feat)) "Taxon" else {
      if ("Feature" %in% colnames(imp_feat)) "Feature" else NA_character_
    }
    if (!is.na(score_col) && !is.na(label_col) && nrow(imp_feat) > 0) {
      plot_df <- head(imp_feat[order(-imp_feat[[score_col]]), , drop = FALSE], topN)
      plot_df[[label_col]] <- as.character(plot_df[[label_col]])
      plot_df$Direction <- as.character(plot_df$Direction %||% rep("Neutral", nrow(plot_df)))
      plot_df$Direction[is.na(plot_df$Direction)] <- "Neutral"

      max_lbl <- max(nchar(plot_df[[label_col]]), na.rm = TRUE)
      imp_height <- max(4, min(12, 1.5 + 0.22 * nrow(plot_df) + 0.015 * max_lbl))

      p_imp <- ggplot2::ggplot(
        plot_df,
        ggplot2::aes_string(x = sprintf("reorder(%s, %s)", label_col, score_col),
                            y = score_col,
                            fill = "Direction")
      ) +
        ggplot2::geom_col() +
        ggplot2::coord_flip() +
        ggplot2::labs(title = sprintf("Top %d Features (%s)", nrow(plot_df), model_label),
                      x = NULL, y = "Permutation Importance") +
        ggplot2::scale_fill_manual(values = mycols, drop = FALSE) +
        ggplot2::theme_bw(base_size = 12) +
        ggplot2::theme(
          text = ggplot2::element_text(size = 11),
          plot.title = ggplot2::element_text(size = 11, hjust = 1),
          axis.text.y = ggplot2::element_text(size = 11, face = "italic"),
          legend.title = ggplot2::element_blank()
        )

        ggplot2::ggsave(
        filename = file.path(imp_dir, sprintf("%s_Importance.pdf", file_tag)),
        plot = p_imp, device = "pdf", width = imp_width, height = imp_height, limitsize = FALSE
      )
    }
  }

  invisible(NULL)
}
