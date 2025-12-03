#' Go_Maaslin2
#'
#' Wrapper around \code{Maaslin2::Maaslin2()} for differential abundance analysis
#' from a \code{phyloseq} object. Handles output foldering/naming, optional
#' relative→integer conversion, ASV→taxon relabeling, factor level ordering,
#' and passes fixed/random effects to MaAsLin2. Supports both a single global
#' multilevel model and pairwise / k-level subset modes.
#'
#' @param psIN A \code{phyloseq} object with OTU/ASV table and sample metadata.
#' @param project Character. Project prefix for dated output folders/files
#'   (e.g., \code{"MyProj"} → \code{"MyProj_250818/"}).
#' @param fixed_effects Character vector of metadata column names used as fixed
#'   effects. Currently the first element (\code{fixed_effects[1]}) is used as
#'   the primary factor for level ordering and combination logic.
#' @param random_effects Optional character vector of metadata column names used
#'   as random effects (not all analyses require this). Default \code{NULL}.
#' @param normalization Character scalar passed to MaAsLin2. One of
#'   \code{"TSS"}, \code{"CLR"}, \code{"CSS"}, \code{"NONE"}. Default \code{"TSS"}.
#' @param transform Logical. If \code{TRUE} and the input looks like relative
#'   abundance, multiply by median library size and round to integers before
#'   running MaAsLin2. Default \code{TRUE}.
#' @param orders Optional character vector specifying the desired level order
#'   for the primary fixed effect (\code{fixed_effects[1]}), e.g.
#'   \code{c("Pre","wk1_Post","mo1_Post")}. Levels not present in the data are
#'   ignored. If \code{NULL}, the existing factor order is kept.
#' @param out_dir Optional output directory to write MaAsLin2 results. If
#'   \code{NULL}, a dated project path is created under
#'   \code{<project_YYMMDD>/table/MaAsLin2/}.
#' @param name Optional character suffix to distinguish runs (e.g.,
#'   \code{"AdjBMI"}). Used in output subdirectory naming
#'   (\code{MaAsLin2.<name>.(FE=... ).(RE=...) }).
#' @param data Optional character flag for relabeling features. If set to
#'   \code{"ASV"}, the function attempts to replace OTU/ASV IDs by taxonomic
#'   labels (Species → Genus_Species → Genus → higher ranks) using
#'   \code{tax_table(psIN)} and \code{make.unique()} to ensure uniqueness.
#'   Any other value or \code{NULL} keeps the original rownames.
#' @param combination Numeric or \code{NULL}. Controls how the primary fixed
#'   effect levels are compared when \code{global = FALSE}:
#'   \itemize{
#'     \item \code{NULL}: Single model using all available levels (no subsetting).
#'     \item \code{2}: Run MaAsLin2 for all pairwise combinations of levels
#'       (A–B, A–C, B–C, ...). Each pair is analyzed as a two-level factor with
#'       a separate output folder.
#'     \item \code{>2}: For a given \code{k}, run MaAsLin2 on all k-level subsets
#'       of the factor levels (e.g., A+B+C, A+B+D, ...), skipping subsets where
#'       the per-level sample size is too small.
#'   }
#'   Ignored when \code{global = TRUE}.
#' @param min_per_level Integer. Minimum number of samples required per level
#'   of the primary fixed effect when running pairwise or k-level subset
#'   analyses (\code{combination = 2} or \code{> 2}). Subsets not meeting this
#'   requirement are skipped with a message. Default \code{3}.
#' @param global Logical. If \code{TRUE}, run a single global multilevel model
#'   including all levels of \code{fixed_effects[1]} in one MaAsLin2 call
#'   (omnibus-style). In this mode, \code{combination} is ignored. If
#'   \code{FALSE}, the behavior is controlled by \code{combination} as described
#'   above. Default \code{FALSE}.
#' @param max_sig Numeric scalar passed to MaAsLin2 as \code{max_significance}.
#'   Controls the q-value threshold used for selecting features to plot in
#'   MaAsLin2's heatmap and scatter plots (and for defining "significant"
#'   features). The default \code{0.25} corresponds to MaAsLin2's original
#'   default. Increasing this value (e.g., \code{0.4}–\code{0.5}) can be useful
#'   when you want plots even when few features pass a strict FDR cutoff.
#' @param seed Integer random seed for reproducibility. Passed to
#'   \code{set.seed()} inside each MaAsLin2 run. Default \code{123}.
#'
#' @details
#' \strong{Relative→integer transform:}
#' The helper heuristically detects relative abundance by the mean of
#' \code{sample_sums(psIN)} (≈100 implies relative). If detected and
#' \code{transform=TRUE}, counts are reconstructed by multiplying the feature
#' table by the median library size and rounding. Set \code{transform=FALSE}
#' to skip this step and use the original counts.
#'
#' \strong{Global multilevel vs. combination modes:}
#' When \code{global = TRUE}, the function runs exactly one MaAsLin2 model
#' using all samples and all levels of \code{fixed_effects[1]} together. This
#' preserves MaAsLin2's multilevel global (omnibus) testing behavior and is
#' typically recommended when you are interested in an overall association
#' across multiple groups (e.g., CP class A/B/C, multiple timepoints, etc.).
#'
#' When \code{global = FALSE}, two types of analyses are supported:
#' \itemize{
#'   \item Single-model mode (\code{combination = NULL}): one MaAsLin2 run
#'     including all levels, similar to the global mode but without any special
#'     handling.
#'   \item Pairwise / k-level subset mode (\code{combination = 2} or \code{> 2}):
#'     MaAsLin2 is run on each level pair or subset separately. This is useful
#'     for detailed pairwise contrasts, but does not provide a single omnibus
#'     global p-value across all levels.
#' }
#'
#' \strong{Ordering factor levels:}
#' If \code{orders} is supplied, the primary fixed effect
#' (\code{fixed_effects[1]}) is refactored to the specified level order. Levels
#' not present in the data are dropped. This affects both the modeling and the
#' interpretation of coefficients/contrasts.
#'
#' Outputs are written under a dated project directory; MaAsLin2’s heatmap and
#' scatter plots are enabled by default, with plotting significance controlled
#' via \code{max_sig}.
#'
#' @return (Invisibly) the \code{Maaslin2} result object (list) from the last
#'   MaAsLin2 run. Side effects: writes tables (TSV+CSV), settings, and plots
#'   to \code{out_dir} and its subdirectories.
#'
#' @seealso \code{\link[Maaslin2]{Maaslin2}}
#'
#' @examples
#' \dontrun{
#' # Basic global multilevel model (recommended when multiple groups)
#' res_global <- Go_Maaslin2(
#'   psIN          = ps,
#'   project       = "IBD",
#'   fixed_effects = c("Group"),
#'   random_effects = c("Subject"),
#'   normalization = "TSS",
#'   transform     = TRUE,
#'   orders        = c("Control","Case"),
#'   name          = "GlobalModel",
#'   global        = TRUE
#' )
#'
#' # Pairwise comparisons between all Group levels, with a relaxed plotting cutoff
#' res_pairwise <- Go_Maaslin2(
#'   psIN          = ps,
#'   project       = "IBD",
#'   fixed_effects = c("Group"),
#'   random_effects = c("Subject"),
#'   normalization = "TSS",
#'   transform     = TRUE,
#'   orders        = c("Control","Case"),
#'   name          = "Pairwise",
#'   combination   = 2,
#'   global        = FALSE,
#'   max_sig       = 0.5
#' )
#' }
#'
#' @importFrom phyloseq sample_data taxa_are_rows otu_table sample_sums tax_table
#' @importFrom stats median
#' @importFrom utils file_test
#' @export

Go_Maaslin2 <- function(psIN,
                        project,
                        fixed_effects,
                        random_effects = NULL,
                        normalization = "TSS",
                        transform = TRUE,
                        orders = NULL,           # c("wk1_Post","mo1_Post",...)
                        out_dir = NULL,
                        name = NULL,             # 상위 태그 (예: "Case"/"Control")
                        data = NULL,             # "ASV"면 taxa 기반 라벨 치환
                        combination = NULL,      # 2 => 모든 페어, >2 => 그 레벨들만 서브셋
                        min_per_level = 3,       # 레벨별 최소 샘플수 안전장치
                        global = FALSE,
                        max_sig = 0.25, # NEW: TRUE면 full multilevel model 1번만 실행
                        seed = 123)              # 재현성
{
  ## ---------------- helpers ----------------
  .safe_tag <- function(x){
    x <- gsub("[^A-Za-z0-9+._-]", "_", x); gsub("__+", "_", x)
  }
  .tag_FE <- function(fx) sprintf("(FE=%s)", .safe_tag(paste(fx, collapse="+")))
  .tag_RE <- function(re) if (is.null(re) || length(re)==0) "(RE=None)" else sprintf("(RE=%s)", .safe_tag(paste(re, collapse="+")))

  .write_settings <- function(dir_path, meta, args){
    fn <- file.path(dir_path, "settings.txt")
    lines <- c(
      sprintf("date: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
      sprintf("n_samples: %d", nrow(meta)),
      sprintf("fixed_effects: %s", paste(args$fixed_effects, collapse = " + ")),
      sprintf("random_effects: %s", ifelse(is.null(args$random_effects),"None", paste(args$random_effects, collapse = " + "))),
      sprintf("normalization: %s", args$normalization),
      sprintf("transform(LOG inside MaAsLin2): TRUE"),
      sprintf("integerize_relative_by_median_depth: %s", isTRUE(args$transform)),
      sprintf("orders: %s", ifelse(is.null(args$orders),"NULL", paste(args$orders, collapse = ", "))),
      sprintf("data_label_mode: %s", ifelse(is.null(args$data),"ASV_id", toupper(args$data))),
      sprintf("combination: %s", ifelse(is.null(args$combination),"NULL", as.character(args$combination))),
      sprintf("global_full_model: %s", ifelse(isTRUE(args$global), "TRUE", "FALSE")),
      sprintf("min_per_level: %s", args$min_per_level),
      sprintf("seed: %s", args$seed)
    )
    writeLines(lines, fn, useBytes = TRUE)
  }

  .save_tsv_as_csv <- function(out_dir_run){
    for (f in c("all_results.tsv","significant_results.tsv")) {
      fp <- file.path(out_dir_run, f)
      if (file.exists(fp)) {
        df <- try(read.table(fp, header = TRUE, sep = "\t", check.names = FALSE), silent = TRUE)
        if (!inherits(df, "try-error")) {
          utils::write.csv(df, file = file.path(out_dir_run, sub("\\.tsv$", ".csv", f)), row.names = FALSE)
        }
      }
    }
  }

  ## ---------------- packages ----------------
  if (!requireNamespace("Maaslin2", quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager", repos = "https://cloud.r-project.org")
    BiocManager::install("Maaslin2", ask = FALSE, update = FALSE)
  }
  suppressPackageStartupMessages(library(Maaslin2))

  ## ---------------- outputs ----------------
  out_root  <- sprintf("%s_%s", project, format(Sys.Date(), "%y%m%d"))
  if (!file_test("-d", out_root)) dir.create(out_root)
  out_table <- file.path(out_root, "table"); if (!file_test("-d", out_table)) dir.create(out_table)
  out_DA    <- file.path(out_table, "MaAsLin2"); if (!file_test("-d", out_DA)) dir.create(out_DA, recursive = TRUE)

  FE_tag <- .tag_FE(fixed_effects)
  RE_tag <- .tag_RE(random_effects)
  base_tag <- if (is.null(name) || !nzchar(name)) "MaAsLin2.Base" else sprintf("MaAsLin2.%s", .safe_tag(name))
  subdir <- sprintf("%s.%s.%s", base_tag, FE_tag, RE_tag)

  out_dir_default <- file.path(out_DA, subdir)
  if (is.null(out_dir)) out_dir <- out_dir_default
  if (!file_test("-d", out_dir)) dir.create(out_dir, recursive = TRUE)
  message(sprintf("[INFO] MaAsLin2 base output -> %s", out_dir))

  ## ---------------- data ----------------
  metadata_df <- as.data.frame(sample_data(psIN))
  otu_mat <- if (taxa_are_rows(psIN)) as.matrix(otu_table(psIN)) else t(as.matrix(otu_table(psIN)))
  otu_mat <- as.data.frame(otu_mat)

  # ASV → taxa 라벨 치환 (원하면)
  if (!is.null(data) && toupper(data) == "ASV" && !is.null(tax_table(psIN, errorIfNULL = FALSE))) {
    tt <- as.data.frame(tax_table(psIN))
    label_vec <- if ("Species" %in% colnames(tt)) as.character(tt$Species) else rownames(tt)
    if (!("Species" %in% colnames(tt)) || all(is.na(label_vec) | label_vec == "")) {
      genus <- if ("Genus" %in% colnames(tt)) as.character(tt$Genus) else NA
      sp    <- if ("Species" %in% colnames(tt)) as.character(tt$Species) else NA
      label_vec <- ifelse(!is.na(genus) & nzchar(genus),
                          ifelse(!is.na(sp) & nzchar(sp), paste0(genus, "_", sp), genus),
                          NA)
    }
    if (all(is.na(label_vec) | label_vec == "")) {
      for (col in c("Genus","Family","Order","Class","Phylum","Kingdom")) {
        if (col %in% colnames(tt)) {
          fill <- as.character(tt[[col]])
          idx <- (is.na(label_vec) | label_vec == "")
          label_vec[idx] <- fill[idx]
        }
        if (!all(is.na(label_vec) | label_vec == "")) break
      }
    }
    old_ids <- rownames(otu_mat)
    label_vec[is.na(label_vec) | label_vec == ""] <- old_ids[is.na(label_vec) | label_vec == ""]
    rownames(otu_mat) <- make.unique(label_vec)
  }

  # 샘플 동기화
  common_samples <- intersect(colnames(otu_mat), rownames(metadata_df))
  if (length(common_samples) < 2) stop("[ERROR] Not enough overlapping samples between OTU and metadata.")
  otu_mat     <- otu_mat[, common_samples, drop = FALSE]
  metadata_df <- metadata_df[ common_samples, , drop = FALSE]

  # 상대 abundance면 정수화(선택)
  detect_abundance_type <- function(physeq) {
    lib_sizes <- sample_sums(physeq); mean_lib <- mean(lib_sizes)
    if (abs(mean_lib - 100) < 0.2) "relative" else if (mean_lib > 1000) "absolute" else "unknown"
  }
  if (detect_abundance_type(psIN) == "relative" && isTRUE(transform)) {
    total_reads <- median(sample_sums(psIN))
    message(sprintf("[INFO] Relative abundance detected. Multiply by median depth (%.0f) & round.", total_reads))
    otu_mat <- round(otu_mat * total_reads)
  } else {
    message("[INFO] Skip integer transform (absolute data or transform=FALSE).")
  }

  # 메타데이터 정리 + orders
  rn <- rownames(metadata_df)
  metadata_df <- as.data.frame(lapply(metadata_df, function(x){
    if (is.character(x)) {
      if (all(grepl("^[-+]?[0-9]*\\.?[0-9]+$", x[!is.na(x)]))) as.numeric(x) else as.factor(x)
    } else if (is.logical(x)) as.factor(x) else x
  }), stringsAsFactors = FALSE)
  rownames(metadata_df) <- rn

  if (length(fixed_effects) == 0) stop("[ERROR] fixed_effects must be provided.")
  fx <- fixed_effects[1]
  if (!fx %in% colnames(metadata_df)) stop("[ERROR] Missing fixed_effect in metadata: ", fx)
  if (!is.factor(metadata_df[[fx]])) metadata_df[[fx]] <- as.factor(metadata_df[[fx]])

  if (!is.null(orders) && length(orders) > 0 && is.vector(orders)) {
    keep <- intersect(orders, levels(metadata_df[[fx]]))
    if (length(keep) >= 2) metadata_df[[fx]] <- factor(metadata_df[[fx]], levels = keep)
  }
  if (!is.null(random_effects)) {
    for (re in random_effects) if (re %in% colnames(metadata_df) && !is.factor(metadata_df[[re]]))
      metadata_df[[re]] <- factor(metadata_df[[re]])
  }
  metadata_df <- metadata_df[colnames(otu_mat), , drop = FALSE]

  # 공통 러너
  .run_one <- function(otu_df, meta_df, fx, out_dir_run){
    if (!file_test("-d", out_dir_run)) dir.create(out_dir_run, recursive = TRUE)
    set.seed(seed)
    fit <- Maaslin2::Maaslin2(
      input_data     = otu_df,
      input_metadata = meta_df,
      output         = out_dir_run,
      fixed_effects  = fx,
      random_effects = random_effects,
      normalization  = normalization,
      transform      = "LOG",
      plot_heatmap   = TRUE,
      plot_scatter   = TRUE,
      max_significance = max_sig
    )
    # settings 기록용 인자 패키징
    args <- list(
      fixed_effects = fx,
      random_effects = random_effects,
      normalization = normalization,
      transform = transform,
      orders = orders,
      data = data,
      combination = combination,
      global = global,
      min_per_level = min_per_level,
      seed = seed
    )
    .write_settings(out_dir_run, meta_df, args)
    .save_tsv_as_csv(out_dir_run)
    invisible(fit)
  }

  ## ---------------- run: global / single / combination ----------------
  levs_all <- levels(metadata_df[[fx]])

  ## 1) Global full-model mode (multilevel factor, omnibus p-value)
  if (isTRUE(global)) {
    message("[INFO] Global mode enabled: running ONE full multilevel model (all levels) for omnibus test.")
    .run_one(otu_mat, metadata_df, fx, out_dir)
    message("[INFO] Global model finished. Omnibus p-values are in all_results.tsv / significant_results.tsv.")
    return(invisible(TRUE))
  }

  ## 2) 기존 모드들 (global = FALSE일 때만)
  if (is.null(combination)) {
    message("[INFO] Running single model (no combination).")
    .run_one(otu_mat, metadata_df, fx, out_dir)

  } else if (is.numeric(combination) && combination == 2) {
    message("[INFO] Pairwise mode: running ALL pairs.")
    # 모든 조합 (A–B, A–C, B–C, ...)
    pairs <- t(combn(levs_all, 2))
    for (i in seq_len(nrow(pairs))) {
      a <- pairs[i, 1]; b <- pairs[i, 2]
      keep_idx <- metadata_df[[fx]] %in% c(a, b)
      meta_sub <- droplevels(metadata_df[keep_idx, , drop = FALSE])
      otu_sub  <- otu_mat[, rownames(meta_sub), drop = FALSE]

      tbl <- table(meta_sub[[fx]])
      if (any(tbl < min_per_level)) {
        message(sprintf("[SKIP] Too few samples: %s vs %s (counts: %s)",
                        a, b, paste(names(tbl), tbl, sep="=", collapse=", ")))
        next
      }
      # 페어 수준에서 레벨 순서 고정
      meta_sub[[fx]] <- factor(meta_sub[[fx]], levels = c(a, b))

      pair_core <- sprintf("MaAsLin2.%s.%s_vs_%s", fx, .safe_tag(a), .safe_tag(b))
      out_dir_pair <- file.path(out_dir, sprintf("%s.%s", pair_core, RE_tag))  # RE=... 반영
      message(sprintf("[INFO] Running %s ...", basename(out_dir_pair)))
      .run_one(otu_sub, meta_sub, fx, out_dir_pair)
    }

  } else if (is.numeric(combination) && combination > 2) {
    message(sprintf("[INFO] k-level subset mode: combination = %d", combination))
    sets <- combn(levs_all, combination, simplify = FALSE)
    for (levset in sets) {
      keep_idx <- metadata_df[[fx]] %in% levset
      meta_sub <- droplevels(metadata_df[keep_idx, , drop = FALSE])
      otu_sub  <- otu_mat[, rownames(meta_sub), drop = FALSE]
      # 안전장치
      tbl <- table(meta_sub[[fx]])
      if (any(tbl < min_per_level)) {
        message(sprintf("[SKIP] Too few samples for levels: %s (counts: %s)",
                        paste(levset, collapse="+"),
                        paste(names(tbl), tbl, sep="=", collapse=", ")))
        next
      }
      meta_sub[[fx]] <- factor(meta_sub[[fx]], levels = levset)
      tag_set <- paste(.safe_tag(levset), collapse = "+")
      out_dir_set <- file.path(out_dir, sprintf("MaAsLin2.%s.%s.%s", fx, tag_set, RE_tag))
      message(sprintf("[INFO] Running %s ...", basename(out_dir_set)))
      .run_one(otu_sub, meta_sub, fx, out_dir_set)
    }

  } else {
    stop("[ERROR] 'combination' must be NULL, 2, or >2 (numeric).")
  }

  message("[INFO] Done.")
  invisible(TRUE)
}
