#' Go_Maaslin2
#'
#' Wrapper around \code{Maaslin2::Maaslin2()} for differential abundance analysis
#' from a \code{phyloseq} object. Handles output foldering/naming, optional
#' relative→integer conversion, factor level ordering, and passes fixed/random
#' effects to MaAsLin2.
#'
#' @param psIN A \code{phyloseq} object with OTU/ASV table and sample metadata.
#' @param project Character. Project prefix for dated output folders/files
#'   (e.g., \code{"MyProj"} → \code{"MyProj_250818/"}).
#' @param fixed_effects Character vector of metadata column names used as fixed effects.
#' @param random_effects Optional character vector of metadata column names used
#'   as random effects (not all analyses require this). Default \code{NULL}.
#' @param normalization Character scalar passed to MaAsLin2. One of
#'   \code{"TSS"}, \code{"CLR"}, \code{"CSS"}, \code{"NONE"}. Default \code{"TSS"}.
#' @param transform Logical. If \code{TRUE} and the input looks like relative
#'   abundance, multiply by median library size and round to integers before
#'   running MaAsLin2. Default \code{TRUE}.
#' @param orders Optional named list of factor orders for metadata variables.
#'   Names should match columns in sample metadata; values are character vectors
#'   giving the desired levels (e.g., \code{list(Timepoint = c("Pre","Post"))}).
#' @param out_dir Optional output directory to write MaAsLin2 results. If \code{NULL},
#'   a dated project path is created (e.g., \code{<project_YYMMDD>/table/MaAsLin2/MaAsLin2.<name>}).
#' @param name Optional character suffix to distinguish runs (e.g., \code{"AdjBMI"}).
#'
#' @details
#' \strong{Relative→integer transform:} The helper detects relative abundance
#' heuristically by the mean of \code{sample_sums(psIN)} (≈100 implies relative).
#' If detected and \code{transform=TRUE}, counts are reconstructed by multiplying
#' the feature table by the median library size and rounding. Set \code{transform=FALSE}
#' to skip this step.
#'
#' \strong{Ordering factors:} If \code{orders} is supplied, matching metadata
#' columns are refactored to the provided level order (levels not present are ignored).
#'
#' Outputs are written under a dated project directory; MaAsLin2’s heatmap and
#' scatter plots are enabled by default.
#'
#' @return (Invisibly) the \code{Maaslin2} result object (list). Side effect:
#'   writes tables/plots to \code{out_dir}.
#'
#' @seealso \code{\link[Maaslin2]{Maaslin2}}
#'
#' @examples
#' \dontrun{
#' # Suppose ps is a phyloseq object with metadata columns Group and Subject
#' res <- Go_Maaslin2(
#'   psIN = ps,
#'   project = "IBD",
#'   fixed_effects = c("Group"),
#'   random_effects = c("Subject"),
#'   normalization = "TSS",
#'   transform = TRUE,
#'   orders = list(Group = c("Control","Case")),
#'   name = "BaseModel"
#' )
#' }
#'
#' @importFrom phyloseq sample_data taxa_are_rows otu_table sample_sums
#' @importFrom stats median
#' @importFrom utils file_test
#' @export

Go_Maaslin2 <- function(psIN,
                        project,
                        fixed_effects,
                        random_effects = NULL,
                        normalization = "TSS",   # "TSS","CLR","CSS","NONE"
                        transform = TRUE,        # 상대 abundance면 median depth 곱해 정수화
                        orders = NULL,           # 예: list(Timepoint=c("Pre","Post"))
                        out_dir = NULL,
                        name = NULL) {           # <- 새 옵션 추가

  ## ---- 0) 패키지 로드/설치 ----
  if (!requireNamespace("Maaslin2", quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
    BiocManager::install("Maaslin2")
  }
  for (pkg in c("methods","stats")) {
    if (!requireNamespace(pkg, quietly = TRUE)) stop(sprintf("Install %s first.", pkg))
  }

  suppressPackageStartupMessages({
    library(Maaslin2)
  })

  ## ---- 1) 출력 경로 ----
  out_root  <- file.path(sprintf("%s_%s", project, format(Sys.Date(), "%y%m%d")))
  if (!file_test("-d", out_root)) dir.create(out_root)
  out_table <- file.path(out_root, "table"); if (!file_test("-d", out_table)) dir.create(out_table)
  out_DA    <- file.path(out_table, "MaAsLin2"); if (!file_test("-d", out_DA)) dir.create(out_DA, recursive = TRUE)

  if (is.null(name)) {
    subdir <- "MaAsLin2.Base"
  } else {
    subdir <- sprintf("MaAsLin2.%s", name)
  }

  out_dir_final <- file.path(out_DA, subdir)
  if (!file_test("-d", out_dir_final)) dir.create(out_dir_final, recursive = TRUE)

  if (is.null(out_dir)) out_dir <- out_dir_final

  message(sprintf("[INFO] MaAsLin2 outputs -> %s", out_dir))

  ## ---- 2) 데이터 추출/정렬 ----
  metadata_df <- as.data.frame(sample_data(psIN))
  otu_mat <- if (taxa_are_rows(psIN)) as.matrix(otu_table(psIN)) else t(as.matrix(otu_table(psIN)))
  otu_mat <- as.data.frame(otu_mat)

  common_samples <- intersect(colnames(otu_mat), rownames(metadata_df))
  if (length(common_samples) < 2) stop("[ERROR] Not enough overlapping samples between OTU and metadata.")
  otu_mat     <- otu_mat[, common_samples, drop = FALSE]
  metadata_df <- metadata_df[ common_samples, , drop = FALSE]

  ## ---- 3) 상대 abundance → 정수 (옵션) ----
  detect_abundance_type <- function(physeq) {
    lib_sizes <- sample_sums(physeq)
    mean_lib <- mean(lib_sizes)
    if (abs(mean_lib - 100) < 0.1) "relative" else if (mean_lib > 1000) "absolute" else "unknown"
  }
  if (detect_abundance_type(psIN) == "relative" && isTRUE(transform)) {
    total_reads <- median(sample_sums(psIN))
    message(sprintf("[INFO] Relative abundance detected. Multiplying by median depth (%s) and rounding.", total_reads))
    otu_mat <- round(otu_mat * total_reads)
  } else {
    message("[INFO] Skip integer transform (absolute data or transform=FALSE).")
  }

  ## ---- 4) metadata 전처리 & orders 적용 ----
  rn <- rownames(metadata_df)
  metadata_df <- as.data.frame(lapply(metadata_df, function(x) {
    if (is.character(x)) {
      if (all(grepl("^[-+]?[0-9]*\\.?[0-9]+$", x[!is.na(x)]))) as.numeric(x) else as.factor(x)
    } else if (is.logical(x)) {
      as.factor(x)
    } else {
      x
    }
  }), stringsAsFactors = FALSE)
  rownames(metadata_df) <- rn

  if (!is.null(orders) && length(orders) > 0 && length(fixed_effects) > 0) {
    for (fx in fixed_effects) {
      if (fx %in% names(orders) && fx %in% colnames(metadata_df)) {
        levs <- orders[[fx]]
        if (is.factor(metadata_df[[fx]]) || is.character(metadata_df[[fx]])) {
          cur <- as.factor(metadata_df[[fx]])
          keep <- levs[levs %in% levels(cur)]
          if (length(keep) >= 2) metadata_df[[fx]] <- factor(cur, levels = keep)
        }
      }
    }
  }

  if (is.null(rownames(metadata_df)) || any(is.na(rownames(metadata_df)))) {
    stop("[ERROR] metadata must have valid rownames matching sample IDs.")
  }
  metadata_df <- metadata_df[colnames(otu_mat), , drop = FALSE]

  ## ---- 5) 효과 항목 확인 ----
  missing_fx <- setdiff(fixed_effects, colnames(metadata_df))
  if (length(missing_fx)) stop("[ERROR] Missing fixed_effects in metadata: ", paste(missing_fx, collapse=", "))

  if (!is.null(random_effects)) {
    missing_re <- setdiff(random_effects, colnames(metadata_df))
    if (length(missing_re)) {
      warning("[WARN] Some random_effects not found and will be dropped: ", paste(missing_re, collapse=", "))
      random_effects <- setdiff(random_effects, missing_re)
      if (length(random_effects) == 0) random_effects <- NULL
    }
  }

  ## ---- 6) MaAsLin2 실행 ----
  fit <- Maaslin2::Maaslin2(
    input_data     = otu_mat,
    input_metadata = metadata_df,
    output         = out_dir,
    fixed_effects  = fixed_effects,
    random_effects = random_effects,
    normalization  = normalization,
    transform      = "LOG",
    plot_heatmap   = TRUE,
    plot_scatter   = TRUE
  )

  message("[INFO] MaAsLin2 finished. Plots/results saved in: ", out_dir)
  invisible(fit)
}
