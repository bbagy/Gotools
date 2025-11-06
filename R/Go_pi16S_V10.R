Go_pi16S_v11.9b <- function(
    psIN,
    project,
    level        = "Genus",
    target       = "all",
    method       = c("simple","nucdiv"),
    aligner      = c("DECIPHER","MAFFT"),
    min_asv      = 3,
    min_abund    = 10,
    clustering_cutoff = 0.995,
    distance_gap = c("exclude","include"),
    weighting    = c("abundance","entropy"),
    seed         = 123,
    n_cores      = 4
){
  suppressPackageStartupMessages({
    library(phyloseq)
    library(DECIPHER)
    library(Biostrings)
    library(pwalign)
    library(ape)
    library(parallel)
  })
  
  start_time <- Sys.time()
  set.seed(seed)
  
  method       <- match.arg(method)
  aligner      <- match.arg(aligner)
  distance_gap <- match.arg(distance_gap)
  weighting    <- match.arg(weighting)
  
  date_tag <- format(Sys.Date(), "%y%m%d")
  dir_base <- sprintf("%s_%s/table/pi_tab", project, date_tag)
  dir.create(dir_base, recursive = TRUE, showWarnings = FALSE)
  
  otu_tab <- as(otu_table(psIN), "matrix")
  tax_tab <- as(tax_table(psIN), "matrix")
  
  ## refseq
  if (!is.null(refseq(psIN))) {
    seqs_all <- as.character(refseq(psIN))
    names(seqs_all) <- taxa_names(psIN)
  } else {
    seqs_all <- taxa_names(psIN)
    if (any(!grepl("^[ACGTN]+$", seqs_all))) {
      warning("⚠️ refseq is missing and ASV names are not DNA; alignment may fail.")
    }
  }
  
  ## taxonomy (공백 표준화)
  tax_raw <- tax_tab[, level]
  tax_raw[is.na(tax_raw)] <- "Unclassified"
  tax_labels <- trimws(tax_raw)
  names(tax_labels) <- taxa_names(psIN)
  
  target_std   <- if (identical(target, "all")) "all" else trimws(target)
  taxa_pool    <- unique(tax_labels)
  taxa_targets <- if (target_std == "all") taxa_pool else target_std
  message("Target taxa: ", paste(taxa_targets, collapse = ", "))
  
  ## 메인 루프 (✅ target 정확 적용)
  results <- mclapply(taxa_targets, function(target_taxon){
    idx <- which(tax_labels == target_taxon)
    if (length(idx) < min_asv) return(NULL)
    
    sub_abund <- otu_tab[idx, , drop = FALSE]   # taxa x samples
    sub_abund <- t(sub_abund)                   # samples x taxa (ASV)
    keep_samples <- rowSums(sub_abund) >= min_abund
    sub_abund <- sub_abund[keep_samples, , drop = FALSE]
    if (nrow(sub_abund) == 0) return(NULL)
    
    ## alignment 서열 준비
    seqs_sub <- seqs_all[colnames(sub_abund)]
    seqs_sub <- seqs_sub[!is.na(seqs_sub)]
    if (any(!grepl("^[ACGTN]+$", seqs_sub))) return(NULL)
    seqs <- DNAStringSet(seqs_sub)
    
    ## alignment
    aln <- tryCatch({
      if (aligner == "MAFFT" && system("mafft --version", ignore.stdout=TRUE, ignore.stderr=TRUE)==0) {
        tmp <- tempfile(fileext = ".fasta")
        writeXStringSet(seqs, filepath = tmp)
        system(paste("mafft --auto --quiet", tmp, ">", tmp))
        readDNAStringSet(tmp)
      } else {
        AlignSeqs(seqs, iterations = 2, refinements = 2, verbose = FALSE)
      }
    }, error = function(e) NULL)
    if (is.null(aln)) return(NULL)
    
    ## 거리/π 행렬
    if (method == "simple") {
      dist_mat <- pwalign::stringDist(aln, method = "hamming")
      dist_mat <- as.matrix(dist_mat)
    } else {
      aln_chr <- as.matrix(aln)
      aln_dna <- as.DNAbin(aln_chr)
      dist_mat <- ape::dist.dna(
        aln_dna,
        model = "raw",
        pairwise.deletion = (distance_gap == "exclude"),
        as.matrix = TRUE
      )
    }
    
    ## 순서 동기화
    ids <- intersect(colnames(sub_abund), colnames(dist_mat))
    if (length(ids) < 2) return(NULL)
    dist_mat  <- dist_mat[ids, ids, drop = FALSE]
    sub_abund <- sub_abund[, ids, drop = FALSE]
    
    L <- lower.tri(dist_mat, diag = FALSE)
    
    ## π 계산 (샘플별)
    calc_pi <- function(a){
      tot <- sum(a)
      if (tot <= 0) return(NA_real_)
      p <- as.numeric(a) / tot
      w <- if (identical(weighting, "entropy")) {
        v <- p * (1 - p)
        if (sum(v) == 0) return(0)
        v / sum(v)
      } else p
      val <- 2 * sum(outer(w, w)[L] * dist_mat[L], na.rm = TRUE)
      if (!is.finite(val)) val <- 0
      val
    }
    
    pi_values  <- sapply(rownames(sub_abund), function(s) calc_pi(sub_abund[s, ]))
    asv_counts <- rowSums(sub_abund > 0)   ## ✅ 샘플별 실제 ASV 수
    
    data.frame(Sample = rownames(sub_abund),
               Taxon = target_taxon,
               ASV_count = asv_counts,
               Pi = pi_values,
               stringsAsFactors = FALSE)
  }, mc.cores = n_cores)
  
  results <- do.call(rbind, results)
  if (is.null(results) || nrow(results) == 0) { 
    message("No valid taxa passed thresholds.")
    return(psIN)
  }
  
  ## target 재확인 (문자 이슈 방지)
  if (target_std != "all") {
    results <- results[trimws(results$Taxon) == target_std, , drop = FALSE]
    if (nrow(results) == 0) {
      message("No rows for requested target after filtering.")
      return(psIN)
    }
  }
  
  ## ---- 안전 집계 유틸 (reshape 방어) ----
  safe_tapply <- function(x, i, j, fun) {
    res <- try(tapply(x, list(i, j), fun), silent = TRUE)
    if (inherits(res, "try-error") || is.null(dimnames(res))) {
      res <- matrix(NA, nrow = length(unique(i)), ncol = length(unique(j)),
                    dimnames = list(unique(i), unique(j)))
    }
    res
  }
  
  ## 중복 제거 후 집계
  if (!all(c("Sample","Taxon") %in% colnames(results))) {
    warning("⚠️ Results missing key columns — skipping reshape.")
    return(psIN)
  }
  uniq_idx <- !duplicated(results[, c("Sample","Taxon"), drop = FALSE])
  results  <- results[uniq_idx, , drop = FALSE]
  
  pi_tab  <- as.data.frame(safe_tapply(results$Pi,        results$Sample, results$Taxon, identity))
  asv_tab <- as.data.frame(safe_tapply(results$ASV_count, results$Sample, results$Taxon, identity))
  
  ## 정렬 및 target 열만 유지
  pi_tab  <- pi_tab[order(rownames(pi_tab)), , drop = FALSE]
  asv_tab <- asv_tab[order(rownames(asv_tab)), , drop = FALSE]
  
  if (target_std != "all") {
    keep_col <- colnames(pi_tab) %in% target_std
    pi_tab   <- pi_tab[,  keep_col, drop = FALSE]
    asv_tab  <- asv_tab[, keep_col, drop = FALSE]
  }
  
  ## 열 이름 정리
  colnames(pi_tab) <- paste0("pi_", colnames(pi_tab))
  
  ## 파일명
  pi_file  <- sprintf("%s/pi_matrix_%s_%s_%s_%s.csv",  dir_base, method, level, target_std, date_tag)
  asv_file <- sprintf("%s/asv_count_matrix_%s_%s_%s.csv", dir_base, level, target_std, date_tag)
  
  ## 저장
  write.csv(pi_tab,  pi_file,  row.names = TRUE, na = "")
  write.csv(asv_tab, asv_file, row.names = TRUE, na = "")
  
  ## 요약
  end_time <- Sys.time()
  cat("\n--------------------------------------------------\n",
      sprintf("✅ π summary (%s/%s): mean=%.5f | sd=%.5f | range=%.5f–%.5f",
              method, weighting,
              suppressWarnings(mean(results$Pi, na.rm = TRUE)),
              suppressWarnings(sd(results$Pi,  na.rm = TRUE)),
              suppressWarnings(min(results$Pi, na.rm = TRUE)),
              suppressWarnings(max(results$Pi, na.rm = TRUE))), "\n",
      sprintf("Seed=%d | Runtime=%.2f min",
              seed, round(as.numeric(difftime(end_time, start_time, units='mins')), 2)), "\n",
      sprintf("π matrix saved: %s", pi_file), "\n",
      sprintf("ASV count matrix saved: %s", asv_file), "\n--------------------------------------------------\n")
  
  return(psIN)
}