Go_pi16S_v6 <- function(psIN,
                        project = "GoPi16S",
                        level = "Species",
                        target = "all",
                        method = c("simple", "nucdiv"),
                        aligner = c("auto", "DECIPHER", "Biostrings", "MAFFT"),
                        iterations = 5,
                        refinements = 3,
                        distance_gap = c("exclude", "include"),
                        distance_model = c("raw", "JC69"),
                        delimiter = ";",
                        min_asv = 2,
                        min_abund = 10,
                        # New: accuracy boosters
                        trim_nt = 8,                 # 좌우 트리밍 nt 수 (0이면 사용 안 함)
                        clustering_cutoff = 0.995,   # intra-taxon 유사도 커트오프(99.5%)
                        refine_pairwise = TRUE,      # 유사쌍 고정밀 pairwise 재정렬
                        refine_threshold = 0.99,     # 재정렬 대상 쌍의 identity 기준
                        weighting = c("abundance", "entropy"),
                        # CI options
                        compute_ci = FALSE,
                        n_boot = 200,
                        n_cores = 4,
                        seed = 123) {
  suppressPackageStartupMessages({
    library(phyloseq)
    library(DECIPHER)
    library(Biostrings)
    library(ape)
    library(parallel)
  })
  start_time <- Sys.time()
  set.seed(seed)
  method <- match.arg(method)
  aligner <- match.arg(aligner)
  distance_gap <- match.arg(distance_gap)
  distance_model <- match.arg(distance_model)
  weighting <- match.arg(weighting)
  
  # --- Taxonomy 유연 파싱 ---
  tax <- as.data.frame(tax_table(psIN))
  if (ncol(tax) == 1) {
    tax_str <- strsplit(as.character(tax[,1]), delimiter)
    tax <- do.call(rbind, lapply(tax_str, function(x) {
      x <- c(x, rep(NA, 7 - length(x)))
      setNames(x, c("Kingdom","Phylum","Class","Order","Family","Genus","Species"))
    }))
  }
  seqs <- refseq(psIN)
  asv_abund <- as(otu_table(psIN), "matrix")
  if (taxa_are_rows(psIN)) asv_abund <- t(asv_abund)
  
  # 대상 taxon 목록
  taxa_targets <- if (target == "all") unique(na.omit(tax[[level]])) else target
  
  sample_df <- as.data.frame(sample_data(psIN))
  out_pi_dir <- file.path(sprintf("%s_%s/table/pi_tab", project, format(Sys.Date(), "%y%m%d")))
  if (!file_test("-d", out_pi_dir)) dir.create(out_pi_dir, recursive = TRUE)
  
  # 유틸: 좌우 트리밍
  trim_alignment <- function(aln, k) {
    if (k <= 0) return(aln)
    # 양끝 k nt 제거 (gap 포함 위치 기준으로 단순 crop)
    w <- width(aln)[1]
    rng <- IRanges(start = 1 + k, end = w - k)
    DNAStringSet(subseq(aln, start = start(rng), end = end(rng)))
  }
  
  # 유틸: identity 계산
  pairwise_identity <- function(s1, s2) {
    aln <- pairwiseAlignment(s1, s2, substitutionMatrix = nucleotideSubstitutionMatrix(),
                             gapOpening = 5, gapExtension = 2)
    pid(aln) / 100
  }
  
  # 유틸: entropy weight
  entropy_weights <- function(p) {
    v <- p * (1 - p)
    if (sum(v) == 0) return(rep(0, length(p)))
    v / sum(v)
  }
  
  # 병렬 루프
  results <- mclapply(taxa_targets, function(target_taxon) {
    idx <- tax[[level]] == target_taxon
    seqs_sub <- seqs[idx]
    abund_sub <- asv_abund[, idx, drop = FALSE]
    if (length(seqs_sub) < min_asv) return(NULL)
    if (sum(abund_sub) < min_abund) return(NULL)
    
    # --- Intra-taxon clustering (오분류/먼 ASV 제거) ---
    # 거리행렬(빠른 raw/Hamming 기반)로 99.5% 유사도 클러스터 유지
    dm <- DECIPHER::DistanceMatrix(DNAStringSet(seqs_sub), includeTerminalGaps = FALSE)
    cl <- DECIPHER::IdClusters(dm, method = "complete", cutoff = 1 - clustering_cutoff, showPlot = FALSE)
    # 최빈 cluster만 채택 (가장 큰 동질군)
    keep_cluster <- names(which.max(table(cl$cluster)))
    keep_ids <- rownames(cl)[cl$cluster == keep_cluster]
    seqs_sub <- seqs_sub[names(seqs_sub) %in% keep_ids]
    abund_sub <- abund_sub[, colnames(abund_sub) %in% keep_ids, drop = FALSE]
    if (length(seqs_sub) < min_asv) return(NULL)
    
    # --- Alignment ---
    aln <- NULL
    aln_path <- tempfile(fileext = ".fasta")
    writeXStringSet(seqs_sub, aln_path)
    
    if (aligner == "MAFFT") {
      if (system("mafft --version", ignore.stdout = TRUE, ignore.stderr = TRUE) != 0) {
        warning("⚠️ MAFFT not found — fallback to DECIPHER::AlignSeqs()")
        aln <- tryCatch(
          AlignSeqs(seqs_sub, iterations = iterations, refinements = refinements, verbose = FALSE),
          error = function(e) NULL
        )
      } else {
        aln_out <- tempfile(fileext = ".fasta")
        cmd <- sprintf("mafft --auto --reorder --quiet %s > %s", aln_path, aln_out)
        system(cmd)
        aln <- readDNAStringSet(aln_out)
      }
    } else if (aligner == "auto") {
      if (length(seqs_sub) <= 5) {
        aln <- AlignPairwiseDNAStrings(DNAStringSet(seqs_sub))
      } else {
        aln <- tryCatch(
          AlignSeqs(seqs_sub, iterations = iterations, refinements = refinements, verbose = FALSE),
          error = function(e) NULL
        )
      }
    } else if (aligner == "DECIPHER") {
      aln <- tryCatch(
        AlignSeqs(seqs_sub, iterations = iterations, refinements = refinements, verbose = FALSE),
        error = function(e) NULL
      )
    } else {
      aln <- AlignPairwiseDNAStrings(DNAStringSet(seqs_sub))
    }
    if (is.null(aln)) return(NULL)
    
    # --- 하이퍼변이 end bias 제거용 trimming ---
    if (trim_nt > 0) aln <- trim_alignment(aln, trim_nt)
    
    # --- (옵션) 유사쌍 pairwise 재정렬로 미세보정 ---
    if (refine_pairwise && length(aln) >= 2) {
      # 고유 쌍만 선별하여 pid > refine_threshold 재정렬 (성능 고려: 서브샘플링 가능)
      ids <- names(aln)
      to_refine <- list()
      if (length(ids) > 1) {
        for (i in 1:(length(ids)-1)) {
          for (j in (i+1):length(ids)) {
            pid_ <- pairwise_identity(seqs[ids[i]], seqs[ids[j]])
            if (is.finite(pid_) && pid_ >= refine_threshold) {
              # 재정렬 후 더 정확한 local alignment로 치환 (여기선 전체 aln 수정 대신 거리 계산 단계에서만 사용)
              # 단, 전체 MSA 반영은 비용↑이므로 아래 거리 계산에서 on-the-fly로 반영
              to_refine[[length(to_refine)+1]] <- c(ids[i], ids[j])
            }
          }
        }
      }
    }
    
    # --- gap 비율 QC ---
    gap_prop <- sum(letterFrequency(aln, "-")) / (width(aln)[1] * length(aln))
    
    # --- π 계산 ---
    if (method == "simple") {
      # 거리행렬 계산: gap 처리/모델 선택
      if (distance_gap == "exclude") {
        aln_bin <- as.DNAbin(as.character(aln))
        dist_mat <- as.matrix(dist.dna(aln_bin,
                                       model = ifelse(distance_model=="raw","N","JC69"),
                                       pairwise.deletion = TRUE))
      } else {
        dist_mat <- as.matrix(stringDist(aln, method = "hamming")) / width(aln)[1]
      }
      
      # (선택) 재정렬 대상 쌍의 distance 미세보정: pairwiseAlignment 기반 재계산
      if (refine_pairwise && length(aln) >= 2 && length(to_refine) > 0) {
        for (pair in to_refine) {
          i <- which(names(aln) == pair[1]); j <- which(names(aln) == pair[2])
          pa <- pairwiseAlignment(seqs[pair[1]], seqs[pair[2]],
                                  substitutionMatrix = nucleotideSubstitutionMatrix(),
                                  gapOpening = 5, gapExtension = 2)
          mm <- nmis(pa); L <- nchar(pattern(pa)) # 간단 근사
          d_ij <- mm / L
          dist_mat[i, j] <- dist_mat[j, i] <- d_ij
        }
      }
      
      # 샘플별 가중치 선택 (abundance vs entropy)
      pi_values <- sapply(rownames(abund_sub), function(sam) {
        abund <- abund_sub[sam, ]
        if (sum(abund) == 0) return(NA_real_)
        p <- abund / sum(abund)
        w <- if (weighting == "abundance") {
          # 쌍별 가중 식: Σ(2 p_i p_j d_ij)
          NULL
        } else {
          # entropy 가중 (ASV 단위 가중치 → 쌍으로 확장: 2 w_i w_j d_ij)
          entropy_weights(p)
        }
        val <- 0
        k <- length(p)
        for (i in 1:(k-1)) {
          for (j in (i+1):k) {
            if (weighting == "abundance") {
              val <- val + 2 * p[i] * p[j] * dist_mat[i, j]
            } else {
              val <- val + 2 * w[i] * w[j] * dist_mat[i, j]
            }
          }
        }
        # (옵션) CI 부트스트랩
        if (compute_ci && k >= 3) {
          boot <- replicate(n_boot, {
            idx <- sample(1:k, replace = TRUE)
            pb <- p[idx]; db <- dist_mat[idx, idx, drop = FALSE]
            if (weighting == "abundance") {
              vb <- 0
              for (ii in 1:(k-1)) for (jj in (ii+1):k) vb <- vb + 2 * pb[ii] * pb[jj] * db[ii, jj]
              vb
            } else {
              wb <- entropy_weights(pb)
              vb <- 0
              for (ii in 1:(k-1)) for (jj in (ii+1):k) vb <- vb + 2 * wb[ii] * wb[jj] * db[ii, jj]
              vb
            }
          })
          attr(val, "ci") <- quantile(boot, c(0.025, 0.975), na.rm = TRUE)
        }
        val
      })
      
      out <- data.frame(Sample = names(pi_values),
                        Taxon = target_taxon,
                        Method = paste0(method, "_", weighting),
                        Pi = as.numeric(pi_values),
                        GapProp = gap_prop,
                        ColName = paste0("pi_", method, "_", weighting, "_", gsub(" ", "_", target_taxon)),
                        stringsAsFactors = FALSE)
      # CI 붙이기(있다면)
      if (compute_ci) {
        ciL <- sapply(pi_values, function(x) if (!is.null(attr(x,"ci"))) attr(x,"ci")[1] else NA_real_)
        ciU <- sapply(pi_values, function(x) if (!is.null(attr(x,"ci"))) attr(x,"ci")[2] else NA_real_)
        out$Pi_CI_low <- as.numeric(ciL); out$Pi_CI_high <- as.numeric(ciU)
      }
      
    } else {
      # nucdiv: per-site heterozygosity
      aln_bin <- as.DNAbin(as.character(aln))
      true_pi <- tryCatch(nuc.div(aln_bin), error = function(e) NA_real_)
      
      pi_unw <- rep(true_pi, nrow(abund_sub)); names(pi_unw) <- rownames(abund_sub)
      pi_w_ab <- sapply(rownames(abund_sub), function(sam) {
        a <- abund_sub[sam, ]; if (sum(a) == 0 || is.na(true_pi)) return(NA_real_)
        (sum(a) / sum(abund_sub)) * true_pi
      })
      pi_w_en <- sapply(rownames(abund_sub), function(sam) {
        a <- abund_sub[sam, ]; if (sum(a) == 0 || is.na(true_pi)) return(NA_real_)
        p <- a / sum(a); w <- entropy_weights(p); # 스칼라로 축약: 총 엔트로피 가중 합의 평균계수
        # 엔트로피 가중의 강도 반영 (정규화된 평균 가중치)
        mean_w <- mean(w[w>0])
        true_pi * mean_w * length(w[w>0])
      })
      
      out <- data.frame(Sample = rownames(abund_sub),
                        Taxon = target_taxon,
                        Method = method,
                        Pi_unweighted = as.numeric(pi_unw[rownames(abund_sub)]),
                        Pi_weighted_abundance = as.numeric(pi_w_ab[rownames(abund_sub)]),
                        Pi_weighted_entropy = as.numeric(pi_w_en[rownames(abund_sub)]),
                        GapProp = gap_prop,
                        Col_unweighted = paste0("pi_", method, "_", gsub(" ", "_", target_taxon)),
                        Col_w_ab = paste0("pi_", method, "_wAb_", gsub(" ", "_", target_taxon)),
                        Col_w_en = paste0("pi_", method, "_wEn_", gsub(" ", "_", target_taxon)),
                        stringsAsFactors = FALSE)
      # (선택) CI는 nucdiv 전체값에 적용 의미가 낮아 기본 비활성
    }
    
    out
  }, mc.cores = n_cores)
  
  results <- do.call(rbind, results)
  if (is.null(results)) { message("No valid taxa passed thresholds."); return(psIN) }
  
  # sample_data에 쓰기
  if (method == "simple") {
    for (col in unique(results$ColName)) {
      sub <- results[results$ColName == col, ]
      vals <- setNames(sub$Pi, sub$Sample)
      sample_df[[col]] <- vals[rownames(sample_df)]
      if (compute_ci) {
        sample_df[[paste0(col,"_CI_low")]]  <- setNames(sub$Pi_CI_low,  sub$Sample)[rownames(sample_df)]
        sample_df[[paste0(col,"_CI_high")]] <- setNames(sub$Pi_CI_high, sub$Sample)[rownames(sample_df)]
      }
    }
  } else {
    for (col in unique(results$Col_unweighted)) {
      sub <- results[results$Col_unweighted == col, ]
      sample_df[[col]] <- setNames(sub$Pi_unweighted, sub$Sample)[rownames(sample_df)]
    }
    for (col in unique(results$Col_w_ab)) {
      sub <- results[results$Col_w_ab == col, ]
      sample_df[[col]] <- setNames(sub$Pi_weighted_abundance, sub$Sample)[rownames(sample_df)]
    }
    for (col in unique(results$Col_w_en)) {
      sub <- results[results$Col_w_en == col, ]
      sample_df[[col]] <- setNames(sub$Pi_weighted_entropy, sub$Sample)[rownames(sample_df)]
    }
  }
  sample_data(psIN) <- sample_data(sample_df)
  
  # 저장
  pi_mat <- sample_df[, grep("^pi_", colnames(sample_df)), drop = FALSE]
  out_file <- file.path(out_pi_dir, paste0("pi_matrix_", method, "_", level, ".csv"))
  write.csv(pi_mat, out_file, row.names = TRUE)
  
  # QC 로그
  pi_vals_all <- as.numeric(as.matrix(pi_mat))
  pi_mean <- mean(pi_vals_all, na.rm = TRUE)
  pi_sd   <- sd(pi_vals_all, na.rm = TRUE)
  pi_range <- range(pi_vals_all, na.rm = TRUE)
  end_time <- Sys.time()
  runtime <- round(as.numeric(difftime(end_time, start_time, units = "mins")), 2)
  
  message("--------------------------------------------------")
  message(sprintf("✅ π summary (%s/%s): mean=%.4f | sd=%.4f | range=%.4f–%.4f",
                  method, ifelse(method=="simple", weighting, "multi"), pi_mean, pi_sd, pi_range[1], pi_range[2]))
  if ("GapProp" %in% names(results)) {
    message(sprintf("✅ Avg alignment gap proportion = %.2f%%", mean(results$GapProp, na.rm = TRUE) * 100))
  }
  message(sprintf("Seed=%s | n_cores=%d | Runtime=%.2f min", seed, n_cores, runtime))
  message(sprintf("Results saved to: %s", out_file))
  message("--------------------------------------------------")
  
  return(psIN)
}