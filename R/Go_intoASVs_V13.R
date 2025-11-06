#' Go_intoASVs
#'
#' Compute intra-taxon nucleotide diversity (π) from ASV-level 16S sequences
#' within each taxon, optionally trimming noisy V3–V4 termini, clustering,
#' and bootstrapping CIs.
#'
#' @author Heekuk Park <hp2523@cumc.columbia.edu>
#' @date 2025-11-01
#'
#' @description
#' For each target taxon (or all taxa) in a \code{phyloseq} object, this function:
#' \itemize{
#'   \item extracts ASV sequences (\code{refseq(psIN)} or rownames),
#'   \item clusters intra-taxon sequences at \code{clustering_cutoff} (default 99.5\%),
#'   \item aligns sequences via \code{DECIPHER::AlignSeqs} or external \code{MAFFT},
#'   \item trims both ends by \code{trim_nt} nt to remove V3–V4 terminal noise,
#'   \item computes pairwise distances (\code{"simple"} = Hamming; or \code{"nucdiv"} = JC69, etc.),
#'   \item calculates per-sample nucleotide diversity (π) weighted by relative abundance
#'         or entropy, and (optionally) bootstrapped confidence intervals,
#'   \item saves π and ASV-count matrices, and merges them into \code{sample_data(psIN)}.
#' }
#'
#' @param psIN A \code{phyloseq} object with ASV-level data and (optionally)
#'   DNA sequences in \code{refseq(psIN)}.
#' @param project Character; project prefix used to create an output folder
#'   \verb{<project_YYMMDD>/table/pi_tab/}.
#' @param level Character; taxonomic rank to group ASVs before computing π
#'   (e.g., \code{"Genus"}, \code{"Family"}). Default \code{"Genus"}.
#' @param target Character; specific taxon name (e.g., \code{"Lactobacillus"}).
#'   Use \code{"all"} to compute across all taxa. Default \code{"all"}.
#' @param method Character; one of \code{c("simple","nucdiv")}.
#'   \code{"simple"} uses per-position Hamming distances; \code{"nucdiv"} uses
#'   substitution models from \code{ape::dist.dna()}.
#' @param aligner Character; alignment engine, \code{"DECIPHER"} (R-based) or
#'   \code{"MAFFT"} (external, faster if installed). Default \code{"DECIPHER"}.
#' @param min_asv Integer; minimum number of ASVs required to compute π. Default \code{3}.
#' @param min_abund Numeric; minimum total abundance per sample for inclusion. Default \code{10}.
#' @param clustering_cutoff Numeric; similarity cutoff for intra-taxon clustering
#'   via \code{DECIPHER::IdClusters}. Default \code{0.995}.
#' @param trim_nt Integer; number of nucleotides to trim from both sequence ends
#'   (useful for removing noisy termini of V3–V4 amplicons). Default \code{8}.
#' @param distance_gap Character; whether to \code{"exclude"} or \code{"include"} gaps
#'   when computing nucleotide distances. Default \code{"exclude"}.
#' @param distance_model Character; substitution model for \code{method="nucdiv"}
#'   (\code{"raw"}, \code{"JC69"}, etc.). Default \code{"raw"}.
#' @param weighting Character; weighting scheme for π:
#'   \code{"abundance"} (relative abundance) or \code{"entropy"} (p·(1–p)). Default \code{"abundance"}.
#' @param compute_ci Logical; if \code{TRUE}, compute 95\% bootstrap confidence intervals
#'   (\code{n_boot} replicates). Default \code{FALSE}.
#' @param n_boot Integer; number of bootstrap replicates for CI. Default \code{200}.
#' @param seed Integer; random seed. Default \code{123}.
#' @param n_cores Integer; parallel cores for per-taxon computation. Default \code{4}.
#'
#' @details
#' π (nucleotide diversity) is computed per sample as:
#' \deqn{π = 2 \sum_{i<j} w_i w_j d_{ij}}
#' where \(w_i\) are normalized weights (abundance- or entropy-based) and
#' \(d_{ij}\) is pairwise sequence distance.
#'
#' When \code{compute_ci=TRUE}, 95% CIs are estimated by bootstrap resampling
#' of ASVs within each sample. Clustering removes redundant ASVs within taxa
#' before alignment to reduce overcounting identical sequences.
#'
#' The function always saves results in:
#' \itemize{
#'   \item \code{pi_matrix_<method>_<level>_<target>_<date>.csv}
#'   \item \code{asv_count_matrix_<level>_<target>_<date>.csv}
#'   \item optional \code{pi_ci_low_matrix_*.csv}, \code{pi_ci_high_matrix_*.csv}
#'   \item text log: \code{pi_log_<method>_<level>_<target>_<date>.txt}
#' }
#' and merges the computed π and ASV counts into \code{sample_data(psIN)} as:
#' new columns prefixed by \code{pi_} and \code{asvN_}.
#'
#' @return A \code{phyloseq} object identical to input but with new columns added
#'   to \code{sample_data(psIN)} for each computed π and ASV count variable.
#'   Also writes CSVs and a log file under the project folder.
#'
#' @examples
#' \dontrun{
#' ps_pi <- Go_intoASVs(
#'   psIN = ps,
#'   project = "V34Study",
#'   level = "Genus",
#'   target = "Lactobacillus",
#'   method = "simple",
#'   aligner = "DECIPHER",
#'   trim_nt = 8,
#'   compute_ci = TRUE,
#'   n_boot = 100,
#'   n_cores = 4
#' )
#' sample_data(ps_pi)[, c("pi_Lactobacillus","asvN_Lactobacillus")]
#' }
#'
#' @importFrom phyloseq otu_table taxa_are_rows tax_table refseq sample_data
#' @importFrom DECIPHER AlignSeqs DistanceMatrix IdClusters
#' @importFrom Biostrings DNAStringSet writeXStringSet readDNAStringSet subseq letterFrequency
#' @importFrom pwalign stringDist
#' @importFrom ape dist.dna as.DNAbin
#' @importFrom parallel mclapply
#' @importFrom stats quantile
#' @importFrom utils write.csv
#' @importFrom methods is
#' @export


Go_intoASVs <- function(
    psIN,
    project,
    level = "Genus",
    target = "all",
    method = c("simple","nucdiv"),
    aligner = c("DECIPHER","MAFFT"),
    min_asv = 3,
    min_abund = 10,
    clustering_cutoff = 0.995,     # intra-taxon clustering
    trim_nt = 8,                   # End trimming to remove noise at the V3–V4 termini
    distance_gap = c("exclude","include"),
    distance_model = c("raw","JC69"),  # kept for compatibility (not used when TN93 is set below)
    weighting = c("abundance","entropy"),
    compute_ci = FALSE,
    n_boot = 200,
    seed = 123,
    n_cores = 4
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

  method         <- match.arg(method)
  aligner        <- match.arg(aligner)
  distance_gap   <- match.arg(distance_gap)
  distance_model <- match.arg(distance_model)
  weighting      <- match.arg(weighting)

  date_tag <- format(Sys.Date(), "%y%m%d")
  dir_base <- sprintf("%s_%s/table/pi_tab", project, date_tag)
  dir.create(dir_base, recursive = TRUE, showWarnings = FALSE)

  # ---------- IO helpers ----------
  sanitize <- function(x) gsub("[^A-Za-z0-9._-]+","_", x)
  safe_target  <- if (identical(target,"all")) "all" else sanitize(trimws(target))
  pi_file  <- sprintf("%s/pi_matrix_%s_%s_%s_%s.csv",  dir_base, method, level, safe_target, date_tag)
  asv_file <- sprintf("%s/asv_count_matrix_%s_%s_%s.csv", dir_base, level, safe_target, date_tag)
  log_file <- sprintf("%s/pi_log_%s_%s_%s_%s.txt", dir_base, method, level, safe_target, date_tag)

  cat(sprintf("[Go_intoASVs] %s | method=%s | level=%s | target=%s\n", Sys.time(), method, level, target),
      file = log_file, append = TRUE)

  # ---------- extract tables ----------
  otu_tab <- as(otu_table(psIN), "matrix")
  if (!taxa_are_rows(psIN)) otu_tab <- t(otu_tab)        # taxa x samples
  tax_tab <- as(tax_table(psIN), "matrix")

  # ---------- sequences ----------
  seqs_all <- tryCatch({
    as.character(refseq(psIN))
  }, error = function(e) {
    rownames(as(otu_table(psIN), "matrix"))
  })
  names(seqs_all) <- taxa_names(psIN)

  if (any(!grepl("^[ACGTN]+$", seqs_all))) {
    warning("The refseq slot is empty or taxa_names are not valid DNA sequences. Alignment may fail.")
    cat("WARN: refseq missing or taxa_names are not valid DNA sequences\n",
        file = log_file, append = TRUE)
  }

  # ---------- taxonomy ----------
  tax_raw <- tax_tab[, level, drop = TRUE]
  tax_raw[is.na(tax_raw)] <- "Unclassified"
  tax_labels <- trimws(tax_raw); names(tax_labels) <- rownames(tax_tab)
  taxa_pool <- unique(tax_labels)
  taxa_targets <- if (identical(target,"all")) taxa_pool else trimws(target)
  message("Target taxa: ", paste(taxa_targets, collapse = ", "))

  # ---------- utils ----------
  trim_alignment <- function(aln, k) {
    if (k <= 0) return(aln)
    w <- width(aln)[1]; if (k >= w/2) return(aln)
    rng <- IRanges(start = 1 + k, end = w - k)
    DNAStringSet(subseq(aln, start = start(rng), end = end(rng)))
  }
  entropy_weights <- function(p) {
    v <- p * (1 - p); if (sum(v) == 0) return(rep(0, length(p))); v / sum(v)
  }
  safe_tapply <- function(x, i, j, fun) {
    res <- try(tapply(x, list(i, j), fun), silent = TRUE)
    if (inherits(res, "try-error") || is.null(dimnames(res))) {
      res <- matrix(NA, nrow = length(unique(i)), ncol = length(unique(j)),
                    dimnames = list(unique(i), unique(j)))
    }
    res
  }

  # ---------- main loop over targets ----------
  results <- mclapply(taxa_targets, function(target_taxon){
    idx <- which(tax_labels == target_taxon)
    if (length(idx) < min_asv) {
      cat(sprintf("Skip %s: < min_asv\n", target_taxon), file = log_file, append = TRUE); return(NULL)
    }

    sub_abund <- otu_tab[idx, , drop = FALSE]      # taxa x samples
    sub_abund <- t(sub_abund)                      # samples x taxa(ASV)
    keep_samples <- rowSums(sub_abund) >= min_abund
    sub_abund <- sub_abund[keep_samples, , drop = FALSE]
    if (nrow(sub_abund) == 0) {
      cat(sprintf("Skip %s: total abundance < min_abund\n", target_taxon), file = log_file, append = TRUE);
      return(NULL)
    }

    # sequences subset
    seqs_sub <- seqs_all[colnames(sub_abund)]
    seqs_sub <- seqs_sub[!is.na(seqs_sub)]
    if (any(!grepl("^[ACGTN]+$", seqs_sub))) {
      cat(sprintf("Skip %s: invalid sequences\n", target_taxon), file = log_file, append = TRUE);
      return(NULL)
    }
    seqs <- DNAStringSet(seqs_sub)

    # ----- intra-taxon clustering (DECIPHER::IdClusters guard) -----
    if ("IdClusters" %in% getNamespaceExports("DECIPHER")) {
      dm <- DECIPHER::DistanceMatrix(seqs, includeTerminalGaps = FALSE)
      cl <- DECIPHER::IdClusters(dm, method = "complete",
                                 cutoff = 1 - clustering_cutoff, showPlot = FALSE)
      keep_cluster <- names(which.max(table(cl$cluster)))
      keep_ids <- rownames(cl)[cl$cluster == keep_cluster]
      seqs <- seqs[names(seqs) %in% keep_ids]
      sub_abund <- sub_abund[, colnames(sub_abund) %in% keep_ids, drop = FALSE]
      if (length(seqs) < min_asv) {
        cat(sprintf("Skip %s: < min_asv after clustering\n", target_taxon), file = log_file, append = TRUE);
        return(NULL)
      }
    } else {
      cat("NOTE: DECIPHER::IdClusters not exported; clustering skipped\n", file = log_file, append = TRUE)
    }

    # ----- alignment -----
    aln <- tryCatch({
      if (aligner == "MAFFT" && system("mafft --version", ignore.stdout = TRUE, ignore.stderr = TRUE) == 0) {
        tmp <- tempfile(fileext = ".fasta"); out <- tempfile(fileext = ".fasta")
        writeXStringSet(seqs, filepath = tmp)
        system(sprintf("mafft --auto --quiet %s > %s", tmp, out))
        readDNAStringSet(out)
      } else {
        AlignSeqs(seqs, iterations = 2, refinements = 2, verbose = FALSE)
      }
    }, error = function(e) {
      cat(sprintf("Align failed %s: %s\n", target_taxon, e$message), file = log_file, append = TRUE)
      NULL
    })
    if (is.null(aln)) return(NULL)

    # ----- trimming & gap QC -----
    if (trim_nt > 0) aln <- trim_alignment(aln, trim_nt)
    gap_prop <- tryCatch({
      sum(letterFrequency(aln, "-")) / (width(aln)[1] * length(aln))
    }, error = function(e) NA_real_)
    cat(sprintf("%s: Gap proportion %.2f%%\n", target_taxon, 100*gap_prop), file = log_file, append = TRUE)

    # ----- distance matrix -----
    if (method == "simple") {
      # Hamming-based proportional distance on aligned columns with optional gap exclusion
      A <- as.matrix(aln)
      n <- nrow(A)
      dist_mat <- matrix(0, n, n, dimnames = list(rownames(A), rownames(A)))
      for (i in seq_len(n - 1)) {
        Ai <- A[i, ]
        for (j in (i + 1):n) {
          Aj <- A[j, ]
          keep <- if (distance_gap == "exclude") (Ai != '-') & (Aj != '-') else rep(TRUE, length(Ai))
          Lij <- sum(keep)
          if (Lij > 0) {
            dist_mat[i, j] <- sum(Ai[keep] != Aj[keep]) / Lij
          } else {
            dist_mat[i, j] <- NA_real_
          }
          dist_mat[j, i] <- dist_mat[i, j]
        }
      }
      diag(dist_mat) <- 0
    } else {
      # Model-based distance with rate heterogeneity
      aln_chr <- as.matrix(aln)
      aln_dna <- as.DNAbin(aln_chr)
      dist_mat <- ape::dist.dna(
        aln_dna,
        model = "TN93",              # strengthened vs plain "raw"
        gamma = TRUE,                # site-rate heterogeneity
        pairwise.deletion = (distance_gap == "exclude"),
        as.matrix = TRUE
      )
    }

    # sync order
    ids <- intersect(colnames(sub_abund), colnames(dist_mat))
    if (length(ids) < 2) {
      cat(sprintf("Skip %s: <2 ids after sync\n", target_taxon), file = log_file, append = TRUE);
      return(NULL)
    }
    dist_mat  <- dist_mat[ids, ids, drop = FALSE]
    sub_abund <- sub_abund[, ids, drop = FALSE]
    L <- lower.tri(dist_mat, diag = FALSE)

    # ----- π calculator (returns a LIST, not numeric with attributes) -----
    calc_pi <- function(a){
      tot <- sum(a)
      if (tot <= 0) return(list(pi = NA_real_, ci_low = NA_real_, ci_high = NA_real_))

      p <- as.numeric(a) / tot
      w <- if (identical(weighting,"entropy")) entropy_weights(p) else p

      # abundance-weighted mean pairwise distance
      pi_val <- 2 * sum(outer(w, w)[L] * dist_mat[L], na.rm = TRUE)
      if (!is.finite(pi_val)) pi_val <- 0

      ci_low <- NA_real_; ci_high <- NA_real_

      if (compute_ci && sum(a > 0) >= 3) {
        boot_vals <- replicate(n_boot, {
          idx <- sample(seq_along(p), replace = TRUE)
          pb <- p[idx]
          db <- dist_mat[idx, idx, drop = FALSE]
          wb <- if (identical(weighting,"entropy")) entropy_weights(pb) else pb
          2 * sum(outer(wb, wb)[lower.tri(db)] * db[lower.tri(db)], na.rm = TRUE)
        })

        finite_any <- any(is.finite(boot_vals))
        if (finite_any) {
          qs <- stats::quantile(boot_vals, c(0.025, 0.975), na.rm = TRUE)
          ci_low <- as.numeric(qs[1]); ci_high <- as.numeric(qs[2])
        }
      }

      list(pi = as.numeric(pi_val), ci_low = ci_low, ci_high = ci_high)
    }

    # compute per-sample
    pi_list   <- lapply(rownames(sub_abund), function(s) calc_pi(sub_abund[s, ]))
    pi_values <- sapply(pi_list, function(x) x$pi)
    ci_low    <- sapply(pi_list, function(x) x$ci_low)
    ci_high   <- sapply(pi_list, function(x) x$ci_high)
    asv_counts <- rowSums(sub_abund > 0)

    df_out <- data.frame(
      Sample    = rownames(sub_abund),
      Taxon     = target_taxon,
      ASV_count = asv_counts,
      Pi        = as.numeric(pi_values),
      stringsAsFactors = FALSE
    )

    if (compute_ci) {
      df_out$Pi_CI_low  <- as.numeric(ci_low)
      df_out$Pi_CI_high <- as.numeric(ci_high)
    }

    df_out
  }, mc.cores = n_cores)

  results <- do.call(rbind, results)
  if (is.null(results) || nrow(results) == 0) {
    message("No valid taxa passed thresholds.")
    cat("No valid taxa passed thresholds\n", file = log_file, append = TRUE)
    return(psIN)
  }

  # target check
  if (!identical(target,"all")) {
    results <- results[trimws(results$Taxon) == trimws(target), , drop = FALSE]
    if (nrow(results) == 0) { message("No rows for requested target."); return(psIN) }
  }

  if (!all(c("Sample","Taxon","Pi","ASV_count") %in% colnames(results))) {
    warning("Results missing key columns — aborting reshape."); return(psIN)
  }

  # Deduplication and safe aggregation
  uniq_idx <- !duplicated(results[, c("Sample","Taxon"), drop = FALSE])
  results  <- results[uniq_idx, , drop = FALSE]

  pi_tab   <- as.data.frame(safe_tapply(results$Pi,        results$Sample, results$Taxon, identity))
  asv_tab  <- as.data.frame(safe_tapply(results$ASV_count, results$Sample, results$Taxon, identity))

  # Keep only the target column when a single target is specified
  if (!identical(target,"all")) {
    keep <- colnames(pi_tab) %in% trimws(target)
    pi_tab  <- pi_tab[,  keep, drop = FALSE]
    asv_tab <- asv_tab[, keep, drop = FALSE]
  }

  # Sorting and naming
  pi_tab  <- pi_tab[order(rownames(pi_tab)), , drop = FALSE]
  asv_tab <- asv_tab[order(rownames(asv_tab)), , drop = FALSE]
  colnames(pi_tab) <- paste0("pi_", colnames(pi_tab))

  # save matrices
  write.csv(pi_tab,  pi_file,  row.names = TRUE, na = "")
  write.csv(asv_tab, asv_file, row.names = TRUE, na = "")

  # Optionally save CI matrices — now taken from results DF (robust)
  if (compute_ci && all(c("Pi_CI_low","Pi_CI_high") %in% names(results))) {
    ci_low_tab  <- as.data.frame(safe_tapply(results$Pi_CI_low,  results$Sample, results$Taxon, identity))
    ci_high_tab <- as.data.frame(safe_tapply(results$Pi_CI_high, results$Sample, results$Taxon, identity))
    if (!identical(target,"all")) {
      keep <- colnames(ci_low_tab) %in% trimws(target)
      ci_low_tab  <- ci_low_tab[,  keep, drop = FALSE]
      ci_high_tab <- ci_high_tab[, keep, drop = FALSE]
    }
    colnames(ci_low_tab)  <- paste0("pi_CI_low_",  colnames(ci_low_tab))
    colnames(ci_high_tab) <- paste0("pi_CI_high_", colnames(ci_high_tab))
    write.csv(ci_low_tab,  sprintf("%s/pi_ci_low_matrix_%s_%s_%s.csv",  dir_base, method, level, date_tag),
              row.names = TRUE, na = "")
    write.csv(ci_high_tab, sprintf("%s/pi_ci_high_matrix_%s_%s_%s.csv", dir_base, method, level, date_tag),
              row.names = TRUE, na = "")
  }

  end_time <- Sys.time()

  summary_text <- sprintf("
--------------------------------------------------
π summary (%s/%s): mean=%.5f | sd=%.5f | range=%.5f–%.5f
Seed=%d | Runtime=%.2f min
π matrix saved: %s
ASV count matrix saved: %s
--------------------------------------------------
", method, weighting,
                          mean(results$Pi, na.rm = TRUE),
                          sd(results$Pi,  na.rm = TRUE),
                          min(results$Pi, na.rm = TRUE),
                          max(results$Pi, na.rm = TRUE),
                          seed,
                          round(as.numeric(difftime(end_time, start_time, units='mins')), 2),
                          pi_file, asv_file)

  cat(summary_text)
  cat(summary_text, file = log_file, append = TRUE)

  # ---- merge π & ASV_count into sample_data ----
  if (!is.null(sample_data(psIN, errorIfNULL = FALSE))) {
    sd <- as.data.frame(sample_data(psIN))

    pi_cols  <- pi_tab[rownames(sd), , drop = FALSE]
    asv_cols <- asv_tab[rownames(sd), , drop = FALSE]
    colnames(asv_cols) <- paste0("asvN_", colnames(asv_cols))

    # numeric coercion (NAs allowed)
    if (ncol(pi_cols) > 0)  pi_cols[]  <- lapply(pi_cols,  function(x) suppressWarnings(as.numeric(x)))
    if (ncol(asv_cols) > 0) asv_cols[] <- lapply(asv_cols, function(x) suppressWarnings(as.numeric(x)))

    # ensure unique names
    make_uniq <- function(df) { colnames(df) <- make.unique(colnames(df)); df }
    pi_cols  <- make_uniq(pi_cols)
    asv_cols <- make_uniq(asv_cols)

    sd_merged <- cbind(sd, pi_cols, asv_cols)
    sample_data(psIN) <- sd_merged

    cat(sprintf("\n[merge] Added to sample_data: %d pi cols, %d asv-count cols\n",
                ncol(pi_cols), ncol(asv_cols)))
  } else {
    warning("sample_data(psIN) is NULL — cannot merge π results into sample metadata.")
  }

  return(psIN)
}
