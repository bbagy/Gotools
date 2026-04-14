# Legacy helper: archived from active Gotools exports on 2026-04-08.
# This file is intentionally kept outside package R/ so it is not exported.
#
# Purpose:
# - cluster ASV sequences into OTU-like bins at a requested identity threshold
# - rebuild a phyloseq object from the merged count table
# - optionally reassign taxonomy with dada2::assignTaxonomy()
#
# Notes:
# - This is retained only for legacy/bridging workflows that still require
#   97%/99%/100% style clustering.
# - Modern analyses should prefer native ASV workflows.

Go_cluster_legacy <- function(psIN, project, db, percent, nproc = 4L, out_dir = ".") {
  pkgs <- c("Biostrings", "DECIPHER", "dada2", "dplyr", "phyloseq", "speedyseq", "tibble")
  missing_pkgs <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_pkgs) > 0) {
    stop("Go_cluster_legacy() requires packages: ", paste(missing_pkgs, collapse = ", "))
  }

  if (!inherits(psIN, "phyloseq")) {
    stop("`psIN` must be a phyloseq object.")
  }
  if (!is.numeric(percent) || length(percent) != 1 || is.na(percent) || percent <= 0 || percent > 100) {
    stop("`percent` must be a single number in (0, 100].")
  }
  if (!is.character(project) || length(project) != 1 || !nzchar(project)) {
    stop("`project` must be a non-empty string.")
  }
  if (!is.character(db) || length(db) != 1 || !nzchar(db)) {
    stop("`db` must be a non-empty taxonomy database path.")
  }

  cutoff <- 1 - percent / 100
  nproc <- as.integer(nproc)
  if (!is.finite(nproc) || nproc < 1) nproc <- 1L

  rds_dir <- file.path(out_dir, "2_rds")
  csv_dir <- file.path(out_dir, "1_out")
  if (!dir.exists(rds_dir)) dir.create(rds_dir, recursive = TRUE)
  if (!dir.exists(csv_dir)) dir.create(csv_dir, recursive = TRUE)

  otu_mat <- as(phyloseq::otu_table(psIN), "matrix")
  taxa_are_rows <- phyloseq::taxa_are_rows(psIN)
  if (taxa_are_rows) {
    otu_mat <- t(otu_mat)
  }
  otu_mat <- as.matrix(otu_mat)

  asv_sequences <- colnames(otu_mat)
  if (is.null(asv_sequences) || length(asv_sequences) == 0) {
    stop("No ASV sequences found in `otu_table(psIN)` column names.")
  }

  dna <- Biostrings::DNAStringSet(asv_sequences)
  aln <- DECIPHER::AlignSeqs(dna, processors = nproc)
  dist_mat <- DECIPHER::DistanceMatrix(aln, processors = nproc)
  clusters <- DECIPHER::IdClusters(
    dist_mat,
    method = "complete",
    cutoff = cutoff,
    processors = nproc
  )

  cluster_tbl <- tibble::add_column(clusters, sequence = asv_sequences)
  merged_seqtab <- t(rowsum(t(otu_mat), clusters$cluster))

  clustered <- dplyr::distinct(cluster_tbl, cluster, .keep_all = TRUE)
  merged_seqtab_t <- as.data.frame(t(merged_seqtab), stringsAsFactors = FALSE)
  merged_seqtab_t$seqs <- clustered$sequence[
    match(as.character(rownames(merged_seqtab_t)), as.character(clustered$cluster))
  ]
  rownames(merged_seqtab_t) <- merged_seqtab_t$seqs
  merged_seqtab_t$seqs <- NULL
  seqtab_clustered <- as.matrix(t(merged_seqtab_t))

  seqs <- speedyseq::getSequences(seqtab_clustered)
  fasta_headers <- paste0(">", seqs)
  fasta <- c(rbind(fasta_headers, seqs))

  stamp <- format(Sys.Date(), "%y%m%d")
  project_tag <- sprintf("%s_%s", project, percent)
  fasta_file <- file.path(csv_dir, sprintf("ASVs%s.%s.%s.seqs.fna", percent, project, stamp))
  write(fasta, file = fasta_file)

  tax <- dada2::assignTaxonomy(
    seqtab_clustered,
    db,
    taxLevels = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
    minBoot = 80,
    verbose = TRUE,
    multithread = TRUE
  )

  ps_clustered <- phyloseq::phyloseq(
    phyloseq::otu_table(seqtab_clustered, taxa_are_rows = FALSE),
    phyloseq::tax_table(tax)
  )

  tax_df <- as.data.frame(phyloseq::tax_table(ps_clustered), stringsAsFactors = FALSE)
  if (all(c("Genus", "Species") %in% names(tax_df))) {
    tax_df$Species <- paste(tax_df$Genus, tax_df$Species)
    phyloseq::tax_table(ps_clustered) <- as.matrix(tax_df)
  }

  message("#-- Before clustering --#")
  print(psIN)
  message(sprintf("#-- After clustering at %s%% --#", percent))
  print(ps_clustered)

  rds_file <- file.path(rds_dir, sprintf("ps%s_filtered.%s.%s.rds", percent, project_tag, stamp))
  saveRDS(ps_clustered, rds_file)

  otu_df <- as.data.frame(t(phyloseq::otu_table(ps_clustered)), stringsAsFactors = FALSE)
  tax_df <- as.data.frame(phyloseq::tax_table(ps_clustered), stringsAsFactors = FALSE)
  otu_tax_table <- cbind(otu_df, tax_df)
  csv_file <- file.path(csv_dir, sprintf("ASVs%s_clustered.%s.%s.csv", percent, project_tag, stamp))
  utils::write.csv(otu_tax_table, quote = FALSE, col.names = NA, row.names = TRUE, file = csv_file)

  invisible(list(
    phyloseq = ps_clustered,
    files = list(
      fasta = fasta_file,
      rds = rds_file,
      csv = csv_file
    ),
    cluster_table = cluster_tbl
  ))
}
