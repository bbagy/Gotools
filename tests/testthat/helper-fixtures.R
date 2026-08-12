## Many Gotools internals call phyloseq/ape/etc. functions unqualified,
## relying on Gotool_dependency() having attached them via library() first
## (as every Download_templates script does). Mirror that here without the
## network-dependent install-checking Gotool_dependency() would otherwise do.
## helper-*.R files are auto-sourced by testthat before tests run, under
## both devtools::test() and R CMD check (unlike tests/testthat.R, which
## devtools::test() does not execute).
for (.gt_pkg in c("ape", "dplyr", "ggplot2", "grid", "patchwork", "phyloseq", "tidyr", "vegan", "circlize", "broom", "glue")) {
  suppressPackageStartupMessages(library(.gt_pkg, character.only = TRUE))
}
rm(.gt_pkg)

## Synthetic input fixtures for the end-to-end report smoke tests.
## Kept small (n_taxa/n_sample) to keep test runtime reasonable.

make_16S_fixture <- function(root) {
  dir.create(file.path(root, "1_out", "TREE"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(root, "3_map"), recursive = TRUE, showWarnings = FALSE)
  old_wd <- getwd()
  on.exit(setwd(old_wd), add = TRUE)
  setwd(root)

  set.seed(1)
  n_asv <- 20
  n_sample <- 10
  today <- format(Sys.Date(), "%y%m%d")

  rand_dna <- function(n, len) {
    vapply(seq_len(n), function(i) paste(sample(c("A", "T", "C", "G"), len, replace = TRUE), collapse = ""), character(1))
  }
  asv_seqs <- unique(rand_dna(n_asv, 420))
  while (length(asv_seqs) < n_asv) asv_seqs <- unique(c(asv_seqs, rand_dna(n_asv - length(asv_seqs), 420)))
  asv_seqs <- asv_seqs[1:n_asv]

  tax <- data.frame(
    Kingdom = "Bacteria",
    Phylum  = sample(c("Firmicutes", "Bacteroidota"), n_asv, replace = TRUE),
    Class   = paste0("Class", sample(1:3, n_asv, replace = TRUE)),
    Order   = paste0("Order", sample(1:3, n_asv, replace = TRUE)),
    Family  = paste0("Family", sample(1:3, n_asv, replace = TRUE)),
    Genus   = sample(c("Bacteroides", "Faecalibacterium", "Blautia"), n_asv, replace = TRUE),
    Species = paste0("sp", 1:n_asv),
    row.names = asv_seqs, check.names = FALSE
  )
  sample_ids <- sprintf("Sample%02d", 1:n_sample)
  counts <- matrix(stats::rpois(n_asv * n_sample, lambda = 300), n_asv, n_sample, dimnames = list(asv_seqs, sample_ids))
  counts[1:3, ] <- counts[1:3, ] + 2000L
  asv_table <- cbind(tax, as.data.frame(counts, check.names = FALSE))

  write.csv(asv_table, sprintf("1_out/TestFixture16S.final_asvTable.%s.csv", today), row.names = TRUE)
  write.csv(data.frame(Sample = sample_ids, RawReads = 30000, FilteredReads = 25000), "1_out/track.csv", row.names = FALSE)
  ape::write.tree(ape::rtree(n_asv, tip.label = asv_seqs), "1_out/TREE/tree.nwk")
  write.csv(
    data.frame(row.names = sample_ids, TreatmentGroup = rep(c("GroupA", "GroupB"), each = n_sample / 2), check.names = FALSE),
    "3_map/MAPPING.csv", row.names = TRUE
  )
  invisible(root)
}

make_MG_fixture <- function(root) {
  dir.create(file.path(root, "1_out"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(root, "3_map"), recursive = TRUE, showWarnings = FALSE)
  old_wd <- getwd()
  on.exit(setwd(old_wd), add = TRUE)
  setwd(root)

  set.seed(2)
  n_sp <- 15
  n_sample <- 10
  sample_ids <- sprintf("Sample%02d", 1:n_sample)
  species_names <- paste0(sample(c("Bacteroides", "Faecalibacterium", "Prevotella"), n_sp, replace = TRUE), "_species", 1:n_sp)
  lineage <- sprintf(
    "d__Bacteria|p__%s|c__Class%d|o__Order%d|f__Family%d|g__%s|s__%s",
    sample(c("Firmicutes", "Bacteroidota"), n_sp, replace = TRUE),
    sample(1:3, n_sp, replace = TRUE), sample(1:3, n_sp, replace = TRUE), sample(1:3, n_sp, replace = TRUE),
    sub("_species.*$", "", species_names), species_names
  )
  write.table(data.frame(ID = lineage, check.names = FALSE), "1_out/KRAKEN2_MPA.txt", sep = "\t", quote = FALSE, row.names = FALSE)

  counts <- matrix(stats::rpois(n_sp * n_sample, lambda = 300), n_sp, n_sample, dimnames = list(NULL, sample_ids))
  counts[1:3, ] <- counts[1:3, ] + 2000L
  write.table(data.frame(ID = paste0("s__", species_names), counts, check.names = FALSE),
              "1_out/BRACKEN_MPA.txt", sep = "\t", quote = FALSE, row.names = FALSE)

  write.csv(data.frame(Sample = sample_ids, RawReads = 30000, FilteredReads = 25000), "1_out/QC_summary.csv", row.names = FALSE)
  write.csv(
    data.frame(row.names = sample_ids, TreatmentGroup = rep(c("GroupA", "GroupB"), each = n_sample / 2), Sex = rep(c("F", "M"), n_sample / 2), check.names = FALSE),
    "3_map/MAPPING.csv", row.names = TRUE
  )

  n_path <- 6
  path_ids <- sprintf("PWY-%04d: Pathway description %d", 1:n_path, 1:n_path)
  path_counts <- matrix(stats::rpois(n_path * n_sample, lambda = 50), n_path, n_sample, dimnames = list(path_ids, sample_ids))
  write.table(path_counts, "1_out/PATHABUNDANCE.tsv", sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)

  n_ko <- 6
  ko_ids <- sprintf("K%05d: KO description %d", 1:n_ko, 1:n_ko)
  ko_counts <- matrix(stats::rpois(n_ko * n_sample, lambda = 50), n_ko, n_sample, dimnames = list(ko_ids, sample_ids))
  write.table(ko_counts, "1_out/GENEFAMILIES.tsv", sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)

  invisible(root)
}

make_RNAseq_fixture <- function(root) {
  dir.create(file.path(root, "1_out"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(root, "3_map"), recursive = TRUE, showWarnings = FALSE)
  old_wd <- getwd()
  on.exit(setwd(old_wd), add = TRUE)
  setwd(root)

  set.seed(4)
  n_gene <- 20
  n_sample <- 10
  sample_ids <- sprintf("Sample%02d", 1:n_sample)
  annot <- data.frame(
    Kingdom = "Gene", Phylum = "Gene", Class = "Gene", Order = "Gene", Family = "Gene",
    Genus = paste0("Gene", 1:n_gene), Species = sprintf("Gene%04d", 1:n_gene),
    row.names = sprintf("Gene%04d", 1:n_gene), check.names = FALSE
  )
  counts <- matrix(stats::rpois(n_gene * n_sample, lambda = 200), n_gene, n_sample, dimnames = list(rownames(annot), sample_ids))
  counts[1:3, ] <- counts[1:3, ] + 1500L
  write.csv(cbind(annot, as.data.frame(counts, check.names = FALSE)), "1_out/RNAseq_expression_table.csv", row.names = TRUE)
  write.csv(data.frame(Sample = sample_ids, RawReads = 30000, FilteredReads = 25000), "1_out/QC_summary.csv", row.names = FALSE)
  write.csv(
    data.frame(row.names = sample_ids, TreatmentGroup = rep(c("GroupA", "GroupB"), each = n_sample / 2), Sex = rep(c("F", "M"), n_sample / 2), check.names = FALSE),
    "3_map/MAPPING.csv", row.names = TRUE
  )
  invisible(root)
}
