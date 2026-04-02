##############################################################################
#  CCM Phase 2 — Trajectory Deep Analysis
#  Author  : Heekuk Park  (generated 2026-03-29, revised 2026-03-30)
#  Purpose : Paper-inspired analyses NOT yet in 20260321 main script
#            Based on Cell 2025 (Usyk) & Comm.Bio 2025 (Banila) methods
#
#  Requires (already loaded from 20260321 main script):
#    ps2         phyloseq, all samples, full sample_data
#    ps2.traj    phyloseq, trajectory label attached to all timepoints
#    traj_meta   data.frame: StudyID, trajectory, studygrp2
#    adiv.traj   alpha diversity data.frame with trajectory
#    orders      character vector of factor level order
#    cols2       color palette
#    project     "CCM2_16S"
#
#  New analyses:
#    1.  Trajectory count/proportion — stacked bar
#    2.  Sankey / alluvial — CST transition V1->V2->V3 by trajectory
#    3.  Hierarchical clustering heatmap — top species x V1 samples
#    4.  Spearman co-occurrence heatmap — significant bacteria pairs
#    5.  MRS (Microbial Risk Score) with covariate — elastic net
#    6.  Logistic regression — OR table + forest plot (Cell Table 2)
#    7.  JSD longitudinal change — compositional stability
#    8.  Dysbiosis score — CST-based molBV proxy (Cell Fig.2)
#    9.  Combined summary figure (patchwork)
#    10. ConDA-dist — trajectory pairwise DA
#    10-1. ConDA-dist — trajectory pairwise DA with Bray
#    11. Beta diversity — trajectory PCoA at V1/V2
##############################################################################
# rm(list=ls())
source("~/Dropbox/04_scripts/R_source/20260326_git_package_upload_helper.R")
git_commit_and_push("Gotools", document = F, commit_message = "Update: v1.2")
git_commit_and_push("intoASV",document = F, commit_message = "Update: v1.2")
git_commit_and_push("ConDAdist",document = F, commit_message = "Update: v0.9")
git_commit_and_push("dotfiles", commit_message = "Update dotfiles")

# library(devtools)
# devtools::install_github("bbagy/Gotools", force = TRUE)
# devtools::install_github("bbagy/intoASV", force = TRUE)
# devtools::install_github("bbagy/ConDAdist", force = TRUE)


suppressPackageStartupMessages({
  library(circlize)        # colorRamp2
  library(broom)           # tidy() for logistic regression
  library(glue)
  library("Gotools")
  Gotool_dependency()
})


# plyr이 로드되어 있으면 dplyr 함수가 마스킹됨 → 명시적으로 우선순위 지정
if ("plyr" %in% loadedNamespaces()) {
  filter  <- dplyr::filter
  select  <- dplyr::select
  mutate  <- dplyr::mutate
  arrange <- dplyr::arrange
  count   <- dplyr::count
  rename  <- dplyr::rename
}

# ── Shared palettes ───────────────────────────────────────────────────────────
traj_cols <- c(
  "WNL_WNL"   = "#2166AC",
  "CIN2_WNL"  = "#74ADD1",
  "WNL_CIN2"  = "#F46D43",
  "CIN2_CIN2" = "#D73027"
)


traj_labels <- c(
  "WNL_WNL"   = "WNL -> WNL",
  "CIN2_WNL"  = "CIN2+ -> WNL",
  "WNL_CIN2"  = "WNL -> CIN2+",
  "CIN2_CIN2" = "CIN2+ -> CIN2+"
)

grp_labels <- c(
  "Grp1" = "Grp1 (Treated)",
  "Grp2" = "Grp2 (Untreated)"
)

cst_cols <- c(
  "I"    = "#2166AC",
  "II"   = "#74ADD1",
  "III"  = "#FEE090",
  "IV-A" = "#FDBB84",
  "IV-B" = "#E31A1C",
  "IV-C" = "#990000",
  "V"    = "#A6DBA0"
)

#------------------------------#
#----         Input         ---#
#------------------------------#
project<-"CCM2_16S"

currentwd <- "~/Columbia University Irving Medical Center/Uhlemann Lab - Microbiome Core/1_Projects/MAPS_Uhl_All/MAPS_Uhl_CCM2/2_Analysis/16S/20250825_SCrub_merged/"

setwd(sprintf("%s",currentwd))

# Go_blastASVs(project, asvsTable="1_out/CCM2_16S.250826.psTotab.asvTable.csv", blastDB="/Users/heekukpark/DB/blastDB/16S_ribosomal_RNA")
ps <- Go_tabTops(csv = "1_out/CCM2_16S.updated_sequences_with_blast_results.250826_cleaned.csv",project);ps

sampledata <- read.csv("3_map/260204.CCM_Phase2.mapping.csv",row.names=1,check.names=F);head(sampledata)


# read tree file
files <- list.files(path = sprintf("%s/1_out/", currentwd), recursive = TRUE, full.names = TRUE)
nwk_file <- grep(".*\\.seqs\\.fna_tree/exported-tree/tree\\.nwk$", files, value = TRUE)
tree <- read_tree(nwk_file, errorIfNULL = T)
ps1 <- merge_phyloseq(ps, sample_data(data.frame(sampledata)), phy_tree(tree));ps1

#===== Rarefaction and filter
cutoff <- 5000 # 7000 5000 2500
ps1.prune = prune_samples(sample_sums(ps1) > cutoff, ps1);ps1.prune

# p <- Go_rare(ps1.prune, step = 1000, color = "Timepoint", xlimit = 10000, label = NULL, se = FALSE) + theme(legend.position="right")
#dir <- Go_path(project, pdf="yes", table="no", path=NULL)

# ggsave(sprintf("%s/1_rarefaction2.%s.%s.pdf", dir$pdf, project ,format(Sys.Date(), "%y%m%d")), plot = p, device = "pdf", width = 8, height = 5, limitsize = FALSE)

Go_SeqLengths(ps1.prune)
ps1.size <- Go_SeqLengths(ps1.prune, from = 397, to= 450)

ps2 <- Go_filter(ps1.size, cutoff = 0.00005)# 1: 0.00005 (713) 2: 0.0001 (496)


# ── Prerequisite check: ps2.traj & traj_meta ─────────────────────────────────
# ps2.traj is built in 20260321 main script.
# If not found, recreate here from ps2 + pattern CSV files.



if (!exists("traj_meta")) {
  cat("  [prereq] traj_meta not found — building from pattern CSV...\n")

  make_trajectory <- function(df, grp_label) {
    df |>
      dplyr::mutate(
        V2c = ifelse(V2 %in% c("WNL", "CIN2+"), V2, NA_character_),
        V3c = ifelse(V3 %in% c("WNL", "CIN2+"), V3, NA_character_),
        trajectory = dplyr::case_when(
          V2c == "WNL"   & V3c == "WNL"   ~ "WNL_WNL",
          V2c == "WNL"   & V3c == "CIN2+" ~ "WNL_CIN2",
          V2c == "CIN2+" & V3c == "WNL"   ~ "CIN2_WNL",
          V2c == "CIN2+" & V3c == "CIN2+" ~ "CIN2_CIN2",
          TRUE ~ NA_character_
        ),
        studygrp2 = grp_label
      ) |>
      dplyr::filter(!is.na(trajectory)) |>
      dplyr::select(StudyID, trajectory, studygrp2)
  }

  pat_grp1 <- read.csv(
    "CCM2_16S_260302/table/pattern/pattern_wide.CINC.CCM2_16S.5000_1.260222.csv",
    stringsAsFactors = FALSE
  )
  pat_grp2 <- read.csv(
    "CCM2_16S_260302/table/pattern/pattern_wide.CINC.CCM2_16S.5000_2.260222.csv",
    stringsAsFactors = FALSE
  )
  traj_meta <- dplyr::bind_rows(
    make_trajectory(pat_grp1, "Grp1"),
    make_trajectory(pat_grp2, "Grp2")
  )
  cat("  [prereq] traj_meta built: ", nrow(traj_meta), " subjects\n")
}

if (!exists("ps2.traj")) {
  cat("  [prereq] ps2.traj not found — building from ps2 + traj_meta...\n")
  if (!exists("ps2")) stop("ps2 not found. Run the main script first.")

  # distinct() prevents duplicated rows when StudyID maps to multiple rows
  traj_key <- traj_meta |>
    dplyr::select(StudyID, trajectory, studygrp2) |>
    dplyr::distinct(StudyID, .keep_all = TRUE) |>
    dplyr::mutate(trajectory = as.character(trajectory))

  sd_traj_full <- data.frame(
    sample_data(ps2), stringsAsFactors = FALSE
  ) |>
    tibble::rownames_to_column(".SampleID") |>
    dplyr::select(-dplyr::any_of(c("trajectory", "studygrp2"))) |>
    dplyr::inner_join(traj_key, by = "StudyID") |>
    tibble::column_to_rownames(".SampleID")

  ps2.traj <- prune_samples(rownames(sd_traj_full), ps2)
  sample_data(ps2.traj) <- sample_data(sd_traj_full)
  cat("  [prereq] ps2.traj built:", nsamples(ps2.traj), "samples\n")
}

# Diagnostic: verify trajectory column is a plain character vector
sd_check <- data.frame(sample_data(ps2.traj), stringsAsFactors = FALSE)
cat("  trajectory class:", class(sd_check$trajectory), "\n")
cat("  Trajectory x Timepoint counts:\n")
print(table(sd_check$Timepoint, sd_check$trajectory, useNA = "ifany"))

# ── Output directories ────────────────────────────────────────────────────────
out_dirs <- Go_path(project, pdf = "yes", table = "yes", path = NULL)
dir_pdf <- out_dirs$pdf
dir_tab <- out_dirs$tab

today <- format(Sys.Date(), "%y%m%d")

make_path <- function(subdir, prefix, ext = "pdf") {
  file.path(subdir, glue("{prefix}.{project}.{today}.{ext}"))
}


##############################################################################
#  1. TRAJECTORY COUNT — stacked bar
##############################################################################
cat("\n=== [1] Trajectory distribution bar chart ===\n")

traj_summary <- as.data.frame(traj_meta) |>
  dplyr::filter(studygrp2 %in% c("Grp1", "Grp2")) |>
  dplyr::mutate(
    trajectory = factor(trajectory, levels = names(traj_labels)),
    studygrp2  = factor(studygrp2,  levels = names(grp_labels))
  ) |>
  dplyr::count(studygrp2, trajectory) |>
  dplyr::group_by(studygrp2) |>
  dplyr::mutate(pct = n / sum(n) * 100, total = sum(n)) |>
  dplyr::ungroup() |>
  dplyr::mutate(traj_label = traj_labels[as.character(trajectory)])

p_traj_bar <- ggplot(
  traj_summary,
  aes(x = studygrp2, y = pct, fill = trajectory)
) +
  geom_bar(
    stat = "identity", width = 0.6,
    color = "white", linewidth = 0.3
  ) +
  geom_text(
    aes(label = paste0(n, "\n(", round(pct, 1), "%)")),
    position = position_stack(vjust = 0.5),
    size = 3.2, color = "white", fontface = "bold"
  ) +
  scale_fill_manual(
    values = traj_cols, labels = traj_labels,
    name = "Trajectory (V2->V3)"
  ) +
  scale_x_discrete(labels = grp_labels) +
  scale_y_continuous(
    labels = function(x) paste0(x, "%"),
    expand = c(0, 0)
  ) +
  labs(
    title    = "Outcome Trajectory Distribution by Group",
    subtitle = "V2 -> V3 CINC outcome pattern",
    x = NULL, y = "Proportion (%)"
  ) +
  theme_classic(base_size = 13) +
  theme(
    legend.position = "right",
    plot.title      = element_text(face = "bold", hjust = 0.2),
    plot.title.position = "plot",
    axis.text.x     = element_text(size = 12, angle = 20, vjust = 0.5, hjust = 0.5)
  )

ggsave(make_path(dir_pdf, "1_trajectory_bar"), p_traj_bar, width = 5, height = 4)
write.csv(traj_summary, make_path(dir_tab, "1_trajectory_counts", "csv"),row.names = FALSE)
cat("  -> saved\n")


##############################################################################
#  2. SANKEY / ALLUVIAL — CST V1 -> V2 -> V3 by trajectory
##############################################################################
cat("\n=== [2] Sankey: CST V1->V2->V3 by trajectory ===\n")


map_full <- data.frame(sample_data(ps2.traj), stringsAsFactors = FALSE)

cst_wide <- map_full |>
  dplyr::filter(Timepoint %in% c("V1", "V2", "V3"), !is.na(simpleCST), !is.na(trajectory), studygrp2 %in% c("Grp1", "Grp2")) |>
  dplyr::select(StudyID, Timepoint, simpleCST, studygrp2, trajectory) |>
  dplyr::distinct(StudyID, Timepoint, .keep_all = TRUE) |>
  tidyr::pivot_wider(names_from = Timepoint, values_from = simpleCST, values_fill = NA_character_, values_fn = dplyr::first) |>
  dplyr::filter(!is.na(V1), !is.na(V2), !is.na(V3))

all_cst <- c("I", "II", "III", "IV", "V")
all_cst <- intersect(all_cst, unique(c(as.character(cst_wide$V1), as.character(cst_wide$V2), as.character(cst_wide$V3))))
write.csv(cst_wide, make_path(dir_tab, "2_cst_transition_wide", "csv"), row.names = FALSE)

cst_wide <- cst_wide |>
  dplyr::mutate(trajectory = factor(trajectory, levels = names(traj_labels)), studygrp2 = factor(studygrp2, levels = names(grp_labels)), V1 = factor(V1, levels = all_cst), V2 = factor(V2, levels = all_cst), V3 = factor(V3, levels = all_cst))

source("/Users/heekukpark/Dropbox/04_scripts/R_source/Gotools/R/Go_alluvialplot_V2.R")
for (grp in c("Grp1", "Grp2")) Go_alluvialplot(
  project = project, data = droplevels(dplyr::filter(cst_wide, studygrp2 == grp)), mode = "transition", axes = c("V1", "V2", "V3"),
  fill_var = "trajectory", orders = all_cst,
  name = paste0("2_CST_transition_by_trajectory_", grp),
  height = 5, width = 6, mycol = traj_cols
)
cat("  -> saved\n")



unique(cst_wide$trajectory)

##############################################################################
#  3. HIERARCHICAL CLUSTERING HEATMAP
#     Top 37 species x V1 samples, sidebar: Group + Trajectory + CST
#     (Comm.Bio 2025 Fig.3 style)
##############################################################################
cat("\n=== [3] Hierarchical heatmap: top species x V1 ===\n")

ps_v1_traj <- subset_samples(ps2.traj, Timepoint == "V1" & studygrp2 %in% c("Grp1", "Grp2") &!is.na(trajectory))

tryCatch({

  # ── Diagnose available ranks & samples ──────────────────────────────────
  tax_ranks_avail <- colnames(tax_table(ps_v1_traj))
  cat("  tax ranks  :", paste(tax_ranks_avail, collapse = ", "), "\n")
  cat("  n samples  :", nsamples(ps_v1_traj), "\n")
  cat("  n taxa     :", ntaxa(ps_v1_traj), "\n")

  if (nsamples(ps_v1_traj) == 0 || ntaxa(ps_v1_traj) == 0)
    stop("ps_v1_traj is empty")

  # ── Choose finest rank available ─────────────────────────────────────────
  rank_pref <- c("Species", "Genus", "Family", "Order", "Class")
  glom_rank <- rank_pref[rank_pref %in% tax_ranks_avail][1]
  cat("  glom rank  :", glom_rank, "\n")

  ps_v1_sp <- tax_glom(ps_v1_traj, taxrank = glom_rank, NArm = FALSE)
  cat("  taxa/glom  :", ntaxa(ps_v1_sp), "\n")

  if (ntaxa(ps_v1_sp) == 0) stop("No taxa after tax_glom")

  ps_v1_rel <- transform_sample_counts(
    ps_v1_sp, function(x) x / sum(x) * 100
  )

  n_top_ht <- min(37L, ntaxa(ps_v1_rel))
  top_taxa  <- names(
    sort(rowMeans(otu_table(ps_v1_rel)), decreasing = TRUE)
  )[seq_len(n_top_ht)]

  ps_v1_top <- prune_taxa(top_taxa, ps_v1_rel)

  otu_mat <- as.matrix(otu_table(ps_v1_top))
  if (!taxa_are_rows(ps_v1_top)) otu_mat <- t(otu_mat)

  tax_df     <- data.frame(tax_table(ps_v1_top), stringsAsFactors = FALSE)
  label_rank <- rank_pref[rank_pref %in% colnames(tax_df)][1]
  genus_rank <- if ("Genus" %in% colnames(tax_df)) "Genus" else label_rank
  rownames(otu_mat) <- ifelse(
    !is.na(tax_df[[label_rank]]),
    paste(tax_df[[genus_rank]], tax_df[[label_rank]]),
    paste(tax_df[[genus_rank]], "sp.")
  )

  sd_v1 <- data.frame(sample_data(ps_v1_top), stringsAsFactors = FALSE)
  sd_v1 <- sd_v1[colnames(otu_mat), , drop = FALSE]

  traj_ann_col <- traj_cols[
    intersect(names(traj_cols), unique(na.omit(sd_v1$trajectory)))
  ]
  cst_ann_col <- cst_cols[
    intersect(names(cst_cols), unique(na.omit(sd_v1$simpleCST)))
  ]
  grp_ann_col <- c("Grp1" = "#E41A1C", "Grp2" = "#377EB8")

  ha <- HeatmapAnnotation(
    Group      = sd_v1$studygrp2,
    Trajectory = sd_v1$trajectory,
    CST        = sd_v1$simpleCST,
    col = list(
      Group      = grp_ann_col,
      Trajectory = traj_ann_col,
      CST        = cst_ann_col
    )
  )

  col_fun <- colorRamp2(
    c(0, 1, 10, 50, 100),
    c("white", "#FEE5D9", "#FC9272", "#DE2D26", "#67000D")
  )

  ht <- Heatmap(
    otu_mat,
    name               = "Rel. Abund. (%)",
    col                = col_fun,
    top_annotation     = ha,
    cluster_rows       = TRUE,
    cluster_columns    = TRUE,
    clustering_distance_rows    = "euclidean",
    clustering_distance_columns = "euclidean",
    clustering_method_rows      = "ward.D2",
    clustering_method_columns   = "ward.D2",
    show_column_names  = FALSE,
    row_names_gp       = gpar(fontsize = 8),
    column_title       = glue(
      "Top {n_top_ht} {glom_rank} x V1 (n={ncol(otu_mat)})"
    ),
    column_title_gp = gpar(fontsize = 12, fontface = "bold")
  )

  pdf(make_path(dir_pdf, "3_heatmap_V1_species"), width = 14, height = 10)
  draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
  dev.off()
  cat("  -> saved\n")

}, error = function(e) {
  message("  [Fig 3] Heatmap skipped: ", conditionMessage(e))
})


##############################################################################
#  4. SPEARMAN CO-OCCURRENCE HEATMAP  (Cell 2025 Fig.3B style)
##############################################################################
cat("\n=== [4] Spearman co-occurrence heatmap ===\n")

compute_cor_heatmap <- function(ps_sub, grp_name, n_top = 20) {
  tryCatch({
    if (nsamples(ps_sub) < 10) {
      message("  [Fig 4] ", grp_name, ": n<10, skip"); return(invisible(NULL))
    }

    # Use Genus if Species not available
    rank_pref2 <- c("Species", "Genus", "Family", "Order", "Class")
    avail2     <- colnames(tax_table(ps_sub))
    gr         <- rank_pref2[rank_pref2 %in% avail2][1]

    ps_sp  <- tax_glom(ps_sub, taxrank = gr, NArm = FALSE)
    ps_rel <- transform_sample_counts(ps_sp, function(x) x / sum(x) * 100)

    otu_m <- as.matrix(otu_table(ps_rel))
    if (!taxa_are_rows(ps_rel)) otu_m <- t(otu_m)

    # Remove all-zero rows before correlation
    otu_m   <- otu_m[rowSums(otu_m) > 0, , drop = FALSE]
    n_sel   <- min(n_top, nrow(otu_m))
    if (n_sel < 3) {
      message("  [Fig 4] ", grp_name, ": <3 taxa, skip"); return(invisible(NULL))
    }
    top_n   <- names(
      sort(rowMeans(otu_m), decreasing = TRUE)
    )[seq_len(n_sel)]
    otu_top <- otu_m[top_n, , drop = FALSE]

    tx <- data.frame(tax_table(ps_rel)[top_n, ], stringsAsFactors = FALSE)
    gr2 <- if ("Genus" %in% colnames(tx)) "Genus" else gr
    rownames(otu_top) <- ifelse(
      !is.na(tx[[gr]]),
      paste(tx[[gr2]], tx[[gr]]),
      paste(tx[[gr2]], "sp.")
    )

    cor_mat  <- cor(t(otu_top), method = "spearman", use = "pairwise.complete.obs")
    n_row    <- nrow(cor_mat)
    pval_mat <- matrix(NA_real_, n_row, n_row, dimnames = dimnames(cor_mat))

    for (i in seq_len(n_row)) {
      for (j in seq_len(n_row)) {
        if (i == j) next
        xi <- as.numeric(otu_top[i, ])
        xj <- as.numeric(otu_top[j, ])
        # Remove positions where either vector has NA
        ok <- !is.na(xi) & !is.na(xj)
        if (sum(ok) < 4) next
        tst <- tryCatch(
          cor.test(xi[ok], xj[ok], method = "spearman", exact = FALSE),
          error = function(e) NULL
        )
        if (!is.null(tst)) pval_mat[i, j] <- tst$p.value
      }
    }

    sig_mask <- pval_mat < 0.05
    sig_mask[is.na(sig_mask)] <- FALSE

    cell_fn <- function(j, i, x, y, width, height, fill) {
      if (!is.na(sig_mask[i, j]) && sig_mask[i, j]) {
        grid.text(
          sprintf("%.2f", cor_mat[i, j]), x, y,
          gp = gpar(fontsize = 7, col = "white")
        )
      }
    }

    col_cor <- colorRamp2(
      c(-1, -0.5, 0, 0.5, 1),
      c("#053061", "#4393C3", "white", "#D6604D", "#67001F")
    )

    ht_cor <- Heatmap(
      cor_mat,
      name            = "Spearman r",
      col             = col_cor,
      cell_fun        = cell_fn,
      cluster_rows    = TRUE,
      cluster_columns = TRUE,
      row_names_gp    = gpar(fontsize = 8),
      column_names_gp = gpar(fontsize = 8),
      column_title    = glue(
        "{grp_name}: Spearman Co-occurrence (V1, p<0.05 labeled)"
      ),
      column_title_gp = gpar(fontsize = 11, fontface = "bold")
    )

    pdf(make_path(dir_pdf, glue("4_cooccurrence_{grp_name}")),
        width = 9, height = 8)
    draw(ht_cor)
    dev.off()
    cat("  ->", grp_name, "saved\n")
    invisible(NULL)

  }, error = function(e) {
    message("  [Fig 4] ", grp_name, " skipped: ", conditionMessage(e))
  })
}

for (g in c("Grp1", "Grp2")) {
  ps_g_v1 <- subset_samples(
    ps2.traj,
    Timepoint == "V1" & studygrp2 == g & !is.na(trajectory)
  )
  compute_cor_heatmap(ps_g_v1, g, n_top = 20)
}
cat("  -> saved\n")


##############################################################################
#  5. MRS — Microbial Risk Score with covariate  (elastic net)
##############################################################################
cat("\n=== [5] MRS — Microbial Risk Score with covariate ===\n")

source("/Users/heekukpark/Dropbox/04_scripts/R_source/Gotools/R/Go_MRS_fit.R")
source("/Users/heekukpark/Dropbox/04_scripts/R_source/Gotools/R/Go_MRS_plot.R")

ps_v1_traj2 <- subset_samples(
  ps2.traj, Timepoint == "V1" & !is.na(trajectory)
)

mrs_covars <- c("simpleCST", "hivstatus")

format_mrs_title <- function(grp_name, comp_label) {
  title_map <- c(
    "failure_vs_success_cov"      = "Success vs treatment failure",
    "recurrence_vs_success_cov"   = "Success vs recurrence",
    "persistent_vs_clearance_cov" = "Clearance vs persistent",
    "relapse_vs_clearance_cov"    = "Clearance vs relapse",
    "WNL_WNL_vs_HPVneg_cov"       = "HPV- vs clearance",
    "WNL_CIN2_vs_HPVneg_cov"      = "HPV- vs recurrence",
    "CIN2_CIN2_vs_HPVneg_cov"     = "HPV- vs persistent/failure"
  )
  title_txt <- title_map[[comp_label]]
  if (is.null(title_txt) || is.na(title_txt) || !nzchar(title_txt)) {
    title_txt <- gsub("_cov$", "", comp_label)
  }
  glue("{grp_name}: {title_txt}")
}

run_mrs_cov <- function(ps_in, grp_name, pos_class, neg_class, comp_label,
                        plot_title, mycol, seed_val = 1234) {
  sd_sub2 <- data.frame(sample_data(ps_in), stringsAsFactors = FALSE)
  n_sub <- nrow(sd_sub2)
  cat(glue("  [MRS-covar] {grp_name} {comp_label}: {neg_class} vs {pos_class}\n"))
  cat(glue("    samples after filter: {n_sub}\n"))
  if (n_sub < 15) {
    message(glue("  {grp_name} {comp_label}: n={n_sub} (<15), skip"))
    return(invisible(NULL))
  }
  rownames(sd_sub2) <- phyloseq::sample_names(ps_in)
  if (!identical(rownames(sd_sub2), phyloseq::sample_names(ps_in))) {
    stop("sample_data rownames do not match ps_in sample_names.")
  }
  sample_data(ps_in) <- sample_data(sd_sub2)
  fit_mrs <- Go_MRS_fit(psIN = ps_in, outcome = "trajectory_mrs", taxrank = "ASV", top_n = 100, meta_vars = intersect(mrs_covars, colnames(sd_sub2)), validation = "oof", seed = seed_val, project = project, name = glue("5_MRS_{grp_name}_{comp_label}"))
  Go_MRS_plot(fit = fit_mrs, title = plot_title, style = "paper", project = project, name = glue("5_MRS_{grp_name}_{comp_label}"), order = c(neg_class, pos_class), mycol = mycol)
  invisible(fit_mrs)
}

traj_key_v1 <- traj_meta |>
  dplyr::distinct(StudyID, .keep_all = TRUE) |>
  dplyr::mutate(trajectory = as.character(trajectory))

sd_v1_grp3ref <- data.frame(sample_data(ps2), stringsAsFactors = FALSE) |>
  tibble::rownames_to_column(".SampleID") |>
  dplyr::filter(Timepoint == "V1", studygrp2 %in% c("Grp1", "Grp2", "Grp3")) |>
  dplyr::select(-dplyr::any_of("trajectory")) |>
  dplyr::left_join(
    traj_key_v1 |> dplyr::select(StudyID, trajectory),
    by = "StudyID"
  ) |>
  tibble::column_to_rownames(".SampleID")

ps_v1_grp3ref <- prune_samples(rownames(sd_v1_grp3ref), ps2)
sample_data(ps_v1_grp3ref) <- sample_data(sd_v1_grp3ref)
ref_label <- "HPVneg_ref"
ref_col <- "#4DAF4A"

# Grp1: CIN2_CIN2 vs WNL_WNL
grp_val <- "Grp1"; pos_val <- "CIN2_CIN2"; neg_val <- "WNL_WNL"; lbl_val <- "failure_vs_success_cov"
sd_tmp <- data.frame(sample_data(ps_v1_traj2), stringsAsFactors = FALSE)
keep_ids <- rownames(sd_tmp)[sd_tmp$studygrp2 == grp_val & sd_tmp$trajectory %in% c(pos_val, neg_val)]
ps_v1_sub <- prune_samples(keep_ids, ps_v1_traj2)
ps_v1_sub <- subset_samples(ps_v1_sub, trajectory %in% c(pos_val, neg_val))
sd_sub2 <- data.frame(sample_data(ps_v1_sub), stringsAsFactors = FALSE)
sd_sub2$trajectory_mrs <- factor(sd_sub2[["trajectory"]], levels = c(neg_val, pos_val))
sample_data(ps_v1_sub) <- sample_data(sd_sub2)
run_mrs_cov(ps_v1_sub, grp_val, pos_val, neg_val, lbl_val, format_mrs_title(grp_val, lbl_val), c("#2166AC", "#D73027"))

# Grp1: WNL_CIN2 vs WNL_WNL
grp_val <- "Grp1"; pos_val <- "WNL_CIN2"; neg_val <- "WNL_WNL"; lbl_val <- "recurrence_vs_success_cov"
sd_tmp <- data.frame(sample_data(ps_v1_traj2), stringsAsFactors = FALSE)
keep_ids <- rownames(sd_tmp)[sd_tmp$studygrp2 == grp_val & sd_tmp$trajectory %in% c(pos_val, neg_val)]
ps_v1_sub <- prune_samples(keep_ids, ps_v1_traj2)
ps_v1_sub <- subset_samples(ps_v1_sub, trajectory %in% c(pos_val, neg_val))
sd_sub2 <- data.frame(sample_data(ps_v1_sub), stringsAsFactors = FALSE)
sd_sub2$trajectory_mrs <- factor(sd_sub2[["trajectory"]], levels = c(neg_val, pos_val))
sample_data(ps_v1_sub) <- sample_data(sd_sub2)
run_mrs_cov(ps_v1_sub, grp_val, pos_val, neg_val, lbl_val, format_mrs_title(grp_val, lbl_val), c("#2166AC", "#D73027"))

# Grp2: CIN2_CIN2 vs WNL_WNL
grp_val <- "Grp2"; pos_val <- "CIN2_CIN2"; neg_val <- "WNL_WNL"; lbl_val <- "persistent_vs_clearance_cov"
sd_tmp <- data.frame(sample_data(ps_v1_traj2), stringsAsFactors = FALSE)
keep_ids <- rownames(sd_tmp)[sd_tmp$studygrp2 == grp_val & sd_tmp$trajectory %in% c(pos_val, neg_val)]
ps_v1_sub <- prune_samples(keep_ids, ps_v1_traj2)
ps_v1_sub <- subset_samples(ps_v1_sub, trajectory %in% c(pos_val, neg_val))
sd_sub2 <- data.frame(sample_data(ps_v1_sub), stringsAsFactors = FALSE)
sd_sub2$trajectory_mrs <- factor(sd_sub2[["trajectory"]], levels = c(neg_val, pos_val))
sample_data(ps_v1_sub) <- sample_data(sd_sub2)
run_mrs_cov(ps_v1_sub, grp_val, pos_val, neg_val, lbl_val, format_mrs_title(grp_val, lbl_val), c("#2166AC", "#D73027"))

# Grp2: WNL_CIN2 vs WNL_WNL
grp_val <- "Grp2"; pos_val <- "WNL_CIN2"; neg_val <- "WNL_WNL"; lbl_val <- "relapse_vs_clearance_cov"
sd_tmp <- data.frame(sample_data(ps_v1_traj2), stringsAsFactors = FALSE)
keep_ids <- rownames(sd_tmp)[sd_tmp$studygrp2 == grp_val & sd_tmp$trajectory %in% c(pos_val, neg_val)]
ps_v1_sub <- prune_samples(keep_ids, ps_v1_traj2)
ps_v1_sub <- subset_samples(ps_v1_sub, trajectory %in% c(pos_val, neg_val))
sd_sub2 <- data.frame(sample_data(ps_v1_sub), stringsAsFactors = FALSE)
sd_sub2$trajectory_mrs <- factor(sd_sub2[["trajectory"]], levels = c(neg_val, pos_val))
sample_data(ps_v1_sub) <- sample_data(sd_sub2)
run_mrs_cov(ps_v1_sub, grp_val, pos_val, neg_val, lbl_val, format_mrs_title(grp_val, lbl_val), c("#2166AC", "#D73027"))

# Grp1: WNL_WNL vs HPV-
grp_val <- "Grp1"; pos_val <- "WNL_WNL"; neg_val <- ref_label; lbl_val <- "WNL_WNL_vs_HPVneg_cov"
sd_tmp <- data.frame(sample_data(ps_v1_grp3ref), stringsAsFactors = FALSE)
keep_ids <- rownames(sd_tmp)[(sd_tmp$studygrp2 == grp_val & sd_tmp$trajectory == pos_val) | sd_tmp$studygrp2 == "Grp3"]
ps_v1_sub <- prune_samples(keep_ids, ps_v1_grp3ref)
sd_sub2 <- data.frame(sample_data(ps_v1_sub), stringsAsFactors = FALSE)
sd_sub2$trajectory_mrs <- ifelse(sd_sub2$studygrp2 == "Grp3", neg_val, as.character(sd_sub2$trajectory))
sd_sub2$trajectory_mrs <- factor(sd_sub2$trajectory_mrs, levels = c(neg_val, pos_val))
sample_data(ps_v1_sub) <- sample_data(sd_sub2)
run_mrs_cov(ps_v1_sub, grp_val, pos_val, neg_val, lbl_val, format_mrs_title(grp_val, lbl_val), c(ref_col, traj_cols[pos_val]))

# Grp1: WNL_CIN2 vs HPV-
grp_val <- "Grp1"; pos_val <- "WNL_CIN2"; neg_val <- ref_label; lbl_val <- "WNL_CIN2_vs_HPVneg_cov"
sd_tmp <- data.frame(sample_data(ps_v1_grp3ref), stringsAsFactors = FALSE)
keep_ids <- rownames(sd_tmp)[(sd_tmp$studygrp2 == grp_val & sd_tmp$trajectory == pos_val) | sd_tmp$studygrp2 == "Grp3"]
ps_v1_sub <- prune_samples(keep_ids, ps_v1_grp3ref)
sd_sub2 <- data.frame(sample_data(ps_v1_sub), stringsAsFactors = FALSE)
sd_sub2$trajectory_mrs <- ifelse(sd_sub2$studygrp2 == "Grp3", neg_val, as.character(sd_sub2$trajectory))
sd_sub2$trajectory_mrs <- factor(sd_sub2$trajectory_mrs, levels = c(neg_val, pos_val))
sample_data(ps_v1_sub) <- sample_data(sd_sub2)
run_mrs_cov(ps_v1_sub, grp_val, pos_val, neg_val, lbl_val, format_mrs_title(grp_val, lbl_val), c(ref_col, traj_cols[pos_val]))

# Grp1: CIN2_CIN2 vs HPV-
grp_val <- "Grp1"; pos_val <- "CIN2_CIN2"; neg_val <- ref_label; lbl_val <- "CIN2_CIN2_vs_HPVneg_cov"
sd_tmp <- data.frame(sample_data(ps_v1_grp3ref), stringsAsFactors = FALSE)
keep_ids <- rownames(sd_tmp)[(sd_tmp$studygrp2 == grp_val & sd_tmp$trajectory == pos_val) | sd_tmp$studygrp2 == "Grp3"]
ps_v1_sub <- prune_samples(keep_ids, ps_v1_grp3ref)
sd_sub2 <- data.frame(sample_data(ps_v1_sub), stringsAsFactors = FALSE)
sd_sub2$trajectory_mrs <- ifelse(sd_sub2$studygrp2 == "Grp3", neg_val, as.character(sd_sub2$trajectory))
sd_sub2$trajectory_mrs <- factor(sd_sub2$trajectory_mrs, levels = c(neg_val, pos_val))
sample_data(ps_v1_sub) <- sample_data(sd_sub2)
run_mrs_cov(ps_v1_sub, grp_val, pos_val, neg_val, lbl_val, format_mrs_title(grp_val, lbl_val), c(ref_col, traj_cols[pos_val]))

# Grp2: WNL_WNL vs HPV-
grp_val <- "Grp2"; pos_val <- "WNL_WNL"; neg_val <- ref_label; lbl_val <- "WNL_WNL_vs_HPVneg_cov"
sd_tmp <- data.frame(sample_data(ps_v1_grp3ref), stringsAsFactors = FALSE)
keep_ids <- rownames(sd_tmp)[(sd_tmp$studygrp2 == grp_val & sd_tmp$trajectory == pos_val) | sd_tmp$studygrp2 == "Grp3"]
ps_v1_sub <- prune_samples(keep_ids, ps_v1_grp3ref)
sd_sub2 <- data.frame(sample_data(ps_v1_sub), stringsAsFactors = FALSE)
sd_sub2$trajectory_mrs <- ifelse(sd_sub2$studygrp2 == "Grp3", neg_val, as.character(sd_sub2$trajectory))
sd_sub2$trajectory_mrs <- factor(sd_sub2$trajectory_mrs, levels = c(neg_val, pos_val))
sample_data(ps_v1_sub) <- sample_data(sd_sub2)
run_mrs_cov(ps_v1_sub, grp_val, pos_val, neg_val, lbl_val, format_mrs_title(grp_val, lbl_val), c(ref_col, traj_cols[pos_val]))

# Grp2: WNL_CIN2 vs HPV-
grp_val <- "Grp2"; pos_val <- "WNL_CIN2"; neg_val <- ref_label; lbl_val <- "WNL_CIN2_vs_HPVneg_cov"
sd_tmp <- data.frame(sample_data(ps_v1_grp3ref), stringsAsFactors = FALSE)
keep_ids <- rownames(sd_tmp)[(sd_tmp$studygrp2 == grp_val & sd_tmp$trajectory == pos_val) | sd_tmp$studygrp2 == "Grp3"]
ps_v1_sub <- prune_samples(keep_ids, ps_v1_grp3ref)
sd_sub2 <- data.frame(sample_data(ps_v1_sub), stringsAsFactors = FALSE)
sd_sub2$trajectory_mrs <- ifelse(sd_sub2$studygrp2 == "Grp3", neg_val, as.character(sd_sub2$trajectory))
sd_sub2$trajectory_mrs <- factor(sd_sub2$trajectory_mrs, levels = c(neg_val, pos_val))
sample_data(ps_v1_sub) <- sample_data(sd_sub2)
run_mrs_cov(ps_v1_sub, grp_val, pos_val, neg_val, lbl_val, format_mrs_title(grp_val, lbl_val), c(ref_col, traj_cols[pos_val]))

# Grp2: CIN2_CIN2 vs HPV-
grp_val <- "Grp2"; pos_val <- "CIN2_CIN2"; neg_val <- ref_label; lbl_val <- "CIN2_CIN2_vs_HPVneg_cov"
sd_tmp <- data.frame(sample_data(ps_v1_grp3ref), stringsAsFactors = FALSE)
keep_ids <- rownames(sd_tmp)[(sd_tmp$studygrp2 == grp_val & sd_tmp$trajectory == pos_val) | sd_tmp$studygrp2 == "Grp3"]
ps_v1_sub <- prune_samples(keep_ids, ps_v1_grp3ref)
sd_sub2 <- data.frame(sample_data(ps_v1_sub), stringsAsFactors = FALSE)
sd_sub2$trajectory_mrs <- ifelse(sd_sub2$studygrp2 == "Grp3", neg_val, as.character(sd_sub2$trajectory))
sd_sub2$trajectory_mrs <- factor(sd_sub2$trajectory_mrs, levels = c(neg_val, pos_val))
sample_data(ps_v1_sub) <- sample_data(sd_sub2)
run_mrs_cov(ps_v1_sub, grp_val, pos_val, neg_val, lbl_val, format_mrs_title(grp_val, lbl_val), c(ref_col, traj_cols[pos_val]))

cat("  -> saved\n")


##############################################################################
#  6. LOGISTIC REGRESSION — OR TABLE + FOREST PLOT  (Cell Table 2 / Fig.3A)
##############################################################################
cat("\n=== [6] Logistic regression — Forest plot ===\n")

source("/Users/heekukpark/Dropbox/04_scripts/R_source/Gotools/R/Go_OR_fit.R")
source("/Users/heekukpark/Dropbox/04_scripts/R_source/Gotools/R/Go_OR_plot.R")

ps_v1_traj3 <- subset_samples(
  ps2.traj,
  Timepoint == "V1" &
    studygrp2 %in% c("Grp1", "Grp2") &
    !is.na(trajectory)
)

# Grp1: CIN2_CIN2 vs others
grp_val <- "Grp1"; pos_val <- "CIN2_CIN2"; lbl_val <- "treatment_failure"
sd_tmp <- data.frame(sample_data(ps_v1_traj3), stringsAsFactors = FALSE)
keep <- sd_tmp$studygrp2 == grp_val & sd_tmp$trajectory %in% c("WNL_WNL", pos_val)
ps_sub <- prune_samples(keep, ps_v1_traj3)
sd_sub <- data.frame(sample_data(ps_sub), stringsAsFactors = FALSE)
cat(glue("  [OR] {grp_val} {lbl_val}: {pos_val} vs WNL_WNL\n"))
cat(glue("    samples after filter: {nrow(sd_sub)}\n"))
if (nrow(sd_sub) >= 15) {
  fit_or <- Go_OR_fit(psIN = ps_sub, outcome = "trajectory", pos_class = pos_val, taxrank = "Species", covars = c("hivstatus", "simpleCST"), transform = "log_rel", top_n = 15, project = project, name = glue("6_OR_{grp_val}_{lbl_val}"))
  Go_OR_plot(fit = fit_or, title = glue("{grp_val}: {lbl_val}"), order = c("WNL_WNL", pos_val), style = "paper", project = project, name = glue("{grp_val}_{lbl_val}"))
} else message(glue("  {grp_val} {lbl_val}: n={nrow(sd_sub)} (<15), skip"))

# Grp1: WNL_CIN2 vs others
grp_val <- "Grp1"; pos_val <- "WNL_CIN2"; lbl_val <- "recurrence"
sd_tmp <- data.frame(sample_data(ps_v1_traj3), stringsAsFactors = FALSE)
keep <- sd_tmp$studygrp2 == grp_val & sd_tmp$trajectory %in% c("WNL_WNL", pos_val)
ps_sub <- prune_samples(keep, ps_v1_traj3)
sd_sub <- data.frame(sample_data(ps_sub), stringsAsFactors = FALSE)
cat(glue("  [OR] {grp_val} {lbl_val}: {pos_val} vs WNL_WNL\n"))
cat(glue("    samples after filter: {nrow(sd_sub)}\n"))
if (nrow(sd_sub) >= 15) {
  fit_or <- Go_OR_fit(psIN = ps_sub, outcome = "trajectory", pos_class = pos_val, taxrank = "Species", covars = c("hivstatus", "simpleCST"), transform = "log_rel", top_n = 15, project = project, name = glue("6_OR_{grp_val}_{lbl_val}"))
  Go_OR_plot(fit = fit_or, title = glue("{grp_val}: {lbl_val}"), order = c("WNL_WNL", pos_val), style = "paper", project = project, name = glue("{grp_val}_{lbl_val}"))
} else message(glue("  {grp_val} {lbl_val}: n={nrow(sd_sub)} (<15), skip"))

# Grp2: CIN2_CIN2 vs others
grp_val <- "Grp2"; pos_val <- "CIN2_CIN2"; lbl_val <- "persistent"
sd_tmp <- data.frame(sample_data(ps_v1_traj3), stringsAsFactors = FALSE)
keep <- sd_tmp$studygrp2 == grp_val & sd_tmp$trajectory %in% c("WNL_WNL", pos_val)
ps_sub <- prune_samples(keep, ps_v1_traj3)
sd_sub <- data.frame(sample_data(ps_sub), stringsAsFactors = FALSE)
cat(glue("  [OR] {grp_val} {lbl_val}: {pos_val} vs WNL_WNL\n"))
cat(glue("    samples after filter: {nrow(sd_sub)}\n"))
if (nrow(sd_sub) >= 15) {
  fit_or <- Go_OR_fit(psIN = ps_sub, outcome = "trajectory", pos_class = pos_val, taxrank = "Species", covars = c("hivstatus", "simpleCST"), transform = "log_rel", top_n = 15, project = project, name = glue("6_OR_{grp_val}_{lbl_val}"))
  Go_OR_plot(fit = fit_or, title = glue("{grp_val}: {lbl_val}"), order = c("WNL_WNL", pos_val), style = "paper", project = project, name = glue("{grp_val}_{lbl_val}"))
} else message(glue("  {grp_val} {lbl_val}: n={nrow(sd_sub)} (<15), skip"))

# Grp2: WNL_CIN2 vs others
grp_val <- "Grp2"; pos_val <- "WNL_CIN2"; lbl_val <- "relapse"
sd_tmp <- data.frame(sample_data(ps_v1_traj3), stringsAsFactors = FALSE)
keep <- sd_tmp$studygrp2 == grp_val & sd_tmp$trajectory %in% c("WNL_WNL", pos_val)
ps_sub <- prune_samples(keep, ps_v1_traj3)
sd_sub <- data.frame(sample_data(ps_sub), stringsAsFactors = FALSE)
cat(glue("  [OR] {grp_val} {lbl_val}: {pos_val} vs WNL_WNL\n"))
cat(glue("    samples after filter: {nrow(sd_sub)}\n"))
if (nrow(sd_sub) >= 15) {
  fit_or <- Go_OR_fit(psIN = ps_sub, outcome = "trajectory", pos_class = pos_val, taxrank = "Species", covars = c("hivstatus", "simpleCST"), transform = "log_rel", top_n = 15, project = project, name = glue("6_OR_{grp_val}_{lbl_val}"))
  Go_OR_plot(fit = fit_or, title = glue("{grp_val}: {lbl_val}"), order = c("WNL_WNL", pos_val), style = "paper", project = project, name = glue("{grp_val}_{lbl_val}"))
} else message(glue("  {grp_val} {lbl_val}: n={nrow(sd_sub)} (<15), skip"))

cat("  -> saved\n")


##############################################################################
#  7. JSD LONGITUDINAL CHANGE  (Cell 2025 Fig.6B style)
##############################################################################
cat("\n=== [7] JSD longitudinal change ===\n")

jsd_distance <- function(p, q) {
  p <- p / sum(p)
  q <- q / sum(q)
  m <- (p + q) / 2
  h_p <- if (any(p > 0)) -sum(p[p > 0] * log2(p[p > 0])) else 0
  h_q <- if (any(q > 0)) -sum(q[q > 0] * log2(q[q > 0])) else 0
  h_m <- if (any(m > 0)) -sum(m[m > 0] * log2(m[m > 0])) else 0
  sqrt(max(0, h_m - 0.5 * h_p - 0.5 * h_q))
}

ps_sp_jsd  <- tax_glom(ps2.traj, taxrank = "Species", NArm = FALSE)
ps_rel_jsd <- transform_sample_counts(ps_sp_jsd, function(x) x / sum(x))
otu_jsd    <- as.matrix(otu_table(ps_rel_jsd))
if (!taxa_are_rows(ps_rel_jsd)) otu_jsd <- t(otu_jsd)

sd_jsd_all <- data.frame(
  sample_data(ps_rel_jsd), stringsAsFactors = FALSE
) |>
  filter(studygrp2 %in% c("Grp1", "Grp2"), !is.na(trajectory))

jsd_rows <- purrr::map_dfr(unique(sd_jsd_all$StudyID), function(sid) {
  rows <- sd_jsd_all[sd_jsd_all$StudyID == sid, ]
  if (!"V1" %in% rows$Timepoint) return(NULL)
  v1_samp <- rownames(rows)[rows$Timepoint == "V1"]
  meta    <- rows[rows$Timepoint == "V1", ]

  purrr::map_dfr(c("V2", "V3"), function(tp2) {
    if (!tp2 %in% rows$Timepoint) return(NULL)
    s2 <- rownames(rows)[rows$Timepoint == tp2]
    jd <- tryCatch(
      jsd_distance(otu_jsd[, v1_samp], otu_jsd[, s2]),
      error = function(e) NA_real_
    )
    data.frame(
      StudyID    = sid,
      To         = tp2,
      JSD        = jd,
      trajectory = meta$trajectory,
      studygrp2  = meta$studygrp2,
      stringsAsFactors = FALSE
    )
  })
})

jsd_df <- jsd_rows |>
  filter(!is.na(JSD)) |>
  mutate(
    trajectory = factor(trajectory, levels = names(traj_labels)),
    Comparison = paste0("V1 -> ", To)
  )

p_jsd <- ggplot(
  jsd_df,
  aes(x = trajectory, y = JSD, fill = trajectory)
) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 1.2, alpha = 0.4) +
  scale_fill_manual(values = traj_cols, labels = traj_labels, name = NULL) +
  scale_x_discrete(labels = traj_labels) +
  facet_grid(
    Comparison ~ studygrp2,
    labeller = labeller(studygrp2 = grp_labels)
  ) +
  labs(
    title    = "Compositional Stability: JSD from V1 per Trajectory",
    subtitle = "Higher JSD = greater change from baseline",
    x = "Trajectory (V2->V3)", y = "JSD Distance (V1 -> Vx)"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title       = element_text(face = "bold"),
    axis.text.x      = element_text(angle = 30, hjust = 1, size = 8),
    legend.position  = "none",
    strip.background = element_rect(fill = "grey90")
  )

ggsave(make_path(dir_pdf, "7_JSD_longitudinal"), p_jsd, width = 9, height = 6)
write.csv(jsd_df, make_path(dir_tab, "7_JSD_by_trajectory", "csv"),
          row.names = FALSE)
cat("  -> saved\n")


##############################################################################
#  8. DYSBIOSIS SCORE — CST-based molBV proxy  (Cell Fig.2 style)
##############################################################################
cat("\n=== [8] Dysbiosis score (CST-proxy molBV) ===\n")

cst_score_map <- c(
  "I"    = 0,
  "II"   = 0,
  "III"  = 1,
  "V"    = 0.5,
  "IV-A" = 2,
  "IV-B" = 3,
  "IV-C" = 3
)

sd_dysbio <- data.frame(
  sample_data(ps2.traj), stringsAsFactors = FALSE
) |>
  filter(
    studygrp2 %in% c("Grp1", "Grp2"),
    Timepoint %in% c("V1", "V2", "V3"),
    !is.na(trajectory),
    !is.na(CST)
  ) |>
  mutate(
    dysbiosis_score = cst_score_map[CST],
    trajectory      = factor(trajectory, levels = names(traj_labels)),
    Timepoint       = factor(Timepoint,  levels = c("V1", "V2", "V3"))
  )

p_dysbio <- ggplot(
  sd_dysbio,
  aes(
    x     = Timepoint,
    y     = dysbiosis_score,
    color = trajectory,
    group = interaction(trajectory, Timepoint)
  )
) +
  geom_boxplot(
    aes(fill = trajectory),
    alpha         = 0.25,
    outlier.shape = NA,
    position      = position_dodge(0.7),
    width         = 0.5
  ) +
  geom_point(
    position = position_jitterdodge(
      jitter.width = 0.1, dodge.width = 0.7
    ),
    size = 1.2, alpha = 0.5
  ) +
  scale_fill_manual(
    values = traj_cols, labels = traj_labels, name = "Trajectory"
  ) +
  scale_color_manual(
    values = traj_cols, labels = traj_labels, name = "Trajectory"
  ) +
  facet_wrap(
    ~studygrp2, labeller = labeller(studygrp2 = grp_labels)
  ) +
  labs(
    title    = "Dysbiosis Score (CST-proxy) over Time by Trajectory",
    subtitle = "Score: I/II=0, V=0.5, III=1, IV-A=2, IV-B/C=3",
    x = "Visit", y = "Dysbiosis Score"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title       = element_text(face = "bold"),
    legend.position  = "bottom",
    strip.background = element_rect(fill = "grey90")
  ) +
  guides(
    fill  = guide_legend(nrow = 2),
    color = guide_legend(nrow = 2)
  )

ggsave(make_path(dir_pdf, "8_dysbiosis_score_trajectory"),
       p_dysbio, width = 10, height = 6)

dysbio_summary <- sd_dysbio |>
  group_by(studygrp2, Timepoint, trajectory) |>
  summarise(
    n    = n(),
    mean = round(mean(dysbiosis_score, na.rm = TRUE), 2),
    sd   = round(sd(dysbiosis_score,   na.rm = TRUE), 2),
    .groups = "drop"
  )
write.csv(
  dysbio_summary,
  make_path(dir_tab, "8_dysbiosis_summary", "csv"),
  row.names = FALSE
)
cat("  -> saved\n")


##############################################################################
#  9. COMBINED SUMMARY FIGURE  (patchwork)
##############################################################################
cat("\n=== [9] Combined summary figure ===\n")

p_dysbio_v1 <- sd_dysbio |>
  filter(Timepoint == "V1") |>
  ggplot(aes(x = trajectory, y = dysbiosis_score, fill = trajectory)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 1.2, alpha = 0.4) +
  scale_fill_manual(values = traj_cols, labels = traj_labels, name = NULL) +
  scale_x_discrete(labels = traj_labels) +
  facet_wrap(~studygrp2, labeller = labeller(studygrp2 = grp_labels)) +
  labs(title = "B. V1 Dysbiosis Score", x = NULL, y = "Dysbiosis Score") +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x     = element_text(angle = 30, hjust = 1, size = 8),
    legend.position = "none",
    plot.title      = element_text(face = "bold", size = 11)
  )

p_jsd_v2 <- jsd_df |>
  filter(To == "V2") |>
  ggplot(aes(x = trajectory, y = JSD, fill = trajectory)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 1.2, alpha = 0.4) +
  scale_fill_manual(values = traj_cols, labels = traj_labels, name = NULL) +
  scale_x_discrete(labels = traj_labels) +
  facet_wrap(~studygrp2, labeller = labeller(studygrp2 = grp_labels)) +
  labs(
    title = "C. Compositional Change V1 -> V2 (JSD)",
    x = NULL, y = "JSD Distance"
  ) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x     = element_text(angle = 30, hjust = 1, size = 8),
    legend.position = "none",
    plot.title      = element_text(face = "bold", size = 11)
  )

combined <- (p_traj_bar + labs(title = "A. Trajectory Distribution")) /
  (p_dysbio_v1 | p_jsd_v2) +
  plot_layout(heights = c(1, 1.5)) +
  plot_annotation(
    title    = "CCM Phase 2: Trajectory-based Microbiome Analysis",
    subtitle = glue(
      "V1 baseline -> V2/V3 outcome  |  {today}"
    ),
    theme = theme(
      plot.title    = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 11)
    )
  )

ggsave(make_path(dir_pdf, "9_combined_summary"), combined,
       width = 14, height = 12)
cat("  -> saved\n")


##############################################################################
#  10. ConDA-dist — trajectory pairwise DA
##############################################################################
cat("\n=== [10] ConDA-dist by trajectory comparisons ===\n")

if (!exists("Go_ConDaDist")) {
  conda_r_dir <- "/Users/heekukpark/Dropbox/04_scripts/ConDA-dist/ConDA-dist/R"
  if (dir.exists(conda_r_dir)) {
    for (f in list.files(conda_r_dir, pattern = "[.]R$", full.names = TRUE)) {
      source(f)
    }
  } else {
    stop("ConDA-dist R directory not found.")
  }
}

da.col <- c("#74ADD1", "#E31A1C")
traj_comps <- list(
  list(grp = "Grp1", g1 = "WNL_WNL", g2 = "CIN2_CIN2", label = "success_vs_failure"),
  list(grp = "Grp1", g1 = "WNL_WNL", g2 = "WNL_CIN2", label = "success_vs_recurrence"),
  list(grp = "Grp1", g1 = "WNL_WNL", g2 = "CIN2_WNL", label = "success_vs_late_resp"),
  list(grp = "Grp2", g1 = "WNL_WNL", g2 = "CIN2_CIN2", label = "clearance_vs_persistent"),
  list(grp = "Grp2", g1 = "WNL_WNL", g2 = "WNL_CIN2", label = "clearance_vs_relapse"),
  list(grp = "Grp2", g1 = "WNL_WNL", g2 = "CIN2_WNL", label = "clearance_vs_late_clear")
)

conda_result <- NULL
for (da in c("deseq2", "ancombc2")) {
  cat(glue("  [ConDA-dist] method = {da}\n"))

  for (tp in c("V1", "V2")) {
    ps_tp <- subset_samples(ps2.traj, Timepoint == tp & !is.na(trajectory))

    for (tc in traj_comps) {
      ps_sub <- subset_samples(
        ps_tp,
        studygrp2 == tc$grp & trajectory %in% c(tc$g1, tc$g2)
      )
      sd_sub <- data.frame(sample_data(ps_sub), stringsAsFactors = FALSE)
      traj_present <- intersect(c(tc$g1, tc$g2), unique(as.character(sd_sub$trajectory)))
      if (length(traj_present) < 2 || nrow(sd_sub) < 10) {
        message(glue("  skip {da} {tp} {tc$grp} {tc$label}: insufficient samples"))
        next
      }

      nm <- glue("{tc$grp}_{tp}")
      cat(glue("    {tp} {tc$grp} {tc$g1} vs {tc$g2} (n={nrow(sd_sub)})\n"))
      conda_result <- Go_ConDaDist(psIN = ps_sub, project = project, group_var = "trajectory", group_1 = tc$g1, group_2 = tc$g2, covariates = intersect(c("hivstatus", "simpleCST", "Batchs"), colnames(sd_sub)), methods = da, distances = NULL, name = nm)
    }
  }
  if (!is.null(conda_result)) {
    Go_volcanoPlot(project = project, result = conda_result, fc = 0, mycols = da.col, name = NULL, font = 3, height = 5, width = 6)
  }
}
cat("  -> saved\n")


##############################################################################
#  10-1. ConDA-dist — trajectory pairwise DA with Bray
##############################################################################
cat("\n=== [10-1] ConDA-dist by trajectory comparisons with Bray ===\n")

conda_result <- NULL
cat("  [ConDA-dist] methods = deseq2 + ancombc2 | distance = bray\n")

for (tp in c("V1", "V2")) {
  ps_tp <- subset_samples(ps2.traj, Timepoint == tp & !is.na(trajectory))

  for (tc in traj_comps) {
    ps_sub <- subset_samples(
      ps_tp,
      studygrp2 == tc$grp & trajectory %in% c(tc$g1, tc$g2)
    )
    sd_sub <- data.frame(sample_data(ps_sub), stringsAsFactors = FALSE)
    traj_present <- intersect(c(tc$g1, tc$g2), unique(as.character(sd_sub$trajectory)))
    if (length(traj_present) < 2 || nrow(sd_sub) < 10) {
      message(glue("  skip {tp} {tc$grp} {tc$label}: insufficient samples"))
      next
    }

    nm <- glue("{tc$grp}_{tp}")
    cat(glue("    {tp} {tc$grp} {tc$g1} vs {tc$g2} (n={nrow(sd_sub)})\n"))
    conda_result <- Go_ConDaDist(psIN = ps_sub, project = project, group_var = "trajectory", group_1 = tc$g1, group_2 = tc$g2, covariates = intersect(c("hivstatus", "simpleCST", "Batchs"), colnames(sd_sub)), methods = c("deseq2", "ancombc2"), distances = "bray", name = nm)
  }
}
if (!is.null(conda_result)) {
  Go_volcanoPlot(project = project, result = conda_result, fc = 0, mycols = da.col, name = NULL, font = 3, height = 5, width = 6)
}
cat("  -> saved\n")


##############################################################################
#  11. Beta diversity — trajectory PCoA at V1/V2
##############################################################################
cat("\n=== [11] Beta diversity by trajectory at V1/V2 ===\n")

for (tp in c("V1", "V2")) {
  ps2.traj.tp <- subset_samples(ps2.traj, Timepoint == tp & !is.na(trajectory))

  for (g in c("Grp1", "Grp2")) {
    ps2.traj.tp.sel <- subset_samples(ps2.traj.tp, studygrp2 == g)
    if (nsamples(ps2.traj.tp.sel) < 5) {
      message(glue("  skip bdiv {g} {tp}: insufficient samples"))
      next
    }

    Go_bdivPM(psIN = ps2.traj.tp.sel, cate.vars = "trajectory", project = project, orders = orders, distance_metrics = c("bray", "wunifrac"), cate.conf = c("hivstatus", "simpleCST", "Batchs"), mycols = cols2, facet = NULL, name = sprintf("%s_%s_trajectory", g, tp), height = 5, width = 5)
  }
}
cat("  -> saved\n")


##############################################################################
#  Done
##############################################################################
