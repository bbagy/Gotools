##############################################################################
#  CCM Phase 2 — Trajectory Deep Analysis
#  Author  : claude-sonnet-4-6  (generated 2026-03-29, revised 2026-03-30)
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
#    5.  MRS (Microbial Risk Score) — weighted microbial predictor
#    6.  Logistic regression — OR table + forest plot (Cell Table 2)
#    7.  JSD longitudinal change — compositional stability
#    8.  Dysbiosis score — CST-based molBV proxy (Cell Fig.2)
#    9.  Combined summary figure (patchwork)
##############################################################################

suppressPackageStartupMessages({
  library(phyloseq)
  library(tidyverse)
  library(ggplot2)
  library(patchwork)
  library(ggalluvial)      # Sankey / alluvial
  library(ComplexHeatmap)  # hierarchical heatmap
  library(circlize)        # colorRamp2
  library(grid)            # grid.text, gpar
  library(vegan)           # JSD
  library(broom)           # tidy() for logistic regression
  library(glue)
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
  "Grp2" = "Grp2 (HPV+, Untreated)"
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
dir_pdf <- "CCM2_16S_traj_deep/pdf"
dir_tab <- "CCM2_16S_traj_deep/table"
dir.create(dir_pdf, showWarnings = FALSE, recursive = TRUE)
dir.create(dir_tab, showWarnings = FALSE, recursive = TRUE)

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
    plot.title      = element_text(face = "bold"),
    axis.text.x     = element_text(size = 12)
  )

ggsave(make_path(dir_pdf, "1_trajectory_bar"), p_traj_bar,
       width = 6, height = 5)
write.csv(traj_summary, make_path(dir_tab, "1_trajectory_counts", "csv"),
          row.names = FALSE)
cat("  -> saved\n")


##############################################################################
#  2. SANKEY / ALLUVIAL — CST V1 -> V2 -> V3 by trajectory
##############################################################################
cat("\n=== [2] Sankey: CST V1->V2->V3 by trajectory ===\n")

map_full <- data.frame(sample_data(ps2.traj), stringsAsFactors = FALSE)

cst_wide <- map_full |>
  dplyr::filter(
    Timepoint %in% c("V1", "V2", "V3"),
    !is.na(simpleCST),
    studygrp2 %in% c("Grp1", "Grp2")
  ) |>
  dplyr::select(StudyID, Timepoint, simpleCST, studygrp2, trajectory) |>
  # one row per StudyID x Timepoint (take first if duplicated)
  dplyr::distinct(StudyID, Timepoint, .keep_all = TRUE) |>
  tidyr::pivot_wider(
    names_from  = Timepoint,
    values_from = simpleCST,
    values_fill = NA_character_,
    values_fn   = dplyr::first   # prevents list-column when duplicates exist
  ) |>
  dplyr::filter(!is.na(V1), !is.na(V2), !is.na(V3))

all_cst <- sort(unique(c(cst_wide$V1, cst_wide$V2, cst_wide$V3)))

for (g in c("Grp1", "Grp2")) {
  cst_g <- cst_wide |>
    dplyr::filter(studygrp2 == g) |>
    dplyr::mutate(
      trajectory = factor(trajectory, levels = names(traj_labels)),
      V1 = factor(V1, levels = all_cst),
      V2 = factor(V2, levels = all_cst),
      V3 = factor(V3, levels = all_cst)
    )

  flow_df <- as.data.frame(cst_g) |>
    dplyr::count(trajectory, V1, V2, V3)

  p_sankey <- ggplot(
    flow_df,
    aes(axis1 = V1, axis2 = V2, axis3 = V3, y = n)
  ) +
    geom_alluvium(
      aes(fill = trajectory),
      width = 0.3, alpha = 0.75, knot.pos = 0.4
    ) +
    geom_stratum(
      width = 0.3, fill = "grey90",
      color = "grey50", linewidth = 0.3
    ) +
    geom_text(
      stat = "stratum",
      aes(label = after_stat(stratum)),
      size = 3.2, fontface = "bold"
    ) +
    scale_fill_manual(
      values = traj_cols, labels = traj_labels,
      name = "Trajectory"
    ) +
    scale_x_discrete(
      limits = c("V1", "V2", "V3"),
      labels = c("Baseline (V1)", "Follow-up (V2)", "Follow-up (V3)"),
      expand = c(0.05, 0.05)
    ) +
    labs(
      title    = glue("CST Transition: {grp_labels[g]}"),
      subtitle = "Each stream = subject subset by trajectory",
      x = "Visit", y = "Number of subjects"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title      = element_text(face = "bold"),
      panel.grid      = element_blank(),
      legend.position = "right"
    )

  ggsave(make_path(dir_pdf, glue("2_sankey_CST_{g}")), p_sankey,
         width = 8, height = 6)
}
cat("  -> saved\n")


##############################################################################
#  3. HIERARCHICAL CLUSTERING HEATMAP
#     Top 37 species x V1 samples, sidebar: Group + Trajectory + CST
#     (Comm.Bio 2025 Fig.3 style)
##############################################################################
cat("\n=== [3] Hierarchical heatmap: top species x V1 ===\n")

ps_v1_traj <- subset_samples(
  ps2.traj,
  Timepoint == "V1" &
    studygrp2 %in% c("Grp1", "Grp2") &
    !is.na(trajectory)
)

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
#  5. MRS — Microbial Risk Score
#     (Cell 2025 Fig.3C — "polygenic risk score analog")
##############################################################################
cat("\n=== [5] MRS — Microbial Risk Score ===\n")

annotate_wilcox <- function(grp_vec, val_vec, neg, pos) {
  sub_df <- data.frame(g = grp_vec, v = val_vec)
  sub_df <- sub_df[sub_df$g %in% c(neg, pos), ]
  wt  <- wilcox.test(v ~ g, data = sub_df)
  p   <- wt$p.value
  lbl <- if (p < 0.001) {
    "p < 0.001"
  } else if (p < 0.05) {
    glue("p = {round(p, 3)}")
  } else {
    glue("p = {round(p, 2)} (NS)")
  }
  annotate("text", x = 1.5, y = max(val_vec, na.rm = TRUE) * 0.95,
           label = lbl, size = 4, fontface = "italic")
}

compute_mrs <- function(ps_v1_sub, outcome_col,
                        pos_class, neg_class,
                        grp_name, comp_label) {

  if (!requireNamespace("pROC", quietly = TRUE)) {
    message("  pROC not installed — skipping"); return(invisible(NULL))
  }

  sd_sub <- data.frame(sample_data(ps_v1_sub), stringsAsFactors = FALSE)
  sd_sub <- sd_sub[sd_sub[[outcome_col]] %in% c(pos_class, neg_class), ]
  if (nrow(sd_sub) < 15) {
    message(glue("  {grp_name} {comp_label}: n<15, skip"))
    return(invisible(NULL))
  }

  ps_v1_sub2 <- prune_samples(rownames(sd_sub), ps_v1_sub)
  ps_sp      <- tax_glom(ps_v1_sub2, taxrank = "Species", NArm = FALSE)
  ps_rel     <- transform_sample_counts(
    ps_sp, function(x) log(x / sum(x) * 100 + 1)
  )

  otu_m <- as.matrix(otu_table(ps_rel))
  if (!taxa_are_rows(ps_rel)) otu_m <- t(otu_m)

  n_top2  <- min(20L, nrow(otu_m))
  top20   <- names(sort(rowMeans(otu_m), decreasing = TRUE))[seq_len(n_top2)]
  feat_df <- data.frame(
    t(otu_m[top20, , drop = FALSE]), check.names = TRUE
  )

  sd_sub$y_bin <- as.integer(sd_sub[[outcome_col]] == pos_class)
  feat_df      <- feat_df[rownames(sd_sub), , drop = FALSE]
  feat_df$.y   <- sd_sub$y_bin

  if (requireNamespace("glmnet", quietly = TRUE)) {
    x_mat   <- as.matrix(feat_df[, seq_len(ncol(feat_df) - 1)])
    y_vec   <- feat_df$.y
    cv_fit  <- glmnet::cv.glmnet(
      x_mat, y_vec, family = "binomial", alpha = 0.5, nfolds = 5
    )
    coefs   <- as.matrix(coef(cv_fit, s = "lambda.min"))[-1, , drop = FALSE]
    weights <- coefs[coefs[, 1] != 0, , drop = FALSE]
  } else {
    feat_names <- names(feat_df)[seq_len(ncol(feat_df) - 1)]
    pvals <- vapply(feat_names, function(v) {
      tryCatch(
        coef(summary(
          glm(feat_df$.y ~ feat_df[[v]], family = binomial)
        ))[2, 4],
        error = function(e) 1
      )
    }, numeric(1))
    sel <- feat_names[pvals < 0.1]
    if (length(sel) == 0) {
      message("  No features"); return(invisible(NULL))
    }
    coef_vals <- vapply(sel, function(v) {
      tryCatch(
        coef(glm(feat_df$.y ~ feat_df[[v]], family = binomial))[2],
        error = function(e) 0
      )
    }, numeric(1))
    weights <- matrix(coef_vals, ncol = 1L, dimnames = list(sel, "s1"))
  }

  if (nrow(weights) == 0) {
    message("  No features in MRS"); return(invisible(NULL))
  }

  feat_mat   <- as.matrix(feat_df[, rownames(weights), drop = FALSE])
  mrs        <- as.numeric(feat_mat %*% weights)
  sd_sub$MRS <- mrs

  roc_obj <- pROC::roc(sd_sub$y_bin, sd_sub$MRS, quiet = TRUE)
  auc_val <- round(as.numeric(pROC::auc(roc_obj)), 3)

  grp_vec <- sd_sub[[outcome_col]]
  p_mrs <- ggplot(
    sd_sub,
    aes(
      x    = factor(.data[[outcome_col]], levels = c(neg_class, pos_class)),
      y    = MRS,
      fill = factor(.data[[outcome_col]], levels = c(neg_class, pos_class))
    )
  ) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    geom_jitter(width = 0.15, size = 1.5, alpha = 0.5) +
    annotate_wilcox(grp_vec, mrs, neg_class, pos_class) +
    scale_fill_manual(
      values = setNames(c("#2166AC", "#D73027"), c(neg_class, pos_class)),
      labels = c(neg_class, pos_class), name = NULL
    ) +
    labs(
      title    = glue("{grp_name}: {comp_label}"),
      subtitle = glue("AUC = {auc_val}  |  n = {nrow(sd_sub)}"),
      x = "Trajectory", y = "Microbial Risk Score (MRS)"
    ) +
    theme_classic(base_size = 12) +
    theme(plot.title = element_text(face = "bold"))

  ggsave(make_path(dir_pdf, glue("5_MRS_{grp_name}_{comp_label}")),
         p_mrs, width = 5, height = 5)

  mrs_tab <- sd_sub |>
    select("StudyID", all_of(outcome_col), "MRS") |>
    mutate(Group = grp_name, Comparison = comp_label, AUC = auc_val)
  write.csv(
    mrs_tab,
    make_path(dir_tab, glue("5_MRS_{grp_name}_{comp_label}"), "csv"),
    row.names = FALSE
  )
  invisible(list(roc = roc_obj, auc = auc_val, weights = weights))
}

compute_mrs_gotools <- function(ps_v1_sub, outcome_col,
                                pos_class, neg_class,
                                grp_name, comp_label) {

  source("/Users/heekukpark/Dropbox/04_scripts/R_source/Gotools/R/Go_MRS_fit.R")
  source("/Users/heekukpark/Dropbox/04_scripts/R_source/Gotools/R/Go_MRS_plot.R")

  sd_sub <- data.frame(sample_data(ps_v1_sub), stringsAsFactors = FALSE)
  sd_sub <- sd_sub[sd_sub[[outcome_col]] %in% c(pos_class, neg_class), , drop = FALSE]
  if (nrow(sd_sub) < 15) {
    message(glue("  {grp_name} {comp_label} [Go_MRS_fit]: n<15, skip"))
    return(invisible(NULL))
  }

  ps_v1_sub2 <- prune_samples(rownames(sd_sub), ps_v1_sub)
  sd_sub2 <- data.frame(sample_data(ps_v1_sub2), stringsAsFactors = FALSE)
  sd_sub2$trajectory_mrs <- factor(sd_sub2[[outcome_col]], levels = c(neg_class, pos_class))
  rownames(sd_sub2) <- rownames(data.frame(sample_data(ps_v1_sub2), stringsAsFactors = FALSE))
  sample_data(ps_v1_sub2) <- sample_data(sd_sub2)

  fit_mrs <- Go_MRS_fit(
    psIN = ps_v1_sub2,
    outcome = "trajectory_mrs",
    taxrank = "Species",
    transform = "log_rel",
    top_n = 20,
    alpha = 0.5,
    nfolds = 5,
    validation = "apparent",
    score_scale = "link"
  )

  p_fit <- Go_MRS_plot(
    fit = fit_mrs,
    plot_type = "score",
    title = glue("{grp_name}: {comp_label} [Go_MRS_fit]"),
    style = "paper",
    order = c(neg_class, pos_class),
    mycol = c("#2166AC", "#D73027")
  )
  ggsave(make_path(dir_pdf, glue("5_MRS_GoMRS_{grp_name}_{comp_label}")), p_fit, width = 5, height = 5)

  pred_tab <- fit_mrs$predictions |>
    dplyr::mutate(
      Group = grp_name,
      Comparison = comp_label,
      observed_label = factor(observed, levels = c(0, 1), labels = c(neg_class, pos_class)),
      AUC = fit_mrs$metrics$auc
    )
  write.csv(
    pred_tab,
    make_path(dir_tab, glue("5_MRS_GoMRS_{grp_name}_{comp_label}"), "csv"),
    row.names = FALSE
  )

  invisible(fit_mrs)
}

ps_v1_traj2 <- subset_samples(
  ps2.traj, Timepoint == "V1" & !is.na(trajectory)
)

mrs_comps <- list(
  list(grp = "Grp1", pos = "CIN2_CIN2", neg = "WNL_WNL",
       label = "failure_vs_success"),
  list(grp = "Grp1", pos = "WNL_CIN2",  neg = "WNL_WNL",
       label = "recurrence_vs_success"),
  list(grp = "Grp2", pos = "CIN2_CIN2", neg = "WNL_WNL",
       label = "persistent_vs_clearance"),
  list(grp = "Grp2", pos = "WNL_CIN2",  neg = "WNL_WNL",
       label = "relapse_vs_clearance")
)

invisible(lapply(mrs_comps, function(mc) {
  grp_val <- mc$grp
  pos_val <- mc$pos
  neg_val <- mc$neg
  lbl_val <- mc$label
  # prune_samples with logical vector avoids subset_samples env issue
  sd_tmp  <- data.frame(sample_data(ps_v1_traj2), stringsAsFactors = FALSE)
  keep    <- sd_tmp$studygrp2 == grp_val &
             sd_tmp$trajectory %in% c(pos_val, neg_val)
  ps_sub  <- prune_samples(keep, ps_v1_traj2)
  compute_mrs(ps_sub, "trajectory", pos_val, neg_val, grp_val, lbl_val)
  compute_mrs_gotools(ps_sub, "trajectory", pos_val, neg_val, grp_val, lbl_val)
}))
cat("  -> saved\n")


##############################################################################
#  6. LOGISTIC REGRESSION — OR TABLE + FOREST PLOT  (Cell Table 2 / Fig.3A)
##############################################################################
cat("\n=== [6] Logistic regression — Forest plot ===\n")

run_logistic <- function(ps_v1_sub, outcome_col, pos_class,
                         covars = c("hivstatus", "simpleCST"),
                         grp_name, comp_label) {

  sd_lr <- data.frame(sample_data(ps_v1_sub), stringsAsFactors = FALSE)
  sd_lr <- sd_lr[!is.na(sd_lr[[outcome_col]]), ]
  if (nrow(sd_lr) < 15) return(invisible(NULL))

  sd_lr$y <- as.integer(sd_lr[[outcome_col]] == pos_class)

  ps_sp  <- tax_glom(ps_v1_sub, taxrank = "Species", NArm = FALSE)
  ps_sp  <- prune_samples(rownames(sd_lr), ps_sp)
  ps_rel <- transform_sample_counts(
    ps_sp, function(x) log(x / sum(x) * 100 + 1)
  )
  otu_m <- as.matrix(otu_table(ps_rel))
  if (!taxa_are_rows(ps_rel)) otu_m <- t(otu_m)

  n_top_lr <- min(15L, nrow(otu_m))
  top15    <- names(
    sort(rowMeans(otu_m), decreasing = TRUE)
  )[seq_len(n_top_lr)]
  tx       <- data.frame(tax_table(ps_rel)[top15, ], stringsAsFactors = FALSE)
  taxa_ids <- make.names(
    ifelse(!is.na(tx$Species),
           paste(tx$Genus, tx$Species),
           paste(tx$Genus, "sp."))
  )

  feat_df <- data.frame(
    t(otu_m[top15, , drop = FALSE]), check.names = FALSE
  )
  colnames(feat_df) <- taxa_ids
  feat_df  <- feat_df[rownames(sd_lr), , drop = FALSE]
  avail_cov <- intersect(covars, colnames(sd_lr))

  results <- purrr::map_dfr(taxa_ids, function(taxa) {
    d <- cbind(
      y     = sd_lr$y,
      taxon = feat_df[[taxa]],
      sd_lr[, avail_cov, drop = FALSE]
    )
    fm <- as.formula(
      paste("y ~ taxon +", paste(avail_cov, collapse = " + "))
    )
    tryCatch({
      m   <- glm(fm, data = d, family = binomial)
      res <- broom::tidy(m, exponentiate = TRUE, conf.int = TRUE)
      res |>
        filter(term == "taxon") |>
        mutate(taxa = taxa, OR = estimate, LCI = conf.low, UCI = conf.high)
    }, error = function(e) NULL)
  })

  if (is.null(results) || nrow(results) == 0) return(invisible(NULL))

  results <- results |>
    mutate(
      sig        = p.value < 0.05,
      taxa_clean = stringr::str_replace_all(taxa, "\\.", " ")
    ) |>
    arrange(OR)

  p_forest <- ggplot(
    results,
    aes(y = reorder(taxa_clean, OR), x = OR, color = sig)
  ) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey60") +
    geom_errorbarh(
      aes(xmin = LCI, xmax = UCI), height = 0.25, linewidth = 0.7
    ) +
    geom_point(aes(size = -log10(p.value + 1e-6)), shape = 18) +
    geom_text(
      aes(label = ifelse(sig, sprintf("%.2f", OR), "")),
      hjust = -0.3, size = 3, fontface = "italic"
    ) +
    scale_color_manual(
      values = c("TRUE" = "#D73027", "FALSE" = "grey50"),
      labels = c("TRUE" = "p < 0.05", "FALSE" = "NS"),
      name   = NULL
    ) +
    scale_size_continuous(range = c(2, 6), name = "-log10(p)") +
    labs(
      title    = glue("{grp_name}: {comp_label}"),
      subtitle = glue(
        "Outcome: {pos_class} vs others  |  adjusted OR (V1 taxa)"
      ),
      x = "Adjusted OR (95% CI)", y = NULL
    ) +
    theme_classic(base_size = 12) +
    theme(
      plot.title      = element_text(face = "bold"),
      axis.text.y     = element_text(face = "italic", size = 9),
      legend.position = "bottom"
    )

  ggsave(
    make_path(dir_pdf, glue("6_forest_{grp_name}_{comp_label}")),
    p_forest, width = 7, height = 6
  )
  write.csv(
    results |>
      select("taxa_clean", "OR", "LCI", "UCI", "p.value", "sig"),
    make_path(
      dir_tab, glue("6_logistic_{grp_name}_{comp_label}"), "csv"
    ),
    row.names = FALSE
  )
  invisible(results)
}

lr_comps <- list(
  list(grp = "Grp1", pos = "CIN2_CIN2", label = "treatment_failure"),
  list(grp = "Grp1", pos = "WNL_CIN2",  label = "recurrence"),
  list(grp = "Grp2", pos = "CIN2_CIN2", label = "persistent"),
  list(grp = "Grp2", pos = "WNL_CIN2",  label = "relapse")
)

ps_v1_traj3 <- subset_samples(
  ps2.traj,
  Timepoint == "V1" &
    studygrp2 %in% c("Grp1", "Grp2") &
    !is.na(trajectory)
)

invisible(lapply(lr_comps, function(lc) {
  grp_val <- lc$grp
  pos_val <- lc$pos
  lbl_val <- lc$label
  sd_tmp  <- data.frame(sample_data(ps_v1_traj3), stringsAsFactors = FALSE)
  keep    <- sd_tmp$studygrp2 == grp_val
  ps_sub  <- prune_samples(keep, ps_v1_traj3)
  run_logistic(
    ps_sub, "trajectory", pos_val,
    covars     = c("hivstatus", "simpleCST"),
    grp_name   = grp_val,
    comp_label = lbl_val
  )
}))
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
    !is.na(simpleCST)
  ) |>
  mutate(
    dysbiosis_score = cst_score_map[simpleCST],
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
#  Done
##############################################################################
cat("\n")
cat(strrep("=", 66), "\n")
cat("  CCM Phase 2 Trajectory Deep Analysis — COMPLETE\n")
cat(glue("  Output: {dir_pdf}/\n\n"))
cat("  Fig 1 : Trajectory distribution bar chart\n")
cat("  Fig 2 : CST transition Sankey (Grp1, Grp2)\n")
cat("  Fig 3 : Hierarchical clustering heatmap (V1 samples)\n")
cat("  Fig 4 : Spearman co-occurrence heatmap (Grp1, Grp2)\n")
cat("  Fig 5 : MRS — Microbial Risk Score (4 comparisons)\n")
cat("  Fig 6 : Logistic regression forest plot (4 comparisons)\n")
cat("  Fig 7 : JSD longitudinal change (V1->V2, V1->V3)\n")
cat("  Fig 8 : Dysbiosis score (CST-proxy) over time\n")
cat("  Fig 9 : Combined summary figure (patchwork)\n")
cat(strrep("=", 66), "\n")


getwd()
