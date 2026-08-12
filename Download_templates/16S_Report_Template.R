############################################################
#=================   16S Analysis + Report Template   =================#
############################################################
# Full pipeline: input -> filter -> analysis -> DA -> report.
# Follow the sections top to bottom. Every Go_*Info() field below is
# listed explicitly (not just the ones a given study needs) so you can
# see every option; leave NA/NULL for anything that doesn't apply.
# Run any Go_*Info() with no arguments in the console to print its guide,
# e.g. Go_expInfo() alone lists kit/prep/spikein numbers.
############################################################

if (!is.null(dev.list())) dev.off()
rm(list = ls())


###########################################
#=========    Library Loading    =========#
###########################################
# First-time setup only:
# devtools::install_github("bbagy/Gotools", force = TRUE)
# devtools::install_github("bbagy/ConDAdist", force = TRUE)
library(magick)
library(Gotools)
library(ConDAdist)
Gotool_dependency()
condadist_dependency(ask = FALSE)


###########################################
#=========    Project Setup      =========#
###########################################
project            <- "FILL_IN_PROJECT_NAME"
masterdir          <- "/PATH/TO/PROJECT_ROOT"
currentwd          <- file.path(masterdir, "2_Analysis", "FILL_IN_CURRENT_ANALYSIS_DIR")
summary_report_dir <- file.path(masterdir, "3_Summary_report")
report_logo_path   <- NULL # NULL = bundled Columbia/Uhlemann Lab logo, or set a custom logo file path

expinfo_kit_number     <- NA # DNeasy PowerSoil Pro Kit = 6, see Go_expInfo() guide
expinfo_pos_number     <- NA # ZymoBIOMICS Microbial Community Standard = 1
expinfo_prep_number    <- NA # Illumina V3V4 = 1
expinfo_spikein_number <- NA # No Spike-in Control = 2

setwd(currentwd)
if (!dir.exists(summary_report_dir)) dir.create(summary_report_dir, recursive = TRUE)


###########################################
#=========    Input and Import   =========#
###########################################
# If mapping is missing, use Go_emptyMap() -> fill in -> save to 3_map/
mapping_file   <- "3_map/MAPPING.csv"
raw_asv_file   <- "1_out/FILL_IN.psTotab.asvTable.csv"
track_file     <- "1_out/FILL_IN_track.csv"
tree_file      <- "1_out/TREE/tree.nwk"

sampledata <- read.csv(mapping_file, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)

final_asv_file <- find_latest_final_asv(project)
if (is.na(final_asv_file) || !file.exists(final_asv_file)) {
  # Go_blastASVs() defaults to the 16S ribosomal RNA DB bundled with Gotools.
  # Pass blastDB = "/path/to/other/blastDB" to use a different database.
  Go_blastASVs(project = project, asvsTable = raw_asv_file)
  final_asv_file <- find_latest_final_asv(project)
}
if (is.na(final_asv_file) || !file.exists(final_asv_file)) stop("No final_asvTable was found for Go_tabTops().")
validate_final_asv_input(final_asv_file, project)

ps.csv <- Go_tabTops(csv = final_asv_file, project = project)
tree   <- ape::read.tree(tree_file)
ps0    <- phyloseq::merge_phyloseq(ps.csv, phyloseq::sample_data(sampledata), phyloseq::phy_tree(tree))


###########################################
#=========  Filtering Workflow   =========#
###########################################
sample_type      <- "fecal" # fecal / vaginal / oral
cutoff           <- switch(sample_type, fecal = 7500, vaginal = 5000, oral = 3000, 5000)
seq_from         <- 397
seq_to           <- 450
abundance_cutoff <- 0.00005

ps1.prune <- phyloseq::prune_samples(phyloseq::sample_sums(ps0) > cutoff, ps0)

p.rare   <- Go_rare(ps1.prune, step = 1000, color = "TreatmentGroup", xlimit = 20000, label = NULL, se = FALSE) + ggplot2::theme(legend.position = "right")
rare_dir <- Go_path(project, pdf = "yes", table = "no", path = NULL)
ggplot2::ggsave(sprintf("%s/1_rarefaction.%s.%s.(Rarefaction).pdf", rare_dir$pdf, project, format(Sys.Date(), "%y%m%d")), plot = p.rare, device = "pdf", width = 8, height = 5, limitsize = FALSE)

Go_SeqLengths(ps1.prune)
ps1.size    <- Go_SeqLengths(ps1.prune, from = seq_from, to = seq_to)
ps_filtered <- Go_filter(ps1.size, cutoff = abundance_cutoff)


###########################################
#=========   Analysis Metadata   =========#
###########################################
meta0 <- data.frame(phyloseq::sample_data(ps_filtered), check.names = FALSE)

group_orders <- c("FILL_IN_Group1", "FILL_IN_Group2")
time_orders  <- c("FILL_IN_Baseline", "FILL_IN_Followup")
all_orders   <- c(group_orders, time_orders)

meta_main <- if ("AnalysisUse" %in% colnames(meta0)) subset(meta0, AnalysisUse == "Include") else meta0
meta_main$TreatmentGroup <- factor(meta_main$TreatmentGroup, levels = group_orders)

ps_main <- phyloseq::prune_samples(rownames(meta_main), ps_filtered)
phyloseq::sample_data(ps_main) <- phyloseq::sample_data(meta_main)

comparison_var    <- "TreatmentGroup"
comparison_orders <- group_orders
comparison_ps     <- ps_main

barchart_tax_ranks <- c("Phylum", "Class", "Order", "Family", "Genus", "Species")
distance_metrics   <- c("bray", "unifrac", "wunifrac")
alpha_metrics      <- c("Chao1", "Shannon")

barplot.col <- Go_myCols(custumCols = "cols3")
alpha.col   <- Go_myCols(piratepal = "basel")
bdiv.col    <- Go_myCols(piratepal = "basel")
da.col      <- c("#026CCBFF", "#F51E02FF")

# Repeated-measures / longitudinal: set to the subject ID column name (e.g. "StudyID")
# Cross-sectional: keep both NULL
longitudinal_id_var <- NULL
box_covariates       <- NULL # set to c("covariate") for ANCOVA; ignored when longitudinal_id_var is set
box_paired           <- longitudinal_id_var # drives Go_boxplot() model selection


############################################################
#=================     Core Analysis Run     =================#
############################################################
#-----------------------------------------#
#---------   Composition Plots   ---------#
#-----------------------------------------#
for (tax_rank in barchart_tax_ranks) {
  Go_barchart(psIN = ps_main, project = project, cutoff = 0.005, simple = FALSE, relative = TRUE,
              taxanames = tax_rank, cate.vars = comparison_var, facet = NULL, orders = comparison_orders,
              x_label = NULL, mycols = barplot.col, legend = "bottom",
              name = ensure_report_token(tax_rank), height = 4, width = 8)
}

Go_pheatmap(psIN = ps_main, project = project, title = NULL,
            group1 = comparison_var, group2 = NULL, group3 = NULL, group4 = NULL,
            Ntax = 50, name = ensure_report_token(comparison_var), col_orders = comparison_orders,
            show_rownames = TRUE, show_colnames = FALSE, cutree_rows = NA, cutree_cols = NA,
            cluster_rows = TRUE, cluster_cols = TRUE, showPhylum = TRUE, width = 10)

#-----------------------------------------#
#---------   Alpha Diversity     ---------#
#-----------------------------------------#
adiv_all <- Go_adiv(psIN = ps_main, project = project, alpha_metrics = alpha_metrics, name = ensure_report_token("AllSamples"))

Go_boxplot(df = adiv_all, project = project, cate.vars = comparison_var, outcomes = alpha_metrics,
           orders = comparison_orders, mycols = alpha.col, covariates = box_covariates, paired = box_paired,
           name = ensure_report_token(comparison_var), cutoff = 0.1, plotCols = 2, plotRows = 1)

#-----------------------------------------#
#---------   Beta Diversity      ---------#
#-----------------------------------------#
Go_bdivPM(psIN = comparison_ps, cate.vars = comparison_var, project = project, orders = comparison_orders,
          distance_metrics = distance_metrics, mycols = bdiv.col, name = ensure_report_token(comparison_var),
          facet = NULL, strata_var = longitudinal_id_var, height = 4.5, width = 7, plotCols = 3, plotRows = 1)


###########################################
#=========  Differential Abundance  =========#
###########################################
if (is.null(longitudinal_id_var)) {
  # Go_ConDaDist() returns a project/date-based output directory (same across
  # methods run the same day) - run all methods first, then call
  # Go_volcanoPlot() once; it scans that directory for every method's CSV.
  da_result_dir <- NULL
  for (method in c("ancombc2", "aldex2", "deseq2")) {
    da_result_dir <- Go_ConDaDist(psIN = comparison_ps, project = project, group_var = comparison_var,
                                   group_1 = comparison_orders[1], group_2 = comparison_orders[2],
                                   covariates = NULL, methods = method, distances = NULL, name = NULL)
  }
  Go_volcanoPlot(project = project, result = da_result_dir, fc = 0, mycols = da.col, name = NULL, font = 3, height = 5, width = 6)
} else {
  da_result_dir <- Go_ConDaDist(psIN = comparison_ps, project = project, group_var = comparison_var,
                                 group_1 = comparison_orders[1], group_2 = comparison_orders[2],
                                 covariates = NULL, methods = "ancombc2", random_effects = c(longitudinal_id_var),
                                 distances = NULL, name = NULL)
  Go_volcanoPlot(project = project, result = da_result_dir, fc = 0, mycols = da.col, name = NULL, font = 3, height = 5, width = 6)
}


###########################################
#=========  Report Objects       =========#
###########################################
dir_root       <- sprintf("%s_%s", project, format(Sys.Date(), "%y%m%d"))
alpha_tab_file <- sprintf("%s/table/adiv/adiv.%s.%s.%s.csv", dir_root, project, ensure_report_token("AllSamples"), format(Sys.Date(), "%y%m%d"))
da_plot_dir    <- file.path(dir_root, "pdf", "DA_plot")
da_labels      <- extract_da_plot_labels(da_plot_dir)

expInfo <- Go_expInfo(
  Project_name        = "FILL_IN human-readable study title",
  Samples_info        = sprintf("%s 16S samples", phyloseq::nsamples(ps_main)),
  Sequencing_date     = format(Sys.Date(), "%m/%d/%Y"),
  Sequencing_platform = "FILL_IN e.g. Illumina MiSeq V3 600 cycles",
  RNAseq_reference    = NA, # not used for 16S
  kit_number          = expinfo_kit_number,
  pos_number          = expinfo_pos_number,
  prep_number         = expinfo_prep_number,
  spikein_number      = expinfo_spikein_number,
  authorName1         = "FILL_IN_Name",
  authorEmail1        = "FILL_IN_email@cumc.columbia.edu",
  authorName2         = NULL, # optional second author
  authorEmail2        = NULL
)

dirInfo <- Go_dirInfo(
  Current_working_dir = currentwd,
  Image_locations      = sprintf("%s/pdf/", dir_root),
  DA_image_location    = sprintf("%s/pdf/DA_plot/", dir_root)
)

tabInfo <- Go_tabInfo(
  Taxa_Tab          = NA, # 16S uses ASVs_Tab instead
  ASVs_Tab          = final_asv_file,
  GTDBTK_Tab        = NA, # MAGs workflow only
  Tract_Tab         = track_file,
  Func_Tab          = NA, # HUMAnN3 functional tables (MG only)
  Other_Tab         = NA,
  Alpha_divTab      = alpha_tab_file,
  Alpha_div_LmerTab = NULL,
  RNAseq            = NA, # not used for 16S
  HumannTab         = NA,
  Tab1              = NA,
  Tab2              = NA,
  PermanovaTab      = NULL
)

imgInfo <- Go_imgInfo(
  Overview       = NA,
  Rarefaction    = "Rarefaction",
  Barchart       = barchart_tax_ranks,
  Bac.heatmap    = comparison_var,
  RNAseq.heatmap = NA, # RNAseq only
  HumannHeatmap  = NA, # MG/HUMAnN3 only
  Adivplot       = c(comparison_var),
  Foreplot       = NA,
  Bdivplot       = c(comparison_var),
  DAplot         = da_labels,
  EBplot         = NA, # MG functional extended barplot tokens
  Network        = NA,
  SHeatmap       = NA, # MG functional stacked heatmap tokens
  MAGs_tree      = NA, # MAGs workflow only
  MAGs_upset     = NA,
  MAGs_heatmap   = NA,
  Extra          = NULL
)

# Fill every field from the actual generated outputs once the run above is complete.
sumInfo <- Go_sumInfo(
  Prep_overview = "FILL_IN - extraction/sequencing QC summary, note any failed samples and final N per group",
  Taxa_overview = "FILL_IN - dominant taxa observed in the barchart",
  Alpha_div     = "FILL_IN - Chao1/Shannon comparison across groups, with p-values",
  Beta_div      = "FILL_IN - PERMANOVA results per distance metric",
  Bacterial_div = "",
  DA_test       = "FILL_IN - taxa flagged by ancombc2/aldex2/deseq2, direction of enrichment",
  Summary       = "",
  Conclusion    = "",
  Useful_info   = ""
)


###########################################
#=========  Render and Verify    =========#
###########################################
Go_renderReport(family = "16S", project = project,
                 expInfo = expInfo, dirInfo = dirInfo, tabInfo = tabInfo, imgInfo = imgInfo, sumInfo = sumInfo,
                 output_dir = summary_report_dir,
                 logo_path = report_logo_path)
