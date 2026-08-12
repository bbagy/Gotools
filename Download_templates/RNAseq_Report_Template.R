############################################################
#=================   RNAseq Analysis + Report Template   =================#
############################################################
# Full pipeline: input -> analysis -> DA -> report.
# Follow the sections top to bottom. Every Go_*Info() field below is
# listed explicitly (not just the ones a given study needs) so you can
# see every option; leave NA/NULL for anything that doesn't apply.
# Run any Go_*Info() with no arguments in the console to print its guide,
# e.g. Go_expInfo() alone lists kit/prep/spikein numbers.
#
# NOTE: no 2026 real analysis has run this report path yet on real data;
# this template is structurally consistent with the 16S/MG templates but
# has not been checked against a completed real RNAseq report.
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
currentwd          <- "/PATH/TO/RNASEQ_ANALYSIS_DIR"
summary_report_dir <- "/PATH/TO/PROJECT_ROOT/3_Summary_report"
report_logo_path   <- NULL # NULL = bundled Columbia/Uhlemann Lab logo, or set a custom logo file path

expinfo_kit_number     <- 1 # ZymoBIOMICS 96 MagBead DNA/RNA Kit
expinfo_pos_number     <- 1 # ZymoBIOMICS Microbial Community Standard
expinfo_prep_number    <- 5 # Zymo Quick-ITS Plus NGS Library Prep Kit
expinfo_spikein_number <- 2 # No Spike-in Control

setwd(currentwd)
if (!dir.exists(summary_report_dir)) dir.create(summary_report_dir, recursive = TRUE)


###########################################
#=========    Input and Import   =========#
###########################################
mapping_file <- "3_map/MAPPING.csv"
rnaseq_table <- "1_out/FILL_IN_RNAseq_expression_table.csv"
track_file   <- "1_out/FILL_IN_QC_summary.csv"

sampledata <- read.csv(mapping_file, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
ps.rna <- Go_tabTops(csv = rnaseq_table, project = project)
ps.rna <- phyloseq::merge_phyloseq(ps.rna, phyloseq::sample_data(data.frame(sampledata)))

comparison_var    <- "TreatmentGroup"
comparison_orders <- c("FILL_IN_Group1", "FILL_IN_Group2")
da.col            <- c("#74ADD1", "#E31A1C")


############################################################
#=================     Core Analysis Run     =================#
############################################################
#-----------------------------------------#
#---------   Expression Heatmap  ---------#
#-----------------------------------------#
Go_pheatmap(psIN = ps.rna, project = project, Ntax = 50, title = NULL, group1 = comparison_var,
            show_rownames = TRUE, show_colnames = FALSE, showPhylum = FALSE,
            cluster_rows = TRUE, cluster_cols = TRUE, cutree_cols = NA, cutree_rows = NA, width = 10)


###########################################
#=========  Differential Abundance  =========#
###########################################
conda_results <- list()
for (method in c("ancombc2", "aldex2", "deseq2")) {
  conda_results[[method]] <- Go_ConDaDist(psIN = ps.rna, project = project, group_var = comparison_var,
                                           group_1 = comparison_orders[1], group_2 = comparison_orders[2],
                                           covariates = NULL, methods = method, distances = NULL, name = NULL)
}
for (method in names(conda_results)) {
  Go_volcanoPlot(project = project, result = conda_results[[method]], mycols = da.col, fc = 0, name = NULL, font = 3, height = 5, width = 6)
}


###########################################
#=========  Report Objects       =========#
###########################################
dir_root    <- sprintf("%s_%s", project, format(Sys.Date(), "%y%m%d"))
da_plot_dir <- file.path(dir_root, "pdf", "DA_plot")
da_labels   <- extract_da_plot_labels(da_plot_dir)

expInfo <- Go_expInfo(
  Project_name        = "FILL_IN human-readable study title",
  Samples_info        = sprintf("%s RNAseq samples", phyloseq::nsamples(ps.rna)),
  Sequencing_date     = format(Sys.Date(), "%m/%d/%Y"),
  Sequencing_platform = "FILL_IN e.g. Illumina NovaSeq",
  RNAseq_reference    = "FILL_IN reference genome/transcriptome build",
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
  Taxa_Tab          = NA, # 16S/MG only
  ASVs_Tab          = NA, # 16S only
  GTDBTK_Tab        = NA, # MAGs workflow only
  Tract_Tab         = track_file,
  Func_Tab          = NA, # MG/HUMAnN3 only
  Other_Tab         = NA,
  Alpha_divTab      = NA,
  Alpha_div_LmerTab = NULL,
  RNAseq            = rnaseq_table,
  HumannTab         = NA,
  Tab1              = NA,
  Tab2              = NA,
  PermanovaTab      = NULL
)

imgInfo <- Go_imgInfo(
  Overview       = NA,
  Rarefaction    = NA,
  Barchart       = NA, # 16S/MG only
  Bac.heatmap    = NA,
  RNAseq.heatmap = c("taxa heatmap"),
  HumannHeatmap  = NA,
  Adivplot       = NA,
  Foreplot       = NA,
  Bdivplot       = NA,
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
  Prep_overview = "FILL_IN - library prep/sequencing QC summary, note any failed samples",
  Taxa_overview = "",
  Alpha_div     = "",
  Beta_div      = "",
  Bacterial_div = "",
  DA_test       = "FILL_IN - genes/transcripts flagged by ancombc2/aldex2/deseq2, direction of change",
  Summary       = "",
  Conclusion    = "",
  Useful_info   = ""
)


###########################################
#=========  Render and Verify    =========#
###########################################
Go_renderReport(family = "RNAseq", project = project,
                 expInfo = expInfo, dirInfo = dirInfo, tabInfo = tabInfo, imgInfo = imgInfo, sumInfo = sumInfo,
                 output_dir = summary_report_dir,
                 logo_path = report_logo_path)
