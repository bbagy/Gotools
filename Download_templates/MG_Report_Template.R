############################################################
#=================   MG Analysis + Report Template   =================#
############################################################
# Full pipeline: input -> filter -> analysis -> DA -> HUMAnN3 -> report.
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
currentwd          <- "/PATH/TO/MG_ANALYSIS_DIR"
summary_report_dir <- "/PATH/TO/PROJECT_ROOT/3_Summary_report"
report_logo_path   <- NULL # NULL = bundled Columbia/Uhlemann Lab logo, or set a custom logo file path

expinfo_kit_number     <- 2 # DNeasy PowerSoil Pro Kit
expinfo_pos_number     <- 1 # ZymoBIOMICS Microbial Community Standard
expinfo_prep_number    <- 6 # Illumina DNA Prep
expinfo_spikein_number <- 2 # No Spike-in Control

setwd(currentwd)
if (!dir.exists(summary_report_dir)) dir.create(summary_report_dir, recursive = TRUE)


###########################################
#=========    Input and Import   =========#
###########################################
mapping_file <- "3_map/MAPPING.csv"
kraken2_file <- "1_out/FILL_IN_KRAKEN2_MPA.txt"
bracken_file <- "1_out/FILL_IN_BRACKEN_MPA.txt"
track_file   <- "1_out/FILL_IN_QC_summary.csv"
kingdom      <- "Bacteria"

sampledata <- read.csv(mapping_file, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)

# Go_kraken2Tops() builds the otu+taxonomy phyloseq directly. Do not rebuild from exported csv files.
ps0     <- Go_kraken2Tops(kraken2 = kraken2_file, bracken = bracken_file, kingdom = kingdom)
ps0     <- phyloseq::merge_phyloseq(ps0, phyloseq::sample_data(data.frame(sampledata)))
ps_main <- Go_filter(ps0, cutoff = 0.00001)


###########################################
#=========   Analysis Metadata   =========#
###########################################
comparison_var    <- "TreatmentGroup"
comparison_orders <- c("FILL_IN_Group1", "FILL_IN_Group2")
comparison_ps     <- ps_main

barchart_tax_ranks <- c("phylum", "class", "order", "family", "genus", "species")
distance_metrics   <- c("jsd", "jaccard")
alpha_metrics      <- c("Chao1", "Shannon")

barplot.col <- Go_myCols(custumCols = "cols3")
alpha.col   <- Go_myCols(piratepal = "basel")
bdiv.col    <- Go_myCols(piratepal = "basel")
da.col      <- c("#74ADD1", "#E31A1C")

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
              taxanames = tax_rank, cate.vars = comparison_var, facet = NULL, legend = "bottom",
              orders = comparison_orders, x_label = NULL, mycols = barplot.col,
              name = sprintf("(%s)", tools::toTitleCase(tax_rank)), height = 4, width = 8)
}

Go_pheatmap(psIN = ps_main, project = project, Ntax = 50, group1 = comparison_var,
            col_orders = comparison_orders, name = ensure_report_token(comparison_var), width = 8)

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
  Go_volcanoPlot(project = project, result = da_result_dir, mycols = da.col, fc = 0, name = NULL, font = 3, height = 5, width = 6)
} else {
  da_result_dir <- Go_ConDaDist(psIN = comparison_ps, project = project, group_var = comparison_var,
                                 group_1 = comparison_orders[1], group_2 = comparison_orders[2],
                                 covariates = NULL, methods = "ancombc2", random_effects = c(longitudinal_id_var),
                                 distances = NULL, name = NULL)
  Go_volcanoPlot(project = project, result = da_result_dir, mycols = da.col, fc = 0, name = NULL, font = 3, height = 5, width = 6)
}


###########################################
#=========   HUMAnN3 Functional  =========#
###########################################
# If HUMAnN3 tables are paired technical columns (*_R1/*_R2), set collapse_pairs = TRUE.
# Go_function2ps() auto-detects _R1/_R2, _r1/_r2, _1/_2 patterns and collapses by sum.
humann_pathabundance_file <- "1_out/FILL_IN_PATHABUNDANCE.tsv"
humann_genefamilies_file  <- "1_out/FILL_IN_GENEFAMILIES.tsv"
humann_collapse_pairs <- FALSE
humann_has_path <- file.exists(humann_pathabundance_file)
humann_has_kegg <- file.exists(humann_genefamilies_file)
humann_enabled  <- humann_has_path || humann_has_kegg

functional_mvar     <- comparison_var
functional_group1   <- comparison_orders[1]
functional_group2   <- comparison_orders[2]
functional_p        <- 0.1
functional_eb_tokens <- character(0)
functional_sh_tokens <- character(0)
ps2.path <- NULL
ps2.kegg <- NULL

if (humann_has_path) {
  ps2.path <- Go_function2ps(tabPath = humann_pathabundance_file, project = project, func.type = "Humann", name = NULL, collapse_pairs = humann_collapse_pairs)
  ps2.path <- phyloseq::merge_phyloseq(ps2.path, phyloseq::sample_data(data.frame(sampledata)))

  Go_extendedBarplot(psIN = ps2.path, project = project, mvar = functional_mvar, group1 = functional_group1, group2 = functional_group2, func = "path.des", wilcox.p = functional_p, name = "PATHWAY", height = 3, width = 10)

  Go_psTotab(ps2.path, project = "PATH")
  path_df <- read.csv(sprintf("1_out/PATH.%s.psTotab.asvTable.csv", format(Sys.Date(), "%y%m%d")), row.names = 1, check.names = FALSE)
  path_df$pathway <- NULL
  path_df$path.des <- NULL

  Go_Groupheatmap(df = path_df, SampleData = sampledata, project = project, Group = functional_mvar, orders = comparison_orders, title = NULL, mycols = NULL, name = ensure_report_token("PATHWAY"), top_n = 50, normalization = "log", width = 11, height = 8, x_label = NULL)

  functional_eb_tokens <- c(functional_eb_tokens, sprintf("%s.vs.%s.PATHWAY", functional_group1, functional_group2))
  functional_sh_tokens <- c(functional_sh_tokens, "PATHWAY")
}

if (humann_has_kegg) {
  ps2.kegg <- Go_function2ps(tabPath = humann_genefamilies_file, project = project, func.type = "Humann", name = NULL, collapse_pairs = humann_collapse_pairs)
  ps2.kegg <- phyloseq::merge_phyloseq(ps2.kegg, phyloseq::sample_data(data.frame(sampledata)))

  Go_extendedBarplot(psIN = ps2.kegg, project = project, mvar = functional_mvar, group1 = functional_group1, group2 = functional_group2, func = "KO.des", wilcox.p = functional_p, name = "KEGG", height = 3, width = 10)

  Go_psTotab(ps2.kegg, project = "KEGG")
  kegg_df <- read.csv(sprintf("1_out/KEGG.%s.psTotab.asvTable.csv", format(Sys.Date(), "%y%m%d")), row.names = 1, check.names = FALSE)
  kegg_df$KO <- NULL
  kegg_df$KO.des <- NULL

  Go_Groupheatmap(df = kegg_df, SampleData = sampledata, project = project, Group = functional_mvar, orders = comparison_orders, title = NULL, mycols = NULL, name = ensure_report_token("KEGG"), top_n = 50, normalization = "log", width = 11, height = 8, x_label = NULL)

  functional_eb_tokens <- c(functional_eb_tokens, sprintf("%s.vs.%s.KEGG", functional_group1, functional_group2))
  functional_sh_tokens <- c(functional_sh_tokens, "KEGG")
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
  Samples_info        = sprintf("%s metagenomic samples", phyloseq::nsamples(ps_main)),
  Sequencing_date     = format(Sys.Date(), "%m/%d/%Y"),
  Sequencing_platform = "FILL_IN e.g. Illumina NovaSeq / Element AVITI 2x150",
  RNAseq_reference    = NA, # not used for MG
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
  Taxa_Tab          = bracken_file,
  ASVs_Tab          = NA, # 16S only
  GTDBTK_Tab        = NA, # MAGs workflow only
  Tract_Tab         = track_file,
  Func_Tab          = if (humann_enabled) c(humann_genefamilies_file, humann_pathabundance_file) else NA,
  Other_Tab         = NA,
  Alpha_divTab      = alpha_tab_file,
  Alpha_div_LmerTab = NULL,
  RNAseq            = NA, # RNAseq only
  HumannTab         = NA,
  Tab1              = NA,
  Tab2              = NA,
  PermanovaTab      = NULL
)

imgInfo <- Go_imgInfo(
  Overview       = NA,
  Rarefaction    = NA,
  Barchart       = c("Phylum", "Class", "Order", "Family", "Genus", "Species"),
  Bac.heatmap    = comparison_var,
  RNAseq.heatmap = NA, # RNAseq only
  HumannHeatmap  = NA,
  Adivplot       = c(comparison_var),
  Foreplot       = NA,
  Bdivplot       = c(comparison_var),
  DAplot         = da_labels,
  EBplot         = if (length(functional_eb_tokens) > 0) functional_eb_tokens else NA,
  Network        = NA,
  SHeatmap       = if (length(functional_sh_tokens) > 0) functional_sh_tokens else NA,
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
Go_renderReport(family = "MG", project = project,
                 expInfo = expInfo, dirInfo = dirInfo, tabInfo = tabInfo, imgInfo = imgInfo, sumInfo = sumInfo,
                 output_dir = summary_report_dir,
                 logo_path = report_logo_path)
