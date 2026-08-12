############################################################
#=================   16S Quickstart (bundled demo data)   =================#
############################################################
# Runs the full 16S pipeline end to end against the tutorial dataset bundled
# with Gotools, so you can see a real report before wiring up your own data.
# The ASV sequences are real public NCBI 16S reference sequences (not
# patient data); counts and sample groups are synthetic.
#
# For your own analysis, copy Download_templates/16S_Report_Template.R
# instead and fill in your own paths.
############################################################

if (!is.null(dev.list())) dev.off()
rm(list = ls())

library(magick)
library(Gotools)
library(ConDAdist)
Gotool_dependency()
condadist_dependency(ask = FALSE)

project    <- "Tutorial16S"
currentwd  <- file.path(path.expand("~"), "Gotools_Quickstart_16S")
summary_report_dir <- file.path(currentwd, "3_Summary_report")

if (!dir.exists(currentwd)) {
  dir.create(currentwd, recursive = TRUE)
  file.copy(system.file("tutorial_data", "16S", package = "Gotools"), currentwd, recursive = TRUE)
  file.rename(file.path(currentwd, "16S", "1_out"), file.path(currentwd, "1_out"))
  file.rename(file.path(currentwd, "16S", "3_map"), file.path(currentwd, "3_map"))
  unlink(file.path(currentwd, "16S"), recursive = TRUE)
}
setwd(currentwd)
if (!dir.exists(summary_report_dir)) dir.create(summary_report_dir, recursive = TRUE)


###########################################
#=========    Input and Import   =========#
###########################################
raw_asv_file <- "1_out/Tutorial16S.psTotab.asvTable.csv"
track_file   <- "1_out/track.csv"
tree_file    <- "1_out/TREE/tree.nwk"

sampledata <- read.csv("3_map/MAPPING.csv", row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)

final_asv_file <- find_latest_final_asv(project)
if (is.na(final_asv_file) || !file.exists(final_asv_file)) {
  Go_blastASVs(project = project, asvsTable = raw_asv_file) # uses the bundled 16S BLAST DB
  final_asv_file <- find_latest_final_asv(project)
}
validate_final_asv_input(final_asv_file, project)

ps.csv <- Go_tabTops(csv = final_asv_file, project = project)
tree   <- ape::read.tree(tree_file)
ps_main <- phyloseq::merge_phyloseq(ps.csv, phyloseq::sample_data(sampledata), phyloseq::phy_tree(tree))


###########################################
#=========   Analysis Metadata   =========#
###########################################
comparison_var    <- "TreatmentGroup"
comparison_orders <- c("Control", "Treatment")
comparison_ps     <- ps_main

barchart_tax_ranks <- c("Phylum", "Genus")
distance_metrics   <- c("bray")
alpha_metrics      <- c("Chao1", "Shannon")

barplot.col <- Go_myCols(custumCols = "cols3")
alpha.col   <- Go_myCols(piratepal = "basel")
bdiv.col    <- Go_myCols(piratepal = "basel")
da.col      <- c("#026CCBFF", "#F51E02FF")


############################################################
#=================     Core Analysis Run     =================#
############################################################
for (tax_rank in barchart_tax_ranks) {
  Go_barchart(psIN = ps_main, project = project, cutoff = 0.005, simple = FALSE, relative = TRUE,
              taxanames = tax_rank, cate.vars = comparison_var, facet = NULL, orders = comparison_orders,
              x_label = NULL, mycols = barplot.col, legend = "bottom",
              name = ensure_report_token(tax_rank), height = 4, width = 8)
}

Go_pheatmap(psIN = ps_main, project = project, title = NULL,
            group1 = comparison_var, group2 = NULL, group3 = NULL, group4 = NULL,
            Ntax = 20, name = ensure_report_token(comparison_var), col_orders = comparison_orders,
            show_rownames = TRUE, show_colnames = FALSE, cutree_rows = NA, cutree_cols = NA,
            cluster_rows = TRUE, cluster_cols = TRUE, showPhylum = TRUE, width = 10)

adiv_all <- Go_adiv(psIN = ps_main, project = project, alpha_metrics = alpha_metrics, name = ensure_report_token("AllSamples"))
Go_boxplot(df = adiv_all, project = project, cate.vars = comparison_var, outcomes = alpha_metrics,
           orders = comparison_orders, mycols = alpha.col, covariates = NULL, paired = NULL,
           name = ensure_report_token(comparison_var), cutoff = 0.1, plotCols = 2, plotRows = 1)

Go_bdivPM(psIN = comparison_ps, cate.vars = comparison_var, project = project, orders = comparison_orders,
          distance_metrics = distance_metrics, mycols = bdiv.col, name = ensure_report_token(comparison_var),
          facet = NULL, strata_var = NULL, height = 4.5, width = 7, plotCols = 1, plotRows = 1)

da_result_dir <- NULL
for (method in c("ancombc2", "aldex2", "deseq2")) {
  da_result_dir <- Go_ConDaDist(psIN = comparison_ps, project = project, group_var = comparison_var,
                                 group_1 = comparison_orders[1], group_2 = comparison_orders[2],
                                 covariates = NULL, methods = method, distances = NULL, name = NULL)
}
Go_volcanoPlot(project = project, result = da_result_dir, fc = 0, mycols = da.col, name = NULL, font = 3, height = 5, width = 6)


###########################################
#=========  Report Objects       =========#
###########################################
dir_root       <- sprintf("%s_%s", project, format(Sys.Date(), "%y%m%d"))
alpha_tab_file <- sprintf("%s/table/adiv/adiv.%s.%s.%s.csv", dir_root, project, ensure_report_token("AllSamples"), format(Sys.Date(), "%y%m%d"))
da_labels      <- extract_da_plot_labels(file.path(dir_root, "pdf", "DA_plot"))

expInfo <- Go_expInfo(
  Project_name = "Gotools 16S Quickstart", Samples_info = sprintf("%s samples (synthetic groups)", phyloseq::nsamples(ps_main)),
  Sequencing_date = format(Sys.Date(), "%m/%d/%Y"), Sequencing_platform = "Demo data",
  kit_number = 6, pos_number = 1, prep_number = 1, spikein_number = 2,
  authorName1 = "Gotools Quickstart", authorEmail1 = NA
)
dirInfo <- Go_dirInfo(Current_working_dir = currentwd, Image_locations = sprintf("%s/pdf/", dir_root), DA_image_location = sprintf("%s/pdf/DA_plot/", dir_root))
tabInfo <- Go_tabInfo(ASVs_Tab = final_asv_file, Tract_Tab = track_file, Alpha_divTab = alpha_tab_file)
imgInfo <- Go_imgInfo(Rarefaction = NA, Barchart = barchart_tax_ranks, Bac.heatmap = comparison_var,
                       Adivplot = comparison_var, Bdivplot = comparison_var, DAplot = da_labels)
sumInfo <- Go_sumInfo(Prep_overview = "Demo dataset: real public 16S reference sequences with synthetic counts and groups.")

out <- Go_renderReport(family = "16S", project = project, expInfo = expInfo, dirInfo = dirInfo,
                        tabInfo = tabInfo, imgInfo = imgInfo, sumInfo = sumInfo, output_dir = summary_report_dir)

message("Quickstart report rendered: ", out)
