############################################################
#=================   MG Quickstart (bundled demo data)   =================#
############################################################
# Runs the full MG pipeline end to end against the synthetic demo dataset
# bundled with Gotools, so you can see a real report before wiring up your
# own data. All taxa, counts, and sample groups here are synthetic.
#
# For your own analysis, copy Download_templates/MG_Report_Template.R
# instead and fill in your own paths.
############################################################

if (!is.null(dev.list())) dev.off()
rm(list = ls())

library(magick)
library(Gotools)
library(ConDAdist)
Gotool_dependency()
condadist_dependency(ask = FALSE)

project    <- "TutorialMG"
currentwd  <- file.path(path.expand("~"), "Gotools_Quickstart_MG")
summary_report_dir <- file.path(currentwd, "3_Summary_report")

if (!dir.exists(currentwd)) {
  dir.create(currentwd, recursive = TRUE)
  file.copy(system.file("tutorial_data", "MG", package = "Gotools"), currentwd, recursive = TRUE)
  file.rename(file.path(currentwd, "MG", "1_out"), file.path(currentwd, "1_out"))
  file.rename(file.path(currentwd, "MG", "3_map"), file.path(currentwd, "3_map"))
  unlink(file.path(currentwd, "MG"), recursive = TRUE)
}
setwd(currentwd)
if (!dir.exists(summary_report_dir)) dir.create(summary_report_dir, recursive = TRUE)


###########################################
#=========    Input and Import   =========#
###########################################
sampledata <- read.csv("3_map/MAPPING.csv", row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)

ps0 <- phyloseq::merge_phyloseq(
  Go_kraken2Tops(kraken2 = "1_out/KRAKEN2_MPA.txt", bracken = "1_out/BRACKEN_MPA.txt", kingdom = "Bacteria"),
  phyloseq::sample_data(sampledata)
)
ps_main <- Go_filter(ps0, cutoff = 0.00001)


###########################################
#=========   Analysis Metadata   =========#
###########################################
comparison_var    <- "TreatmentGroup"
comparison_orders <- c("Control", "Treatment")

barchart_tax_ranks <- c("phylum", "genus")
distance_metrics   <- c("jaccard")
alpha_metrics      <- c("Chao1", "Shannon")

barplot.col <- Go_myCols(custumCols = "cols3")
alpha.col   <- Go_myCols(piratepal = "basel")
bdiv.col    <- Go_myCols(piratepal = "basel")
da.col      <- c("#74ADD1", "#E31A1C")


############################################################
#=================     Core Analysis Run     =================#
############################################################
for (tax_rank in barchart_tax_ranks) {
  Go_barchart(psIN = ps_main, project = project, cutoff = 0.005, simple = FALSE, relative = TRUE,
              taxanames = tax_rank, cate.vars = comparison_var, facet = NULL, legend = "bottom",
              orders = comparison_orders, x_label = NULL, mycols = barplot.col,
              name = sprintf("(%s)", tools::toTitleCase(tax_rank)), height = 4, width = 8)
}

Go_pheatmap(psIN = ps_main, project = project, Ntax = 20, group1 = comparison_var,
            col_orders = comparison_orders, name = ensure_report_token(comparison_var), width = 8)

adiv_all <- Go_adiv(psIN = ps_main, project = project, alpha_metrics = alpha_metrics, name = ensure_report_token("AllSamples"))
Go_boxplot(df = adiv_all, project = project, cate.vars = comparison_var, outcomes = alpha_metrics,
           orders = comparison_orders, mycols = alpha.col, covariates = NULL, paired = NULL,
           name = ensure_report_token(comparison_var), cutoff = 0.1, plotCols = 2, plotRows = 1)

Go_bdivPM(psIN = ps_main, cate.vars = comparison_var, project = project, orders = comparison_orders,
          distance_metrics = distance_metrics, mycols = bdiv.col, name = ensure_report_token(comparison_var),
          facet = NULL, strata_var = NULL, height = 4.5, width = 7, plotCols = 1, plotRows = 1)

da_result_dir <- NULL
for (method in c("ancombc2", "aldex2", "deseq2")) {
  da_result_dir <- Go_ConDaDist(psIN = ps_main, project = project, group_var = comparison_var,
                                 group_1 = comparison_orders[1], group_2 = comparison_orders[2],
                                 covariates = NULL, methods = method, distances = NULL, name = NULL)
}
Go_volcanoPlot(project = project, result = da_result_dir, mycols = da.col, fc = 0, name = NULL, font = 3, height = 5, width = 6)


###########################################
#=========  Report Objects       =========#
###########################################
dir_root       <- sprintf("%s_%s", project, format(Sys.Date(), "%y%m%d"))
alpha_tab_file <- sprintf("%s/table/adiv/adiv.%s.%s.%s.csv", dir_root, project, ensure_report_token("AllSamples"), format(Sys.Date(), "%y%m%d"))
da_labels      <- extract_da_plot_labels(file.path(dir_root, "pdf", "DA_plot"))

expInfo <- Go_expInfo(
  Project_name = "Gotools MG Quickstart", Samples_info = sprintf("%s samples (synthetic groups)", phyloseq::nsamples(ps_main)),
  Sequencing_date = format(Sys.Date(), "%m/%d/%Y"), Sequencing_platform = "Demo data",
  kit_number = 2, pos_number = 1, prep_number = 6, spikein_number = 2,
  authorName1 = "Gotools Quickstart", authorEmail1 = NA
)
dirInfo <- Go_dirInfo(Current_working_dir = currentwd, Image_locations = sprintf("%s/pdf/", dir_root), DA_image_location = sprintf("%s/pdf/DA_plot/", dir_root))
tabInfo <- Go_tabInfo(Taxa_Tab = "1_out/BRACKEN_MPA.txt", Tract_Tab = "1_out/QC_summary.csv", Alpha_divTab = alpha_tab_file)
imgInfo <- Go_imgInfo(Barchart = barchart_tax_ranks, Bac.heatmap = comparison_var,
                       Adivplot = comparison_var, Bdivplot = comparison_var, DAplot = da_labels)
sumInfo <- Go_sumInfo(Prep_overview = "Demo dataset: fully synthetic taxa, counts, and groups.")

out <- Go_renderReport(family = "MG", project = project, expInfo = expInfo, dirInfo = dirInfo,
                        tabInfo = tabInfo, imgInfo = imgInfo, sumInfo = sumInfo, output_dir = summary_report_dir)

message("Quickstart report rendered: ", out)
