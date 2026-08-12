############################################################
#=================   RNAseq Quickstart (bundled demo data)   =================#
############################################################
# Runs the full RNAseq pipeline end to end against the synthetic demo dataset
# bundled with Gotools, so you can see a real report before wiring up your
# own data. All genes, counts, and sample groups here are synthetic.
#
# For your own analysis, copy Download_templates/RNAseq_Report_Template.R
# instead and fill in your own paths.
############################################################

if (!is.null(dev.list())) dev.off()
rm(list = ls())

library(magick)
library(Gotools)
library(ConDAdist)
Gotool_dependency()
condadist_dependency(ask = FALSE)

project    <- "TutorialRNAseq"
currentwd  <- file.path(path.expand("~"), "Gotools_Quickstart_RNAseq")
summary_report_dir <- file.path(currentwd, "3_Summary_report")

if (!dir.exists(currentwd)) {
  dir.create(currentwd, recursive = TRUE)
  file.copy(system.file("tutorial_data", "RNAseq", package = "Gotools"), currentwd, recursive = TRUE)
  file.rename(file.path(currentwd, "RNAseq", "1_out"), file.path(currentwd, "1_out"))
  file.rename(file.path(currentwd, "RNAseq", "3_map"), file.path(currentwd, "3_map"))
  unlink(file.path(currentwd, "RNAseq"), recursive = TRUE)
}
setwd(currentwd)
if (!dir.exists(summary_report_dir)) dir.create(summary_report_dir, recursive = TRUE)


###########################################
#=========    Input and Import   =========#
###########################################
rnaseq_table <- "1_out/RNAseq_expression_table.csv"
sampledata <- read.csv("3_map/MAPPING.csv", row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
ps.rna <- phyloseq::merge_phyloseq(Go_tabTops(csv = rnaseq_table, project = project), phyloseq::sample_data(sampledata))

comparison_var    <- "TreatmentGroup"
comparison_orders <- c("Control", "Treatment")
da.col            <- c("#74ADD1", "#E31A1C")


############################################################
#=================     Core Analysis Run     =================#
############################################################
Go_pheatmap(psIN = ps.rna, project = project, Ntax = 20, title = NULL, group1 = comparison_var,
            show_rownames = TRUE, show_colnames = FALSE, showPhylum = FALSE,
            cluster_rows = TRUE, cluster_cols = TRUE, cutree_cols = NA, cutree_rows = NA, width = 8)

da_result_dir <- NULL
for (method in c("ancombc2", "aldex2", "deseq2")) {
  da_result_dir <- Go_ConDaDist(psIN = ps.rna, project = project, group_var = comparison_var,
                                 group_1 = comparison_orders[1], group_2 = comparison_orders[2],
                                 covariates = NULL, methods = method, distances = NULL, name = NULL)
}
Go_volcanoPlot(project = project, result = da_result_dir, mycols = da.col, fc = 0, name = NULL, font = 3, height = 5, width = 6)


###########################################
#=========  Report Objects       =========#
###########################################
dir_root    <- sprintf("%s_%s", project, format(Sys.Date(), "%y%m%d"))
da_labels   <- extract_da_plot_labels(file.path(dir_root, "pdf", "DA_plot"))

expInfo <- Go_expInfo(
  Project_name = "Gotools RNAseq Quickstart", Samples_info = sprintf("%s samples (synthetic groups)", phyloseq::nsamples(ps.rna)),
  Sequencing_date = format(Sys.Date(), "%m/%d/%Y"), Sequencing_platform = "Demo data", RNAseq_reference = "Demo",
  kit_number = 1, pos_number = 1, prep_number = 5, spikein_number = 2,
  authorName1 = "Gotools Quickstart", authorEmail1 = NA
)
dirInfo <- Go_dirInfo(Current_working_dir = currentwd, Image_locations = sprintf("%s/pdf/", dir_root), DA_image_location = sprintf("%s/pdf/DA_plot/", dir_root))
tabInfo <- Go_tabInfo(RNAseq = rnaseq_table, Tract_Tab = "1_out/QC_summary.csv")
imgInfo <- Go_imgInfo(RNAseq.heatmap = c("taxa heatmap"), DAplot = da_labels)
sumInfo <- Go_sumInfo(Prep_overview = "Demo dataset: fully synthetic genes, counts, and groups.")

out <- Go_renderReport(family = "RNAseq", project = project, expInfo = expInfo, dirInfo = dirInfo,
                        tabInfo = tabInfo, imgInfo = imgInfo, sumInfo = sumInfo, output_dir = summary_report_dir)

message("Quickstart report rendered: ", out)
