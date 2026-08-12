############################################################
# RNAseq Summary Report Template
#
# Order to follow:
#   1. Go_expInfo()   - project / sequencing metadata
#   2. Go_dirInfo()   - where images live
#   3. Go_tabInfo()   - RNAseq results table, track table
#   4. Go_imgInfo()   - which plots to pull into the report
#   5. Go_sumInfo()   - narrative text for each report section
#   6. Go_renderReport() - renders the bundled RNAseq Rmd to HTML
#
# Run each Go_*Info() with no arguments to print its option guide,
# e.g. Go_expInfo() alone in the console lists kit/prep/spikein numbers.
#
# NOTE: no 2026 real analysis has run the RNAseq report path yet, so
# this template mirrors GoAutoBot's build_gab_rnaseq_report_info()
# convention (GAB_report_helpers.R) rather than a completed real run.
# Update this note once a real RNAseq report has been generated and
# this template has been checked against it.
############################################################

library(Gotools)
Gotool_dependency()

project   <- "FILL_IN_PROJECT_NAME"
currentwd <- "FILL_IN_ANALYSIS_DIRECTORY"
setwd(currentwd)

dir_out       <- sprintf("%s_%s", project, format(Sys.Date(), "%y%m%d"))
rnaseq_table  <- "FILL_IN path to DESeq2/edgeR results table"
track_file    <- "FILL_IN path to sequencing QC track table"


############################################################
#=================   Report Objects   =================#
############################################################

expInfo <- Go_expInfo(
  Project_name        = "FILL_IN human-readable study title",
  Samples_info         = "FILL_IN e.g. 12 RNAseq samples",
  Sequencing_date      = format(Sys.Date(), "%m/%d/%Y"),
  Sequencing_platform  = "FILL_IN e.g. Illumina NovaSeq",
  RNAseq_reference     = "FILL_IN reference genome/transcriptome build",
  kit_number           = 1, # ZymoBIOMICS 96 MagBead DNA/RNA Kit
  pos_number           = 1, # ZymoBIOMICS Microbial Community Standard
  prep_number          = 5, # Zymo Quick-ITS Plus NGS Library Prep Kit
  spikein_number       = 2, # No Spike-in Control
  authorName1          = "Heekuk Park",
  authorName2          = ""
)

dirInfo <- Go_dirInfo(
  Current_working_dir = currentwd,
  Image_locations     = sprintf("%s/pdf/", dir_out),
  DA_image_location   = sprintf("%s/pdf/DA_plot/", dir_out)
)

tabInfo <- Go_tabInfo(
  RNAseq            = rnaseq_table,
  Tract_Tab         = track_file,
  Alpha_divTab      = NA,
  Alpha_div_LmerTab = NULL,
  PermanovaTab      = NULL
)

imgInfo <- Go_imgInfo(
  Overview       = NA,
  Rarefaction    = NA,
  Barchart       = NA,
  Bac.heatmap    = NA,
  RNAseq.heatmap = c("taxa heatmap")
  # add DAplot tokens once differential expression is run
)

sumInfo <- Go_sumInfo(
  Prep_overview = "FILL_IN - extraction/sequencing QC summary, note any failed samples"
  # Taxa_overview, Alpha_div, Beta_div, DA_test, Conclusion optional
)


############################################################
#=================   Render HTML Report   =================#
############################################################

Go_renderReport(
  family     = "RNAseq",
  project    = project,
  expInfo    = expInfo,
  dirInfo    = dirInfo,
  tabInfo    = tabInfo,
  imgInfo    = imgInfo,
  sumInfo    = sumInfo,
  output_dir = file.path(currentwd, "3_Summary_report")
)
