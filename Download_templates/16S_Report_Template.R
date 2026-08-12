############################################################
# 16S Summary Report Template
#
# Order to follow:
#   1. Go_expInfo()   - project / sequencing metadata
#   2. Go_dirInfo()   - where images live
#   3. Go_tabInfo()   - ASV table, track table, alpha-div table
#   4. Go_imgInfo()   - which plots to pull into the report
#   5. Go_sumInfo()   - narrative text for each report section
#   6. Go_renderReport() - renders the bundled 16S Rmd to HTML
#
# Run each Go_*Info() with no arguments to print its option guide,
# e.g. Go_expInfo() alone in the console lists kit/prep/spikein numbers.
#
# Based on: 20260706_MAPS_Zach_16S_CRC_GAB_analysis.R
############################################################

library(Gotools)
Gotool_dependency()

project   <- "FILL_IN_PROJECT_NAME"
currentwd <- "FILL_IN_ANALYSIS_DIRECTORY"
setwd(currentwd)

dir_out        <- sprintf("%s_%s", project, format(Sys.Date(), "%y%m%d"))
asv_table_file <- "1_out/FILL_IN.asvTable.csv"
track_file     <- "1_out/FILL_IN.track.csv"
alpha_tab_file <- sprintf("%s/table/adiv/adiv.%s.FILL_IN_COMPARISON.%s.csv",
                          dir_out, project, format(Sys.Date(), "%y%m%d"))

barchart_tax_ranks <- c("Phylum", "Class", "Order", "Family", "Genus", "Species")


############################################################
#=================   Report Objects   =================#
############################################################

expInfo <- Go_expInfo(
  Project_name        = "FILL_IN human-readable study title",
  Samples_info         = "FILL_IN e.g. Total 92 samples; 47 Group A, 45 Group B",
  Sequencing_date      = format(Sys.Date(), "%m/%d/%Y"),
  Sequencing_platform  = "FILL_IN e.g. Illumina MiSeq V3 600 cycles",
  kit_number           = NA, # DNeasy PowerSoil Pro Kit = 6, see Go_expInfo() guide
  pos_number           = NA, # ZymoBIOMICS Microbial Community Standard = 1
  prep_number          = NA, # Illumina V3V4 = 1
  spikein_number       = NA, # No Spike-in Control = 2
  authorName1          = "Heekuk Park",
  authorName2          = ""
)

dirInfo <- Go_dirInfo(
  Current_working_dir = currentwd,
  Image_locations     = sprintf("%s/pdf/", dir_out),
  DA_image_location   = sprintf("%s/pdf/DA_plot/", dir_out)
)

tabInfo <- Go_tabInfo(
  ASVs_Tab          = asv_table_file,
  Tract_Tab         = track_file,
  Alpha_divTab      = alpha_tab_file,
  Alpha_div_LmerTab = NULL,
  PermanovaTab      = NULL
)

imgInfo <- Go_imgInfo(
  Rarefaction = NA,
  Overview    = NA,
  Barchart    = barchart_tax_ranks
  # add Adivplot / Bdivplot / DAplot tokens once those analyses are run
)

sumInfo <- Go_sumInfo(
  Prep_overview = "FILL_IN - extraction/sequencing QC summary, note any failed samples",
  Taxa_overview = "FILL_IN - dominant taxa observed in the barchart"
  # Alpha_div, Beta_div, DA_test, Conclusion optional - add as those analyses are run
)


############################################################
#=================   Render HTML Report   =================#
############################################################

Go_renderReport(
  family     = "16S",
  project    = project,
  expInfo    = expInfo,
  dirInfo    = dirInfo,
  tabInfo    = tabInfo,
  imgInfo    = imgInfo,
  sumInfo    = sumInfo,
  output_dir = file.path(currentwd, "3_Summary_report")
)
