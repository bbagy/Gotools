############################################################
# MG (Metagenomics) Summary Report Template
#
# Order to follow:
#   1. Go_expInfo()   - project / sequencing metadata
#   2. Go_dirInfo()   - where images live
#   3. Go_tabInfo()   - Kraken2/Bracken taxa table, track table, alpha-div, HUMAnN3 func tables
#   4. Go_imgInfo()   - which plots to pull into the report
#   5. Go_sumInfo()   - narrative text for each report section
#   6. Go_renderReport() - renders the bundled MG Rmd to HTML
#
# Run each Go_*Info() with no arguments to print its option guide,
# e.g. Go_expInfo() alone in the console lists kit/prep/spikein numbers.
#
# Based on: 20260421_DEAPIM30_NV_MG_GAB_analysis.R (Kraken2 + HUMAnN3)
############################################################

library(Gotools)
Gotool_dependency()

project   <- "FILL_IN_PROJECT_NAME"
currentwd <- "FILL_IN_ANALYSIS_DIRECTORY"
setwd(currentwd)

dir_out        <- sprintf("%s_%s", project, format(Sys.Date(), "%y%m%d"))
bracken_file   <- "1_out/FILL_IN_bracken_mpa_filled.txt"
track_file     <- "1_out/FILL_IN_kraken2_log.xlsx"
alpha_tab_file <- sprintf("%s/table/adiv/adiv.%s.FILL_IN_COMPARISON.%s.csv",
                          dir_out, project, format(Sys.Date(), "%y%m%d"))

# HUMAnN3 functional tables - set humann_enabled <- FALSE if not run
humann_enabled            <- TRUE
humann_genefamilies_file  <- "humann3_final_out/kegg-orthology/FILL_IN_merged_genefamilies_ko_cpm.txt"
humann_pathabundance_file <- "humann3_final_out/FILL_IN_merged_pathabundance.txt"


############################################################
#=================   Report Objects   =================#
############################################################

expInfo <- Go_expInfo(
  Project_name        = "FILL_IN human-readable study title",
  Samples_info         = "FILL_IN e.g. 17 human stool samples",
  Sequencing_date      = format(Sys.Date(), "%m/%d/%Y"),
  Sequencing_platform  = "FILL_IN e.g. Illumina NovaSeq / Element AVITI 2x150",
  kit_number           = 2, # DNeasy PowerSoil Pro Kit
  pos_number           = 1, # ZymoBIOMICS Microbial Community Standard
  prep_number          = 6, # Illumina DNA Prep
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
  Taxa_Tab          = bracken_file,
  Func_Tab          = if (humann_enabled) c(humann_genefamilies_file, humann_pathabundance_file) else NA,
  Tract_Tab         = track_file,
  Alpha_divTab      = alpha_tab_file,
  Alpha_div_LmerTab = NULL,
  PermanovaTab      = NULL
)

imgInfo <- Go_imgInfo(
  Rarefaction = NA,
  Overview    = NA,
  Barchart    = c("Phylum", "Class", "Order", "Family", "Genus", "Species")
  # add Adivplot / Bdivplot / DAplot / EBplot / SHeatmap tokens once those analyses are run
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
  family     = "MG",
  project    = project,
  expInfo    = expInfo,
  dirInfo    = dirInfo,
  tabInfo    = tabInfo,
  imgInfo    = imgInfo,
  sumInfo    = sumInfo,
  output_dir = file.path(currentwd, "3_Summary_report")
)
