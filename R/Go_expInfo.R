
#' Go_expInfo - Experiment Information Retrieval Function
#'
#' This function retrieves experiment information based on provided kit, prep, and spike-in numbers.
#' If no arguments are provided, it prints available options for each parameter.
#'
#' @param kit_number An optional integer parameter for selecting the kit.
#' Available kit numbers are:
#' 1. ZymoBIOMICS 96 MagBead DNA/RNA Kit (Catalog# D4300)
#' 2. DNeasy® 96 PowerSoil® Pro QIAcube® HT Kit (Catalog# 47021)
#' 3. QIAamp® 96 Virus QIAcube® HT (Catalog# 57731)
#' 4. QIAamp® 96 DNA QIAcube® HT (Catalog# 51331)
#' 5. ZymoBIOMICS DNA Miniprep Kit (Catalog# R2002)
#' 6. DNeasy PowerSoil Pro Kit (Catalog# 47016)
#' 7. DNeasy DNA and Blood kit (Catalog# 69506)
#'
#' @param prep_number An optional integer parameter for selecting the prep method.
#' Available prep numbers are:
#' 1. Illumina V3V4
#' 2. Illumina V1V2
#' 3. Zymo Quick-16S/ V3V4 Plus NGS Library prep kit (#catalog D6420)
#'
#' @param spikein_number An optional integer parameter for selecting the spike-in control.
#' Available spike-in numbers are:
#' 1. / ZymoBIOMICSTM Spike-in Control II (Low Microbial Load) (D6321)
#' 2. No Spike-in Control
#'
#' @return A list containing the selected kit, prep method, and spike-in control based on the input parameters.
#' If no arguments are provided, information about available options is printed and the function returns `NULL`.
#'
#' @examples
#' Go_expInfo() # prints available options
#' exp_info <- Go_expInfo(kit_number = 6, prep_number = 1, spikein_number = 2)
#' print(exp_info)

Go_expInfo <- function(kit_number = NA, prep_number = NA, spikein_number = NA) {
  # Check if all arguments are missing and print options if they are
  if (is.na(kit_number) && is.na(prep_number) && is.na(spikein_number)) {
    cat("Available kit numbers and their descriptions:\n",
        "1: ZymoBIOMICS 96 MagBead DNA/RNA Kit (Catalog# D4300)\n",
        "2: DNeasy® 96 PowerSoil® Pro QIAcube® HT Kit (Catalog# 47021)\n",
        "3: QIAamp® 96 Virus QIAcube® HT (Catalog# 57731)\n",
        "4: QIAamp® 96 DNA QIAcube® HT (Catalog# 51331)\n",
        "5: ZymoBIOMICS DNA Miniprep Kit (Catalog# R2002)\n",
        "6: DNeasy PowerSoil Pro Kit (Catalog# 47016)\n",
        "7: DNeasy DNA and Blood kit (Catalog# 69506)\n\n",
        "Available prep numbers and their descriptions:\n",
        "1: Illumina V3V4\n",
        "2: Illumina V1V2\n",
        "3: Zymo Quick-16S/ V3V4 Plus NGS Library prep kit (#catalog D6420)\n\n",
        "Available spike-in numbers and their descriptions:\n",
        "1: / ZymoBIOMICSTM Spike-in Control II (Low Microbial Load) (D6321)\n",
        "2: No Spike-in Control\n",
        "Available Technicians and email address:\n",
        "Sofia Moscovitz : szm2110@cumc.columbia.edu\n",
        "Dalia Moallem : dhm2127@cumc.columbia.edu\n",
        "Dwayne Seeram : ds4057@cumc.columbia.edu\n",
        "Kristen Lewis : kl2954@cumc.columbia.edu\n")
    return(invisible())
  }
  
  # Define the kit based on kit_number
  kit <- switch(as.character(kit_number),
                "1" = "ZymoBIOMICS 96 MagBead DNA/RNA Kit (Catalog# D4300)",
                "2" = "DNeasy® 96 PowerSoil® Pro QIAcube® HT Kit (Catalog# 47021)",
                "3" = "QIAamp® 96 Virus QIAcube® HT (Catalog# 57731)",
                "4" = "QIAamp® 96 DNA QIAcube® HT (Catalog# 51331)",
                "5" = "ZymoBIOMICS DNA Miniprep Kit (Catalog# R2002)",
                "6" = "DNeasy PowerSoil Pro Kit (Catalog# 47016)",
                "7" = "DNeasy DNA and Blood kit (Catalog# 69506)",
                "Unknown kit")
  
  # Define the prep based on prep_number
  prep <- switch(as.character(prep_number),
                 "1" = "Illumina V3V4",
                 "2" = "Illumina V1V2",
                 "3" = "Zymo Quick-16S/ V3V4 Plus NGS Library prep kit (#catalog D6420)",
                 "Unknown prep")
  
  # Define the spikein based on spikein_number
  spikein <- switch(as.character(spikein_number),
                    "1" = "/ ZymoBIOMICSTM Spike-in Control II (Low Microbial Load) (D6321)",
                    "2" = NULL,  # Representing no spike-in control
                    "Unknown Spike-in")
  
  # Return a list containing the selected options
  return(list(kit = kit, prep = prep, spikein = spikein))
}