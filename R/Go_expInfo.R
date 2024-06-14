
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

Go_expInfo <- function(Project_name=NA,
                       Samples_info=NA,
                       Sequencing_date=NA,
                       Sequencing_platform=NA,
                       RNAseq_reference=NA,
                       kit_number = NA, pos_number=NA,prep_number = NA, spikein_number = NA, authorName1=NA ,authorName2 = NULL) {
  # Check if all arguments are missing and print options if they are
  if (is.na(kit_number) && is.na(prep_number) && is.na(spikein_number)) {
    cat(
        "Project_name: Add the project name.\n",
        "Samples_info: Add the sample information incouding number of samples. \n",
        "Sequencing_date: Add the sequenicng date.\n",
        "Sequencing_platform: Add the sequenicng platform and sequencing kit. \n\n",
        "Available kit numbers and their descriptions:\n",
        "1: ZymoBIOMICS 96 MagBead DNA/RNA Kit (Catalog# D4300)\n",
        "2: DNeasy® 96 PowerSoil® Pro QIAcube® HT Kit (Catalog# 47021)\n",
        "3: QIAamp® 96 Virus QIAcube® HT (Catalog# 57731)\n",
        "4: QIAamp® 96 DNA QIAcube® HT (Catalog# 51331)\n",
        "5: ZymoBIOMICS DNA Miniprep Kit (Catalog# R2002)\n",
        "6: DNeasy PowerSoil Pro Kit (Catalog# 47016)\n",
        "7: DNeasy DNA and Blood kit (Catalog# 69506)\n\n",
        "8: RNeasy Mini Kit (Catalog# 74104)\n\n",

        "Available positive control numbers and their descriptions:\n",
        "1: ZymoBIOMICS Microbial Community Standard (D6300)\n\n",

        "Available prep numbers and their descriptions:\n",
        "1: Illumina V3V4\n",
        "2: Illumina V1V2\n",
        "3: Zymo Quick-16S/ V3V4 Plus NGS Library prep kit (#catalog D6420)\n",
        "4: Illumina DNA Prep (#catalog 20060059)\n\n",
        "5: QIAseq Stranded RNA Lib Kit UDI ((#catalog 180450)\n",
        "6: QIAseq UPXome RNA Library Kits ((#catalog 334702)\n",
        "7: \n\n",

        "Available spike-in numbers and their descriptions:\n",
        "1: / ZymoBIOMICSTM Spike-in Control II (Low Microbial Load) (D6321)\n",
        "2: No Spike-in Control\n\n",
        "Available Technicians and email address:\n",
        "Heekuk Park : hp2523@cumc.columbia.edu\n",
        "Djamila Eliby : dj2711@cumc.columbia.edu\n",
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
                "8" = "RNeasy Mini Kit (Catalog# 74104)",
                "Unknown kit")


  # Define the positive control based on prep_number
  pos <- switch(as.character(pos_number),
                 "1" = "ZymoBIOMICS Microbial Community Standard (D6300)",
                 "Unknown prep")

  # Define the prep based on prep_number
  prep <- switch(as.character(prep_number),
                 "1" = "Illumina V3V4",
                 "2" = "Illumina V1V2",
                 "3" = "Zymo Quick-16S/ V3V4 Plus NGS Library prep kit (#catalog D6420)",
                 "4" = "Illumina DNA Prep (#catalog 20060059)",
                 "5" = "QIAseq Stranded RNA Lib Kit UDI ((#catalog 180450)",
                 "6" = "QIAseq UPXome RNA Library Kits ((#catalog 334702)",
                 "7" = "",

                 "Unknown prep")

  # Define the spikein based on spikein_number
  spikein <- switch(as.character(spikein_number),
                    "1" = "/ ZymoBIOMICSTM Spike-in Control II (Low Microbial Load) (D6321)",
                    "2" = NULL,  # Representing no spike-in control
                    "Unknown Spike-in")

  # Define email addresses based on author names
  getEmail <- function(name) {
    switch(name,
           "Heekuk Park" = "hp2523@cumc.columbia.edu",
           "Djamila Eliby" = "dj2711@cumc.columbia.edu",
           "Sofia Moscovitz" = "szm2110@cumc.columbia.edu",
           "Dalia Moallem" = "dhm2127@cumc.columbia.edu",
           "Dwayne Seeram" = "ds4057@cumc.columbia.edu",
           "Kristen Lewis" = "kl2954@cumc.columbia.edu",
           NA)
  }

  authorEmail1 <- getEmail(authorName1)
  authorEmail2 <- if (!is.null(authorName2)) {
    getEmail(authorName2)
  } else {
    NULL
  }

  contact_email <- c(authorEmail1, authorEmail2)
  contact_email <- contact_email[!is.na(contact_email) & !is.null(contact_email)]  # Remove NA and NULL values

  # Return a list containing the selected options, including author names and their email addresses
  return(list(
    project = Project_name,
    samples = Samples_info,
    date = Sequencing_date,
    platform = Sequencing_platform,
    kit = kit,
    pos = pos,
    prep = prep,
    spikein = spikein,
    rnaseq_reference=RNAseq_reference,
    authorName1 = authorName1,
    authorName2 = authorName2,
    contact_email = contact_email
  ))
}
