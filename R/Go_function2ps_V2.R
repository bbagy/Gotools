
#' Convert Functional Data to Phyloseq Object
#'
#' This function converts functional data from PICRUSt or HUMAnN into a `phyloseq` object for further analysis.
#'
#' @param tabPath The path to the functional data file.
#' @param project (Optional) A string indicating the project name for file naming.
#' @param func.type A string specifying the type of functional data ("PICRUSt" or "HUMAnN").
#' @param name (Optional) A string to add as a suffix to the file name.
#'
#' @details
#' The function reads the functional data (e.g., KO or pathway abundances) and converts it into a format
#' suitable for integration into a `phyloseq` object. The function automatically detects the data type
#' (PICRUSt or HUMAnN) and processes the data accordingly.
#'
#' @return
#' A `phyloseq` object containing the functional data. The function also saves an RDS file of the `phyloseq` object.
#'
#' @examples
#' Go_function2ps(tabPath = "path/to/functional_data.txt",
#'                project = "MyProject",
#'                func.type = "PICRUSt",
#'                name = "Analysis1")
#'
#' @export

Go_function2ps <- function(tabPath, project = NULL, func.type, name = NULL) {
  library(phyloseq)
  library(stringr)

  # Read tab
  func.tab <- try(read.delim(tabPath, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE), silent = TRUE)

  if (inherits(func.tab, "try-error")) {
    print("Error reading the table. Check the file path and format.")
    return(NULL)
  }

  # Identify data type and perform necessary preprocessing
  func.type <- tolower(func.type)

  if (func.type %in% c("picrust", "picrust2", "picrust2")) {
    func.type <- "PICRUSt2"

    # ID (KO, PWY, EC 등) 설정
    if (!is.null(rownames(func.tab))) {
      func.tab$ID <- rownames(func.tab)
    } else if ("function" %in% colnames(func.tab)) {
      func.tab$ID <- func.tab$`function`
    } else {
      stop("Cannot identify function ID column.")
    }

    # Sample matrix와 tax table 분리
    NumOfSample <- ncol(func.tab) - 1
    otu <- as.matrix(func.tab[, 1:NumOfSample])
    tax <- as.matrix(func.tab[, NumOfSample + 1, drop = FALSE])
    rownames(otu) <- func.tab$ID
    rownames(tax) <- func.tab$ID

    # ID가 KO, Pathway, 또는 EC인지 판단
    if (any(grepl("^K\\d{5}$", rownames(tax)))) {
      func <- "KEGG"
      colnames(tax) <- "KO"

    } else if (any(grepl("^PWY", rownames(tax)))) {
      func <- "pathway"
      colnames(tax) <- "Pathway"

    } else if (any(grepl("^\\d+\\.\\d+\\.\\d+\\.\\d+$", rownames(tax)))) {
      func <- "EC"
      colnames(tax) <- "EC"

    } else {
      func <- "Unknown"
      colnames(tax) <- "Function"
    }

  } else if (func.type %in% c("humann", "humann2", "humann3")) {
    func.type <- "Humann3"

    # 특정 단어 포함된 행 제거
    func.tab <- func.tab[!grepl("ribosomal protein|UNGROUPED|unclassified|UNMAPPED|UNINTEGRATED", rownames(func.tab)), ]

    # ID와 Description 분리
    func.tab1 <- data.frame(do.call(rbind, str_split(rownames(func.tab), ": ", n = 2)))
    colnames(func.tab1) <- c("ID", "Description")
    func.tab1$Description[is.na(func.tab1$Description)] <- "Unknown"

    # 병합 및 재구성
    func.tab <- cbind(func.tab1, func.tab)
    rownames(func.tab) <- func.tab$ID

    # Sample matrix와 tax table 분리
    NumOfSample <- ncol(func.tab)
    otu <- as.matrix(func.tab[, 3:NumOfSample])
    tax <- as.matrix(func.tab[, 1:2])
    rownames(otu) <- tax[, 1]
    rownames(tax) <- tax[, 1]

    # KO, PWY, EC 판별
    if (any(grepl("^K\\d{5}", tax[, 1]))) {
      colnames(tax) <- c("KO", "KO.des")
      func <- "KEGG"

    } else if (any(grepl("^PWY", tax[, 1]))) {
      colnames(tax) <- c("pathway", "path.des")
      func <- "pathway"

    } else if (any(grepl("^\\d+\\.\\d+\\.\\d+\\.\\d+$", tax[, 1]))) {
      colnames(tax) <- c("EC", "EC.des")
      func <- "EC"

    } else {
      colnames(tax) <- c("Function", "Description")
      func <- "Unknown"
    }

  } else {
    stop("Please define the correct data type: PICRUSt2 or Humann3")
  }

  # Merge phyloseq
  ps <- phyloseq(otu_table(otu, taxa_are_rows = TRUE), tax_table(tax))

  print(ps)

  # Saving file
  rds <- file.path("2_rds")
  if (!dir.exists(rds)) dir.create(rds, recursive = TRUE)

  saveRDS(ps, file.path(rds, sprintf("ps.%s.%s%s%s%s.rds",
                                     func.type,
                                     ifelse(is.null(func), "", paste0(func, ".")),
                                     ifelse(is.null(project), "", paste0(project, ".")),
                                     ifelse(is.null(name), "", paste0(name, ".")),
                                     format(Sys.Date(), "%y%m%d"))))
}
