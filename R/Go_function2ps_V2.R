
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
  if (tolower(func.type) %in% c("picrust", "picrust2", "picrustt2")) {
    func.type <- "PICRUSt2"
  } else if (tolower(func.type) %in% c("humann", "humann2", "humann3")) {
    func.type <- "Humann3"

    # HUMAnN 데이터 전처리: 특정 단어 포함된 행 제거
    func.tab <- func.tab[!grepl("ribosomal protein|UNGROUPED|unclassified|UNMAPPED|UNINTEGRATED", rownames(func.tab)), ]

    # KO ID와 설명을 분리
    func.tab1 <- data.frame(do.call(rbind, str_split(rownames(func.tab), ": ", n = 2)))
    colnames(func.tab1) <- c("ID", "Description")

    # Description이 없는 경우 처리
    func.tab1$Description[is.na(func.tab1$Description)] <- "Unknown"

    # 원본 데이터에 병합
    func.tab <- cbind(func.tab1, func.tab)
    rownames(func.tab) <- func.tab$ID
  } else {
    print("Please define the data type: PICRUSt2 or Humann3")
    return(NULL)
  }

  # Define number of samples and split data frame
  NumOfSample <- ncol(func.tab)
  otu <- as.matrix(func.tab[, 3:NumOfSample])
  tax <- as.matrix(func.tab[, 1:2])
  rownames(otu) <- tax[, 1]
  rownames(tax) <- tax[, 1]

  # Define KEGG or pathway
  if (any(grepl("^K\\d{5}", tax[, 1]))) {
    colnames(tax) <- c("KO", "KO.des")
    func <- "KEGG"
  } else if (any(grepl("^PWY", tax[, 1]))) {
    colnames(tax) <- c("pathway", "path.des")
    func <- "pathway"
  } else {
    func <- NULL
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
