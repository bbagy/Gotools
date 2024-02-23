
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

version2;
Go_function2ps <- function(tabPath, project = NULL, func.type, name = NULL) {
  # Read tab
  func.tab <- try(read_tsv(tabPath, col_types = cols()), silent = TRUE)

  if (inherits(func.tab, "try-error")) {
    print("Error reading the table. Check the file path and format.")
    return(NULL)
  }

  # Identify data type and perform necessary preprocessing
  if (any(grepl(func.type, c("picrust", "Picrust", "Picrust2", "PICRUSt", "PICRUSTt2")))) {
    func.type <- "PICRUSTt2"
  } else if (any(grepl(func.type, c("Human", "Humann", "humann", "human", "Human2", "Humann2", "humann2", "human2")))) {
    func.type <- "Humann2"
    # Additional preprocessing specific to Humann2
    func.tab <- func.tab[!grepl("ribosomal protein|UNGROUPED|unclassified|UNMAPPED|UNINTEGRATED", rownames(func.tab)), ]
    rownames(func.tab) <- gsub(",", ".", rownames(func.tab))
    func.tab1 <- data.frame(str_split(rownames(func.tab), ": ", simplify = T))
    func.tab <- cbind(description = func.tab1$X2, func.tab)
    rownames(func.tab) <- func.tab1$X1
  } else {
    print("Please define the data type: PICRUSt or Humann")
    return(NULL)
  }

  # Define number of samples and split data frame
  NumOfSample <- dim(func.tab)[2]
  otu <- as.matrix(func.tab[, 3:NumOfSample])
  tax <- as.matrix(func.tab[, 1:2])
  rownames(otu) <- tax[, 1]
  rownames(tax) <- tax[, 1]

  # Define KEGG or pathway
  if (any(grepl("K0", rownames(tax)))) {
    colnames(tax) <- c("KO", "KO.des")
    func <- "KEGG"
  } else if (any(grepl("PWY", rownames(tax)))) {
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
  if (!file_test("-d", rds)) dir.create(rds)
  saveRDS(ps, sprintf("%s/ps.%s.%s.%s%s%s.rds", rds,
                      func.type,
                      ifelse(is.null(func), "", paste(func, ".", sep = "")),
                      ifelse(is.null(project), "", paste(project, ".", sep = "")),
                      ifelse(is.null(name), "", paste(name, ".", sep = "")),
                      format(Sys.Date(), "%y%m%d")))
}

