
#' Convert Functional Annotation Data to Phyloseq Object
#'
#' @param tabPath The path to the functional annotation table file.
#' @param project Optional project name for labeling purposes.
#' @param func.type The type of functional annotation, either PICRUSt or HUMAnN.
#' @param name Optional name for the output object.
#'
#' @details
#' This function converts a functional annotation table into a Phyloseq object for further analysis in ecological and bioinformatic contexts. It supports both PICRUSt and HUMAnN annotation types. The function reads the annotation table, cleans it by removing certain unwanted categories, and structures it into a format compatible with Phyloseq.
#'
#' @return
#' A Phyloseq object containing the functional annotation data, ready for ecological and bioinformatic analyses.
#'
#' @examples
#' ps_object <- Go_function2ps(tabPath = "path/to/annotation_table.csv",
#'                             project = "MyProject",
#'                             func.type = "PICRUSt",
#'                             name = "FunctionalAnalysis")
#'
#' @export

Go_function2ps <- function(tabPath,
                           project=NULL, func.type,
                           name=NULL){
  
  tt<- try(class(func.type),T)
  
  if(class(tt) == "try-error"){
    print("Please define the data type. PICRUSt or Hummann")
    break
  }else if(any(grepl(func.type, c("picrust","Picrust","Picrust2","PICRUSt","PICRUSTt2")))){
    func.type <- "PICRUSTt2"
    # Read tab
    func.tab <- read_tsv(tabPath,col_types = cols())
  }else if(any(grepl(func.type, c("Human","Humann","humann","human","Human2","Humann2","humann2","human2")))){
    func.type <- "Humann2"
    
    func.tab <- read.csv(sprintf("%s",tabPath), header=T, as.is=T, sep="\t", row.names=1, comment.char="", quote="")
    
    for (rev in c("ribosomal protein","UNGROUPED","unclassified","UNMAPPED","UNINTEGRATED")){ # ,"NO_NAME"
      func.tab <- subset(func.tab, !(grepl(rev, rownames(func.tab))))
      print(dim(func.tab))
    }
    
    
    rownames(func.tab) <- gsub(",", ".", rownames(func.tab));head(rownames(func.tab))
    func.tab1 <- data.frame(str_split(rownames(func.tab), ": ", simplify = T));head(func.tab)
    func.tab$X3 <- NULL; head(func.tab)
    
    colnames(func.tab)
    
    
    rownames(func.tab) <- func.tab1$X1
    
    func.tab <- cbind(description =func.tab1$X2, func.tab)
    
  }
  
  NumOFsample <- dim(func.tab)[2]
  
  # split a data frame
  otu <- as.matrix(func.tab[,3:NumOFsample])
  tax <- as.matrix(func.tab1[,1:2])
  
  tax.df <- data.frame(tax)
  rownames(tax.df) <- tax.df$X1
  tax <- as.matrix(tax.df)
  
  if (func.type == "PICRUSTt2"){
    rownames(otu) <- tax[,1];dim(otu)
    rownames(tax) <- tax[,1];dim(tax)
  }
  
  
  
  # define kegg or pathway
  if (any(grepl("K0", rownames(tax)))){
    colnames(tax) <- c("KO","KO.des")
    func <- "KEGG"
  }else if (any(grepl("PWY", rownames(tax)))){
    colnames(tax) <- c("pathway","path.des")
    func <- "pathway"
  }else{
    func <- NULL
  }
  
  #merge phyloseq
  ps <- phyloseq(otu_table(otu, taxa_are_rows=T), tax_table(tax));ps
  
  print(ps)
  
  
  # saving file
  rds <- file.path("2_rds") 
  if(!file_test("-d", rds)) dir.create(rds)
  saveRDS(ps, sprintf("%s/ps.%s.%s.%s%s%s.rds",rds,
                      func.type,
                      project,
                      ifelse(is.null(func), "", paste(func, ".", sep = "")), 
                      ifelse(is.null(name), "", paste(name, ".", sep = "")), 
                      format(Sys.Date(), "%y%m%d")))
}
