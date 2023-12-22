
#' Convert Functional Data to Phyloseq Object
#'
#' @param project Name of the project.
#' @param tool Name of the tool used (e.g., 'humann2' or 'picrust2').
#' @param koTOpath Path to the file containing KEGG Orthology to Pathway mappings.
#' @param koSta Path to the file containing KEGG Orthology stratified data.
#' @param name Optional name for the saved Phyloseq object.
#' @param alpha Threshold value for filtering data based on row sums.
#'
#' @details
#' This function processes the output of functional prediction tools like HUMAnN2 or PICRUSt2 and converts it into a Phyloseq object. It includes steps for data cleaning, subsetting, and reformatting.
#'
#' @return
#' A Phyloseq object representing the processed functional data. This object is also saved as an RDS file.
#'
#' @examples
#' Go_huamnn2ps(project = "MyMicrobiomeProject",
#'              tool = "humann2",
#'              koTOpath = "KO_to_Pathway.csv",
#'              koSta = "KO_Stratified.csv",
#'              name = "FunctionalData",
#'              alpha = 50)
#'
#' @export

Go_huamnn2ps<-function(project, tool, koTOpath, koSta, name, alpha){ # alpha, for rowsum
  # out dir
  out <- file.path("2_rds") 
  if(!file_test("-d", out)) dir.create(out)
  
  print("version2")
  
  # input file for index match file koTOpath
  koTOpath <- read.csv(sprintf("%s",koTOpath),header=T,as.is=T,row.names=1,check.names=F)
  
  # input files kegg-orthology stratified
  if (tool == "humann2"){
    ko.sta <- read.csv(sprintf("%s",koSta), header=T, as.is=T, sep="\t", row.names=1, comment.char="", quote="")
  }else if(tool == "picrust2"){
    ko.sta <- read.csv(sprintf("%s",koSta),header=T,as.is=T,row.names=1,check.names=F)
  }
  

  
  
  
  #ko.sta <- round(1000*ko.sta)
  
  # remove na column
  all_na <- function(x) any(!is.na(x))
  ko.sta.na <- ko.sta %>% select_if(all_na)
  print(sprintf("remove column na %s to %s",dim(ko.sta)[2], dim(ko.sta.na)[2]))
 
  ## rowsum 
  # remove pathways with <50 reads, detected in less than 10 samples, scale to relative abundance
  # inds_to_remove <- which(rowSums(ko.sta.na) < alpha) 
  # ko.sta.na <- ko.sta.na[setdiff(1:nrow(ko.sta.na), inds_to_remove),]
  # print(dim(ko.sta.na))
  
  
  ## 각 수치가 > 2 10개 미만이면 빼기
   n <- round(length(colnames(ko.sta.na))/10);n
   inds_to_remove <- which(rowSums(ko.sta.na > alpha) < n) 
   ko.sta.na <- ko.sta.na[setdiff(1:nrow(ko.sta.na), inds_to_remove),]
  
  print(dim(ko.sta.na))
  # remove less important
  for (rev in c("NO_NAME", "ribosomal protein","UNGROUPED","unclassified","UNMAPPED")){
    ko.sta.na <- subset(ko.sta.na, !(grepl(rev, rownames(ko.sta.na))))
    print(dim(ko.sta.na))
  }
  
  

  if (tool == "humann2"){
    # split names of kolist
    head(rownames(ko.sta.na))
    rownames(ko.sta.na) <- gsub(",", ".", rownames(ko.sta.na));head(rownames(ko.sta.na))
    koList.sta1 <- data.frame(str_split(rownames(ko.sta.na), ": ", simplify = T));head(koList.sta1)
    koList.sta1$X3 <- NULL; head(koList.sta1)
    
    # koList.sta2 <- data.frame(str_split(koList.sta1$X2, "\\|", simplify = T));head(koList.sta2)
    # koList.sta3 <- data.frame(str_split(koList.sta2$X2, "\\.", simplify = T));head(koList.sta3)
    
    # koList.sta1$X2 <- koList.sta2$X1; head(koList.sta1)
    # koList.sta1$X3 <- as.character(koList.sta2$X2); head(koList.sta1)
    # koList.sta1$X3 <- koList.sta3$X1; head(koList.sta1)
    # koList.sta1$X4 <- koList.sta3$X2; head(koList.sta1)
    
    # Add path information
    koList.sta1$Path <- factor(koTOpath$Path[match(koList.sta1$X1, koTOpath$KO)]);head(koList.sta1)
    koList.sta1$Path.des <- factor(koTOpath$Path.description[match(koList.sta1$Path, koTOpath$Path)]);head(koList.sta1)
    
    
    # add header names
    headers <- vector(dim(koList.sta1)[1], mode="character")
    for (i in 1:dim(koList.sta1)[1]) {
      headers[i] <- paste("func", i, sep="_")
    }
    
    rownames(koList.sta1) <- headers
    colnames(koList.sta1) <- c("KO", "KO.des","Path","Path.des");head(koList.sta1)
    koList.sta1$funcNo <- headers;head(koList.sta1)
    rownames(ko.sta.na) <- headers #;head(ko.sta)
    
  }else if(tool == "picrust2"){
    koList.sta1 <- data.frame(rownames(ko.sta.na))
    koList.sta1$KO <- rownames(ko.sta.na);head(koList.sta1)
    koList.sta1$KO.des <- factor(koTOpath$KO.description[match(koList.sta1$KO, koTOpath$KO)]);head(koList.sta1)
    koList.sta1$Path <- factor(koTOpath$Path[match(rownames(ko.sta.na), koTOpath$KO)]);head(koList.sta1)
    koList.sta1$Path.des <- factor(koTOpath$Path.description[match(koList.sta1$Path, koTOpath$Path)]);head(koList.sta1)
    
    # add header names
    rownames(koList.sta1) <- koList.sta1$rownames.ko.sta.na.
    koList.sta1$rownames.ko.sta.na. <-NULL
    #colnames(koList.sta1) <- c("KO", "KO.des","Path","Path.des");head(koList.sta1)
  }
  

  
  # save information
  write.csv(koList.sta1, quote = F, col.names = F, row.names=T,
            file=sprintf("%s/koList.sta.%s.csv", out, format(Sys.Date(), "%y%m%d"),sep="\t"))
  
  # column 정리 
  #colnames(ko.sta)

  
  #--- create phyloseq file ---#
  
  otu.sta <- as.matrix(ko.sta.na)
  tax.sta <- as.matrix(koList.sta1)
  
  
  OTU.sta <- otu_table(otu.sta, taxa_are_rows = TRUE);dim(OTU.sta)
  TAX.sta <- tax_table(tax.sta);dim(TAX.sta)



  #ps1
  kps1.sta <- phyloseq(otu_table(OTU.sta, taxa_are_rows=FALSE), tax_table(TAX.sta))
  
  print(kps1.sta)
  if (length(name) == 1) {
    saveRDS(kps1.sta, sprintf("%s/ps.funtion.%s.%s.%s.rds", out, project, name, format(Sys.Date(), "%y%m%d")))
  } else{
    saveRDS(kps1.sta, sprintf("%s/ps.funtion.%s.%s.rds", out, project,format(Sys.Date(), "%y%m%d")))
  }
}
