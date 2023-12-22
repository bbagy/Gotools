
#' Stratified Conversion of KEGG Orthology to Phyloseq Object
#'
#' This function converts KEGG Orthology data into a phyloseq object with stratification and filtering.
#'
#' @param project A string representing the project name, used for file naming and directory creation.
#' @param koTOpath The file path for the KEGG Orthology to Pathway mapping data.
#' @param koSta The file path for the KEGG Orthology stratified data.
#' @param alpha A threshold for filtering pathways based on read count.
#' @param name (Optional) A name for labeling the output.
#'
#' @details
#' The function reads KEGG Orthology data, applies filtering based on the provided threshold, and
#' converts it into a phyloseq object. It includes additional information about the pathways and
#' stratifies the data for in-depth analysis.
#'
#' @return
#' Generates a phyloseq object from stratified KEGG Orthology data and saves it as an RDS file.
#'
#' @examples
#' # Example usage:
#' Go_huamnn2ps_stratified(project = "MyProject",
#'                        koTOpath = "path_to_ko_to_path.csv",
#'                        koSta = "path_to_ko_stratified.csv",
#'                        alpha = 50,
#'                        name = "MyAnalysis")
#'
#' @export

Go_huamnn2ps_stratified<-function(project, koTOpath, koSta,alpha, name){
  # out dir
  out <- file.path("2_rds") 
  if(!file_test("-d", out)) dir.create(out)
  
  # input file for index match file koTOpath
  koTOpath <- read.csv(sprintf("%s",koTOpath),header=T,as.is=T,row.names=1,check.names=F)
  
  # input files kegg-orthology stratified
  ko.sta <- read.csv(sprintf("%s",koSta), header=T, as.is=T, sep="\t", row.names=1, comment.char="", quote="")
  
  #ko.sta <- round(1000*ko.sta)
  
  # remove na column
  all_na <- function(x) any(!is.na(x))
  ko.sta.na <- ko.sta %>% select_if(all_na)
  print(sprintf("remove column na %s to %s",dim(ko.sta)[2], dim(ko.sta.na)[2]))

  # remove pathways with <50 reads, detected in less than 10 samples, scale to relative abundance

  
  inds_to_remove <- which(rowSums(ko.sta.na) < alpha) 
  ko.sta.na <- ko.sta.na[setdiff(1:nrow(ko.sta.na), inds_to_remove),]
  print(dim(ko.sta.na))
  
  n <- round(length(colnames(ko.sta.na))/10)
  inds_to_remove <- which(rowSums(ko.sta.na > alpha) < n) # 각 수치가 > 2 10개 미만이면 빼기
  ko.sta.na <- ko.sta.na[setdiff(1:nrow(ko.sta.na), inds_to_remove),]
  
  print(dim(ko.sta.na))
  
  # remove less important #,"unclassified"
  for (rev in c("NO_NAME", "ribosomal protein","UNGROUPED")){
    ko.sta.na <- subset(ko.sta.na, !(grepl(rev, rownames(ko.sta.na))))
    print(dim(ko.sta.na))
  }
  
  
#view(head(ko.sta.na))
  
  # split names of kolist
  head(rownames(ko.sta.na))
  rownames(ko.sta.na) <- gsub(",", ".", rownames(ko.sta.na));head(rownames(ko.sta.na))
  koList.sta1 <- data.frame(str_split(rownames(ko.sta.na), ": ", simplify = T));head(koList.sta1)
  koList.sta2 <- data.frame(str_split(koList.sta1$X2, "\\|", simplify = T));head(koList.sta2)
  koList.sta3 <- data.frame(str_split(koList.sta2$X2, "\\.", simplify = T));head(koList.sta3)
  
  koList.sta1$X2 <- koList.sta2$X1; head(koList.sta1)
  koList.sta1$X3 <- as.character(koList.sta2$X2); head(koList.sta1)
  koList.sta1$X3 <- koList.sta3$X1; head(koList.sta1)
  koList.sta1$X4 <- koList.sta3$X2; head(koList.sta1)
  
  # Add path information
  koList.sta1$Path <- factor(koTOpath$Path[match(koList.sta1$X1, koTOpath$KO)]);head(koList.sta1)
  koList.sta1$Path.des <- factor(koTOpath$Path.description[match(koList.sta1$Path, koTOpath$Path)]);head(koList.sta1)
  
  # add header names
  headers <- vector(dim(koList.sta1)[1], mode="character")
  for (i in 1:dim(koList.sta1)[1]) {
    headers[i] <- paste("func", i, sep="_")
  }
  
  rownames(koList.sta1) <- headers
  colnames(koList.sta1) <- c("KO", "KO.des","Genus","Species","Path","Path.des");head(koList.sta1)
  koList.sta1$funcNo <- headers;head(koList.sta1)
  
  # save information
  write.csv(koList.sta1, quote = F, col.names = F, row.names=T,
            file=sprintf("%s/koList.sta.%s.csv", out, format(Sys.Date(), "%y%m%d"),sep="\t"))
  
  # column 정리 
  #colnames(ko.sta)
  rownames(ko.sta.na) <- headers#;head(ko.sta)
  
  #--- create phyloseq file ---#
  tax.sta <- as.matrix(koList.sta1)
  otu.sta <- as.matrix(ko.sta.na)
  
  OTU.sta <- otu_table(otu.sta, taxa_are_rows = TRUE);head(OTU.sta)
  TAX.sta <- tax_table(tax.sta);head(TAX.sta)
  
  #ps1
  kps1.sta <- phyloseq(otu_table(OTU.sta, taxa_are_rows=FALSE), tax_table(TAX.sta))
  print(kps1.sta)
  if (length(name) == 1) {
    saveRDS(kps1.sta, sprintf("%s/ps.sta.funtion.%s.%s.%s.rds", out, project, name, format(Sys.Date(), "%y%m%d")))
  } else{
    saveRDS(kps1.sta, sprintf("%s/ps.sta.funtion.%s.%s.rds", out, project,format(Sys.Date(), "%y%m%d")))
  }
}
