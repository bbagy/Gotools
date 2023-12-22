
#' Convert RPKM Data to Phyloseq Object
#'
#' This function converts Reads Per Kilobase of transcript, per Million mapped reads (RPKM) data to a Phyloseq object.
#'
#' @param project A string indicating the project name used for file naming.
#' @param speciesTab A data frame containing species information.
#' @param genoemSizeTab A data frame containing genome sizes.
#' @param taxaTab A data frame containing taxonomic information.
#'
#' @details
#' The function cleans and processes species names, matches genome sizes, and calculates RPKM values.
#' It then integrates this information into a Phyloseq object, which is useful for microbiome data analysis.
#'
#' @return
#' A Phyloseq object containing RPKM data structured for microbiome analysis.
#'
#' @examples
#' Go_rpkmToPs(project = "MyProject",
#'             speciesTab = species_dataframe,
#'             genoemSizeTab = genome_size_dataframe,
#'             taxaTab = taxonomic_dataframe)
#'
#' @export

Go_rpkmToPs <- function(project, speciesTab, genoemSizeTab, taxaTab) {
  # out path
  tem_path <- file.path(sprintf("%s", "RPKM_tem"))
  if(!file_test("-d", tem_path)) dir.create(tem_path)
  rds_path <- file.path(sprintf("%s", "2_rds"))
  if(!file_test("-d", rds_path)) dir.create(rds_path)
  

  bacTab <- speciesTab
  
  #-- get simple species and genus names --#  species 이름이 길다.
  bacTab$Species <- rownames(bacTab)
  cleanSpecies <- data.frame(do.call('rbind', strsplit(as.character(bacTab$Species),'_',fixed=TRUE)))
  cleanSpecies$Species <- paste(cleanSpecies$X1, cleanSpecies$X2, sep="_")
  bacTab$Genus <- cleanSpecies$X1
  bacTab$Species <- cleanSpecies$Species
  bacTab$SubSpecies <- rownames(bacTab)
  
  # read genome size table
  gSizeTab <- genoemSizeTab
  colnames(gSizeTab) <- c("AccessionVersion", "TaxId", "Completeness", "GenomeLen", "Title")
  cleanSpecies_gSizeTab <- data.frame(do.call('rbind', strsplit(as.character(gSizeTab$Title),' ',fixed=TRUE)))
  cleanSpecies_gSizeTab$Species <- paste(cleanSpecies_gSizeTab$X1, cleanSpecies_gSizeTab$X2, sep="_")
  cleanSpecies_gSizeTab$Genus <- cleanSpecies_gSizeTab$X1
  gSizeTab$Genus <- cleanSpecies_gSizeTab$X1
  gSizeTab$Species <- cleanSpecies_gSizeTab$Species
  
  
  # create new data frame and exrtact genome size for calculation
  genomeFinal <- data.frame(cbind(bacTab$Species, bacTab$Genus)) # 두번째 col은 이름이 안따라온다. 다음 명령어를 사용해야 한다.
  colnames(genomeFinal) <- c("Species","Genus")
  genomeFinal$Genus <- bacTab$Genus
  genomeFinal$GenomeLen <- gSizeTab$GenomeLen[match(bacTab$Genus, gSizeTab$Genus)] # index match
  
  
  write.csv(genomeFinal, quote = FALSE, #col.names = NA, #row.names = FALSE,
            file=sprintf("%s/%s.genomeFinal2.%s.csv", tem_path, project, format(Sys.Date(), "%y%m%d"),sep="/"))
  
  
  bacTab$Genus <- NULL
  bacTab$Species <- NULL
  bacTab$SubSpecies <- NULL
  ## RPKM
  geneLength <- 1
  bacNor <- data.frame(sapply(bacTab , function(column) 10^9 * column / geneLength / sum(column)))
  rownames(bacNor) <- rownames(bacTab); head(rownames(bacNor))
  
  #brac.rpkmS <- bracNor*100/genomeFinal$SlenS
  #brac.rpkmS <- ceiling(brac.rpkmS[-c(99),]) # 반올림
  #brac.rpkmS[is.na(brac.rpkmS)] <- 0 # remove NA
  
  bacTab.rpkmG <- bacNor*100/genomeFinal$GenomeLen
  bacTab.rpkmG <- ceiling(bacTab.rpkmG[-c(99),]) # 반올림
  bacTab.rpkmG[is.na(bacTab.rpkmG)] <- 0 # remove NA
  
  ## add LCA
  cleanSpecies <- data.frame(do.call('rbind', strsplit(as.character(rownames(bacTab.rpkmG)),'_',fixed=TRUE)))
  cleanSpecies$Species <- paste(cleanSpecies$X1, cleanSpecies$X2, sep="_")
  bacTab.rpkmG$Species <- cleanSpecies$Species
  bacTab.rpkmG$Genus <- cleanSpecies$X1
  
  Taxa <- taxaTab
  
  ranks <- c("Rank1", "Phylum", "Class", "Order", "Family")
  for(taxa in ranks){
    bacTab.rpkmG[,taxa]  <- Taxa[,taxa][match(bacTab.rpkmG$Genus, Taxa$Genus)]
  }
  
  bacTab.rpkmG.sel <- subset(bacTab.rpkmG, Rank1 == "Bacteria")
  colnames(bacTab.rpkmG.sel) <-  gsub("X", "", colnames(bacTab.rpkmG.sel));head(colnames(bacTab.rpkmG.sel))
  
  ##########################
  #---      Make ps     ---#
  ##########################
  #------------- tax table -------------#
  headers <- rownames(bacTab.rpkmG.sel);head(headers)
  
  tax <- data.frame(bacTab.rpkmG.sel$Rank1, bacTab.rpkmG.sel$Phylum, bacTab.rpkmG.sel$Class, bacTab.rpkmG.sel$Order, bacTab.rpkmG.sel$Family, bacTab.rpkmG.sel$Genus, bacTab.rpkmG.sel$Species)
  rownames(tax) <- headers
  colnames(tax) <- c("Kingdom","Phylum","Class", "Order", "Family", "Genus", "Species") ;tax
  
  #------------- otu table -------------#
  for(i in colnames(tax)){
    bacTab.rpkmG.sel$Rank1 <- NULL
    bacTab.rpkmG.sel[,i] <- NULL
  }
  
  
  otu <- bacTab.rpkmG.sel
  rownames(otu) <- headers
  otu <- otu[rowSums(otu[, -1])>0, ]
  
  write.csv(otu, quote = FALSE, col.names = NA, #row.names = FALSE,
            file=sprintf("%s/otu_%s.genome_size.%s.csv", tem_path, project, format(Sys.Date(), "%y%m%d"),sep="/"))
  
  
  #--- create phyloseq file ---#
  tax<-as.matrix(tax)
  otu<-as.matrix(otu)
  
  OTU <- otu_table(otu, taxa_are_rows = TRUE);head(OTU)
  TAX <- tax_table(tax);head(TAX)
  ps1 <- phyloseq(otu_table(OTU, taxa_are_rows=FALSE), tax_table(TAX))
  saveRDS(ps1, sprintf("%s/ps_rpkmG.%s.%s.rds",rds_path, project, format(Sys.Date(), "%y%m%d")))
  
  return(ps1)
}




