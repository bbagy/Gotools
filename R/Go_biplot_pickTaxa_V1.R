
#' Select Best Fitting Taxa for Biplot Analysis
#'
#' This function identifies the best combination of taxa for biplot analysis based on their similarity scores.
#'
#' @param psIN A phyloseq object containing microbiome data.
#' @param project A string representing the project name, used for file naming.
#' @param TaxaLevel A string specifying the taxonomic level to be analyzed.
#' @param data_type A string indicating the type of data (e.g., 'dada2', 'Nephele').
#' @param bestPick The number of best taxa combinations to be considered.
#'
#' @details
#' The function aggregates data at the specified taxonomic level and uses bioenv stepwise selection to find the best fitting taxa combinations based on their similarity to environmental variables or other factors. The function is useful for preparing data for biplot analysis in microbiome studies.
#'
#' @return
#' Returns a vector of taxa names that represent the best combination for biplot analysis. The function also prints out similarity scores for each combination.
#'
#' @examples
#' # Example usage:
#' best_taxa <- Go_pickTaxa(psIN = my_phyloseq_object,
#'                          project = "MyMicrobiomeStudy",
#'                          TaxaLevel = "Genus",
#'                          data_type = "16S",
#'                          bestPick = 5)
#'
#' @export

Go_pickTaxa <- function(psIN, project, TaxaLevel, data_type, bestPick){
  # ------------- Aggregate
  cat("#--  Getting Taxa  --#\n")
  ranks <- c("Phylum", "Class", "Order", "Family", "Genus", "Species")
  if (data_type == "dada2" | data_type == "DADA2") {
    otu.filt <- as.data.frame(t(otu_table(psIN))) # for dada2 
  } else if (data_type == "Nephele" | data_type == "nephele" | data_type == "Other" | data_type == "other") {
    otu.filt <- as.data.frame((otu_table(psIN))) # not for dada2 
  }
  
  otu.filt$Genus <- getTaxonomy(otus=rownames(otu.filt), tax_tab=tax_table(psIN),taxRanks=ranks, level= TaxaLevel)
  agg <- aggregate(. ~ Genus, otu.filt, sum, na.action=na.pass) #add na.action=na.pass if have error "no rows to aggregate"
  genera <- agg$Genus
  agg <- agg[,-1]
  rownames(agg) <- genera
  agg_t <-t(agg)
  dim(agg)
  
  
  # ------------- biplot Genus 선택 
  

  #-----------------Parameters----------------#
  cmethod<-"spearman" 
  #Correlation method to use: pearson, spearman, kendall
  fmethod<-"bray" 
  #Fixed distance method: euclidean, manhattan, gower, altGower, canberra, bray, kulczynski, morisita,horn, binomial, and cao
  vmethod<-"bray" 
  #Variable distance method: euclidean, manhattan, gower, altGower, canberra, bray, kulczynski, morisita,horn, binomial, and cao
  nmethod<-"bray" 
  #NMDS distance method:  euclidean, manhattan, gower, altGower, canberra, bray, kulczynski, morisita,horn, binomial, and cao jaccard
  
  best5 <- bestPick #100
  
  if (data_type == "dada2" | data_type == "DADA2") {
    res.bv.step.biobio <- bv.step(wisconsin(agg_t), wisconsin(agg_t), 
                                  fix.dist.method=fmethod, var.dist.method=vmethod,correlation.method=cmethod,
                                  scale.fix=FALSE, scale.var=FALSE, 
                                  max.rho=0.95, min.delta.rho=0.001,
                                  random.selection=TRUE,
                                  prop.selected.var=0.3,
                                  num.restarts= best5,
                                  output.best=best5,
                                  var.always.include=NULL) 
  } else if (data_type == "Nephele" | data_type == "nephele" | data_type == "Other" | data_type == "other") {
    res.bv.step.biobio <- bv.step(wisconsin(agg), wisconsin(agg), 
                                  fix.dist.method=fmethod, var.dist.method=vmethod,correlation.method=cmethod,
                                  scale.fix=FALSE, scale.var=FALSE, 
                                  max.rho=0.95, min.delta.rho=0.001,
                                  random.selection=TRUE,
                                  prop.selected.var=0.3,
                                  num.restarts= best5,
                                  output.best=best5,
                                  var.always.include=NULL) 
  }
  
  taxaNames<-colnames(agg_t);taxaNames
  bestTaxaFit <- ""
  for(i in (1:length(res.bv.step.biobio$order.by.best$var.incl))){
    bestTaxaFit[i]<-paste(paste(taxaNames[as.numeric(unlist(strsplit(res.bv.step.biobio$order.by.best$var.incl[i], split=",")))],collapse=' + '), " = ",res.bv.step.biobio$order.by.best$rho[i],sep="")
  }
  
  bestTaxaFit <- data.frame(bestTaxaFit);bestTaxaFit 
  colnames(bestTaxaFit)<-"Best combination of taxa with similarity score"
  
  bio.keep <- as.numeric(unlist(strsplit(res.bv.step.biobio$order.by.best$var.incl[1], ",")));bio.keep
  
  taxaNames <- colnames(agg_t);taxaNames
  
  cat("\n Use this numbers for biplot;  \n")

  
  print(bio.keep)
  for (taxa in bio.keep){
    cat(sprintf("%s : %s \n", taxa, taxaNames[taxa]))
  }
  return(taxaNames)
}
