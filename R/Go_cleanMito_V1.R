
#' Remove Mitochondrial Reads from Phyloseq Object
#'
#' @param psIN Phyloseq object containing microbiome data.
#' @param project Name of the project, used for output file naming.
#'
#' @details
#' This function processes a phyloseq object to remove mitochondrial reads. It specifically targets reads assigned to the order "Rickettsiales", often considered as mitochondrial in microbiome datasets. The function also removes taxa with NA in the Phylum level and filters out taxa with very low counts.
#'
#' @return
#' The function generates a phyloseq object with mitochondrial reads removed, a corresponding FASTA file of sequences, and a CSV file containing the processed ASV table. These files are saved in specified directories.
#'
#' @examples
#' Go_cleanMito(psIN = my_phyloseq_object,
#'              project = "MyMicrobiomeProject")
#'
#' @export

Go_cleanMito <- function(psIN, project){
  out <- file.path("1_out")
  if(!file_test("-d", out)) dir.create(out)
  
  rds <- file.path("2_rds") 
  if(!file_test("-d", rds)) dir.create(rds)
  
  #====== step 1 removing mitochondria reads
  seqtab.nochim <- as.matrix(otu_table(psIN))
  tax <- as.matrix(tax_table(psIN))

  # is.mito <- tax[,"Order"] %in% "Rickettsiales" 
  is.a <- tax[,"Order"] %in% "Rickettsiales" 
  seqtab.a <- seqtab.nochim[,!is.a];dim(seqtab.a)
  tax.a <- tax[!is.a,]
  
  is.NA <- tax.a[,"Phylum"] %in% NA 
  seqtab.noNA <- seqtab.a[,!is.NA] ;dim(seqtab.noNA)
  tax.noNA <- tax.a[!is.NA,]
  # remove na column sum
  seqtab.noNA <-data.frame(t(seqtab.noNA))
  seqtab.noNA.sum <- seqtab.noNA[,colSums(seqtab.noNA) > 1];dim(seqtab.noNA.sum)
  seqtab.Nomito <- as.matrix(t(seqtab.noNA.sum))
  
  seqs <- getSequences(seqtab.Nomito)
  headers <- paste(">", seqs, sep="")
  fasta <- c(rbind(headers, seqs))
  
  write(fasta, file=sprintf("%s/%s.%s.No_mito.seqs.fna",out, project, format(Sys.Date(), "%y%m%d"),sep="/"))
  
  
  #====== step 2 merge phyloseq
  ps.Nomito <- phyloseq(otu_table(seqtab.Nomito, taxa_are_rows=FALSE), tax_table(tax.noNA));ps.Nomito
  sample_names(ps.Nomito)
  sample_names(ps.Nomito) <- gsub("X","",sample_names(ps.Nomito));sample_names(ps.Nomito)

  #sample_names(ps.Nomito) <- gsub("\\_.*","",sample_names(ps.Nomito));sample_names(ps.Nomito)
  
  
  
  
  #====== step 3 get the table
  otu <- as.data.frame(t(otu_table(ps.Nomito)));dim(otu)
  tax <- tax_table(ps.Nomito);dim(tax)
  
  otuTable <- cbind(otu,tax)
  
  write.csv(otuTable, quote = FALSE,col.names = NA,#row.names = FALSE, 
            file=sprintf("%s/%s.%s.asvTable_No_mito.csv",out,project,format(Sys.Date(), "%y%m%d"), sep="/"))
  
  saveRDS(ps.Nomito, sprintf("%s/ps_No_mito.%s.%s.rds", rds, project,format(Sys.Date(), "%y%m%d")))
  
  
  print(psIN)
  print(ps.Nomito)  
}
