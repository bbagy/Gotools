
#' Cluster ASVs to Form OTUs
#'
#' This function clusters ASVs (Amplicon Sequence Variants) to form OTUs (Operational Taxonomic Units) based on a specified similarity percentage. It then reclassifies the clustered sequences and creates a new phyloseq object.
#'
#' @param psIN A phyloseq object containing ASV data.
#' @param project A string representing the name of the project.
#' @param db A string indicating the database file path for taxonomy assignment.
#' @param percent A numeric value representing the similarity percentage for clustering.
#'
#' @return
#' A new phyloseq object with clustered OTUs and reclassified taxonomy.
#'
#' @details
#' The function performs sequence alignment, calculates distance matrices, and identifies clusters of ASVs using the DECIPHER package. It then merges ASV counts within the same cluster and reassigns taxonomy using the specified database. The result is a new phyloseq object with updated OTU counts and taxonomy.
#'
#' @examples
#' # Example usage:
#' ps_clustered <- Go_cluster(psIN = ps, project = "MyProject", db = "path/to/db", percent = 97)
#'
#' @export


Go_cluster <- function(psIN, project,db, percent){
  library(tibble)
  library(dplyr)
  # install.packages("remotes")
  # remotes::install_github("mikemc/speedyseq")
  library(speedyseq)
  # Packages that are required but not loaded:
  # library(DECIPHER)
  # library(Biostrings)
  
  # out dir
  rds <- file.path("2_rds") 
  if(!file_test("-d", rds)) dir.create(rds)
  
  out <- file.path("1_out") 
  if(!file_test("-d", out)) dir.create(out)
  
  # ----- Input ------#
  project.name <-sprintf("%s_%s",project,percent);project.name
  
  ps <- psIN
  x <- 1-percent/100;x
  
  seqtab <- otu_table(psIN)
  
  
  nproc <- 4 # set to number of cpus/processors to use for the clustering
  
  asv_sequences <- colnames(seqtab);head(asv_sequences)
  sample_names <- rownames(seqtab);(sample_names)
  dna <- Biostrings::DNAStringSet(asv_sequences)
  
  
  ## Find clusters of ASVs to form the new OTUs
  aln <- DECIPHER::AlignSeqs(dna, processors = nproc)
  d <- DECIPHER::DistanceMatrix(aln, processors = nproc)
  clusters <- DECIPHER::IdClusters(
    d, 
    method = "complete",
    cutoff = x , # use `cutoff = 0.03` for a 97% OTU 
    processors = nproc)
  
  ## Use dplyr to merge the columns of the seqtab matrix for ASVs in the same OTU
  # prep by adding sequences to the `clusters` data frame
  cluster <- clusters %>%
    add_column(sequence = asv_sequences)
  
  merged_seqtab <- seqtab %>% 
    t %>%
    rowsum(clusters$cluster) %>%
    t
  
  
  # rebuilt ASVs table 
  clustered <- distinct(cluster, cluster, .keep_all = TRUE)

  merged_seqtab.t <- data.frame(t(merged_seqtab))
  merged_seqtab.t$seqs <- factor(clustered$sequence[match(as.factor(rownames(merged_seqtab.t)), as.factor(clustered$cluster))])
  
  rownames(merged_seqtab.t) <- merged_seqtab.t$seqs
  merged_seqtab.t$seqs <- NULL
  
  seqtab <- as.matrix(t(merged_seqtab.t))
  
  #----- save seqs.fna for tree  -----#
  seqs <- getSequences(seqtab)
  headers <- paste(">", seqs, sep="")
  fasta <- c(rbind(headers, seqs))
  write(fasta, file=sprintf("%s/ASVs%s.%s.%s.seqs.fna",  out, percent,project, format(Sys.Date(), "%y%m%d"),sep="/"))
  

  #----- re classification  -----#
  tax <- assignTaxonomy(seqtab, db, 
                        taxLevels = c("Kingdom","Phylum", "Class", "Order", "Family", "Genus","Species"), 
                        minBoot = 80, verbose = TRUE, multithread = TRUE)
  
  
  #----- merge  -----# 
  ps_clustered <- phyloseq(otu_table(seqtab, taxa_are_rows=FALSE), tax_table(tax));ps_clustered
  
  #----- setting fix species names  -----# 
  
  tax <- data.frame(tax_table(ps_clustered))
  tax$Species.1 <- paste(tax$Genus,tax$Species)
  
  tax$Species <- NULL
  colnames(tax) <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
  tax_table(ps_clustered) <- as.matrix(tax)
  
  
  
  
  cat("#--  Before clustered  --#\n")
  print(psIN)
  cat("\n")
  cat(sprintf("#--  After clustered by %s   --#\n",percent))
  print(ps_clustered)
  
  saveRDS(ps_clustered, sprintf("%s/ps%s_filtered.%s.%s.rds", rds, percent, project.name,format(Sys.Date(), "%y%m%d")))
  
  otu <- as.data.frame(t(otu_table(ps_clustered)));dim(otu)
  tax <- tax_table(ps_clustered);dim(tax)
  otuTable <- cbind(otu,tax)
  
  write.csv(otuTable, quote = FALSE,col.names = NA,row.names = T, file=sprintf("%s/ASVs%s_clustered.%s.%s.csv", out, percent, project.name,format(Sys.Date(), "%y%m%d"), sep="/"))
  
  
}



