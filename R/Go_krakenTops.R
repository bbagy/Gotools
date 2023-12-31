
#' Process Kraken MPA Output for Top Taxa and Prepare for RPKM Analysis
#'
#' This function processes Kraken MPA output, focusing on top taxa, and prepares data for RPKM analysis.
#'
#' @param project Name of the project or analysis.
#' @param mpa File path of the Kraken MPA output.
#'
#' @return The function processes the Kraken MPA data, extracts the top taxa, and saves the relevant
#'         information for further RPKM analysis. It saves a Phyloseq object and other required data as .csv and .rds files.
#'
#' @details
#' The function performs a series of steps to process Kraken MPA output. It reads the MPA file,
#' extracts top taxa information, formats the data, and saves it for downstream RPKM analysis.
#' It generates a Phyloseq object and exports relevant tables for species names and OTU.
#'
#' @examples
#' # Example usage:
#' Go_krakenTops(project = "MyProject",
#'               mpa = "path/to/kraken_mpa_output.txt")
#'
#' @export


Go_krakenTops <- function(project,
                          mpa){
  
  # out dir
  # RPKM_tem <- file.path(sprintf("%s","RPKM_tem"));if(!file_test("-d", RPKM_tem)) dir.create(RPKM_tem)
  RPKM_input <- file.path(sprintf("%s","RPKM_input"));if(!file_test("-d", RPKM_input)) dir.create(RPKM_input)
  rds <- file.path(sprintf("%s","2_rds"));if(!file_test("-d", rds)) dir.create(rds)
  
  # read kraken mpa
  mpatable <- read.table(mpa, header=T, as.is=T, sep="\t", comment.char="", row.names=1, quote="");head(mpatable) 
  
  #L1 <- subset(mpatable, grepl("k__Fungi", rownames(mpatable)))
  L1 <- subset(mpatable, grepl("d__Bacteria", rownames(mpatable)))
  L2 <- subset(L1, grepl("p__", rownames(L1)))
  L3 <- subset(L2, grepl("c__", rownames(L2)))
  L4 <- subset(L3, grepl("o__", rownames(L3)))
  L5 <- subset(L4, grepl("f__", rownames(L4)))
  L6 <- subset(L5, grepl("g__", rownames(L5)))
  L7 <- subset(L6, grepl("s__", rownames(L6)))
  
  data <- L7[, setdiff(1:ncol(L7), grep(".Bacterial.kraken.1", colnames(L7)))]
  
  tt <- try(colnames(data) <- gsub("_mpa", "", gsub("_out.txt", "", colnames(data))),T)

  if(class(tt) =="try-error"){
    colnames(L7) <- gsub("_mpa", "", gsub("_out.txt", "", colnames(L7)));colnames(L7)
    rownames(L7) <- gsub("\\|", ";", rownames(L7))
    colnames(L7) <-  gsub("X", "", colnames(L7));head(colnames(L7))
    data <- L7
  }else{
    colnames(data) <- gsub("_mpa", "", gsub("_out.txt", "", colnames(data)));colnames(data)
    rownames(data) <- gsub("\\|", ";", rownames(data))
    colnames(data) <-  gsub("X", "", colnames(data));head(colnames(data))
  }
  

  
  
  
  
  
  
  
  #ranklist <- c("Rank1", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  ranklist <- c("Rank1", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  taxlist <- lapply(rownames(data), function(x) parse_taxonomy_qiime(x))
  taxa <- matrix(NA, nrow=nrow(data), ncol=length(ranklist)); colnames(taxa) <- ranklist
  
  for (i in 1:length(taxlist)) {
    taxa[i, names(taxlist[[i]])] <- taxlist[[i]]
  }
  
  tt <- tax_table(taxa); rownames(tt) <- rownames(data)
  ps <- phyloseq(otu_table(data, taxa_are_rows=T), tt);ps
  
  
  #random_tree = rtree(ntaxa(ps), rooted=TRUE, tip.label=taxa_names(ps))
  #ps5 <- merge_phyloseq(ps, random_tree);ps5
  #ps.relative <- transform_sample_counts(ps, function(x) x / sum(x) )
  #ps.rarefied <- rarefy_even_depth(ps, sample.size=40224, rngseed=nsamples(ps))
  
  
  #------------------------------#
  #---       Kraken RPKM      ---#
  #------------------------------#
  # make a list species names
  taxaTab <- data.frame(taxa)
  taxaTab$Species <-  gsub(" ", "_", taxaTab$Species);head(taxaTab$Species)
  
  write.table(taxaTab$Species, quote=F, sep="\t", row.names=F, col.names=F,
              file=sprintf("%s/%s.speciesName.txt", RPKM_input, project))
  
  write.csv(data.frame(taxa), quote = FALSE,col.names = NA, #row.names = FALSE,
            file=sprintf("%s/%s.taxaTable_for_rpkm.%s.csv", RPKM_input, project, format(Sys.Date(), "%y%m%d"),sep="/"))
  
  # run GoGenomeSize.sh speciesName.txt genomeSize.txt  / 이름을 Gofriend로 할까나.
  
  # make a list species names
  otu <- otu_table(data, taxa_are_rows=T)
  rownames(otu) <- taxaTab$Species
  otu <- data.frame(otu)
  
  saveRDS(ps, sprintf("%s/ps.%s.%s.rds",rds, project, format(Sys.Date(), "%y%m%d")))
  
  write.csv(otu, quote = FALSE,col.names = NA, #row.names = FALSE,
            file=sprintf("%s/%s.speciesTab_for_rpkm.%s.csv", RPKM_input, project, format(Sys.Date(), "%y%m%d"),sep="/"))
  
  
}
