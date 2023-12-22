
#' Generate Color Palette for Taxonomic Data
#'
#' @param cdf Data frame containing taxonomic information (phylum level) and corresponding colors.
#' @param taxaName Vector of taxonomic names for which colors are assigned.
#'
#' @details
#' This function creates a color palette for taxonomic data, primarily at the phylum level. It assigns distinct colors to
#' each phylum and ensures visual differentiation between various taxonomic groups.
#'
#' @return
#' A list containing two elements:
#'   - `$color_table`: A data frame mapping each taxa to its corresponding HSV color code.
#'   - `$coloring`: A named vector where names are taxa and values are color codes.
#'
#' @examples
#' color_data <- Go_color(cdf = taxonomic_dataframe, taxaName = vector_of_taxa_names)
#' # Access the color table
#' color_table <- color_data$color_table
#' # Access the coloring vector
#' coloring_vector <- color_data$coloring
#'
#' @export

Go_contam <- function(psIN, project,
                      name=NULL,
                      contaminant=NULL,
                      neg=NULL,
                      taxa=NULL) {


  out <- file.path("1_out") 
  if(!file_test("-d", out)) dir.create(out)
  
  rds <- file.path("2_rds") 
  if(!file_test("-d", rds)) dir.create(rds)

  
  map <- data.frame(sample_data(psIN)) 

  #====== step 1 define contaminated ps and ucontaminated ps
  if(!is.null(neg)){

    contams <- psIN
    map.contams <-  subset(map, map[,contaminant] %in% neg)
    sample_data(contams) <- map.contams

    samples <- psIN
    
    map.samples <-  subset(map, !(map[,contaminant] %in% neg))
    sample_data(samples) <- map.samples

    
    #contams <- subset_samples(psIN, map[,contaminant] == neg);contams
    #samples <- subset_samples(psIN, map[,contaminant] != neg);samples
    #====== step 2 substrate ucontaminated ps -  contaminated ps
    contams.n1 <- prune_taxa(taxa_sums(contams)>=1, contams) 
    samples.n1 <- prune_taxa(taxa_sums(samples)>=1, samples) 
    allTaxa <- names(sort(taxa_sums(psIN),TRUE))
    negtaxa <- names(sort(taxa_sums(contams.n1),TRUE))
    taxa.noneg <- allTaxa[!(allTaxa %in% negtaxa)]
    ps.decontam <- prune_taxa(taxa.noneg,samples.n1)
    print("step2")
  }else{
    ps.decontam <- psIN
  }



  if(!is.null(taxa)){
    seqtab.nochim <- as.matrix(otu_table(ps.decontam))
    tax <- as.matrix(tax_table(ps.decontam))
    is.a <- tax[,"Species"] %in% taxa 
    
    
    tt <- try(seqtab.a <- seqtab.nochim[,!is.a],T)
    if(class(tt) == "try-error"){
      seqtab.nochim <- as.matrix(t(otu_table(ps.decontam)))
      tax <- as.matrix(tax_table(ps.decontam))
      is.a <- tax[,"Species"] %in% taxa 
      seqtab.a <- seqtab.nochim[,!is.a];dim(seqtab.a)
    }else{
      seqtab.a <- seqtab.nochim[,!is.a];dim(seqtab.a)
    }
    
    tax.a <- tax[!is.a,]
    ps.decontam <- phyloseq(otu_table(seqtab.a, taxa_are_rows=FALSE), tax_table(tax.a));ps.decontam
    sample_names(ps.decontam)
    sample_names(ps.decontam) <- gsub("X","",sample_names(ps.decontam));sample_names(ps.decontam)
  }else{
    ps.decontam <- ps.decontam
  }


  #====== step 3 get decontam.fna
  seqs <- getSequences(otu_table(ps.decontam))
  headers <- paste(">", seqs, sep="")
  fasta <- c(rbind(headers, seqs))
  print("step3")
  write(fasta, file=sprintf("%s/%s%s.%s.decontam.fna",out, 
                            ifelse(is.null(name), "", paste(name, ".", sep = "")), 
                            project, 
                            format(Sys.Date(), "%y%m%d"),sep="/"))
  
  #====== step 4 get the table
  otu <- as.data.frame(t(otu_table(ps.decontam)));dim(otu)
  tax <- tax_table(ps.decontam);dim(tax)
  
  otuTable <- cbind(otu,tax)
  
  write.csv(otuTable, quote = FALSE,col.names = NA,#row.names = FALSE, 
            file=sprintf("%s/%s%s.%s.asvTable_decontam.csv",out,
                         ifelse(is.null(name), "", paste(name, ".", sep = "")), 
                         project,format(Sys.Date(), "%y%m%d"), sep="/"))
  
  saveRDS(ps.decontam, sprintf("%s/ps_decontam.%s%s.%s.rds", rds, 
                               ifelse(is.null(name), "", paste(name, ".", sep = "")), 
                               project,format(Sys.Date(), "%y%m%d")))
  
  
  print(psIN)
  print(ps.decontam)  
  
  return(ps.decontam)
}
