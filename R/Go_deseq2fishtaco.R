
#' Convert DESeq2 Results for FishTaco Analysis
#'
#' @param psIN The Phyloseq object containing the microbiome data.
#' @param project The name of the project, used for naming output files.
#' @param file_path Path to the folder containing DESeq2 output files.
#'
#' @details
#' This function processes DESeq2 results and prepares them for FishTaco analysis. It reads the DESeq2 output files, extracts relevant taxonomic information from the Phyloseq object, and aggregates data at the KEGG Orthology (KO) level. The function then filters and exports the processed data for each DESeq2 file, making it compatible with FishTaco input requirements.
#'
#' @return
#' For each DESeq2 output file, a corresponding CSV file is created in the specified output directory. These files are formatted for use with the FishTaco tool.
#'
#' @examples
#' Go_deseq2fishtaco(psIN = my_phyloseq_object,
#'                   project = "MyMicrobiomeProject",
#'                   file_path = "path/to/deseq2/results")
#'
#' @export

Go_deseq2fishtaco <- function(psIN, project, file_path){
  
  # out dir
  out <- file.path(sprintf("%s_%s",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out)) dir.create(out)
  out_path <- file.path(sprintf("%s_%s/pdf",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out_path)) dir.create(out_path)
  out_deseq2 <- file.path(sprintf("%s_%s/table/deseq2",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out_deseq2)) dir.create(out_deseq2)
  out_fishtaco <-file.path(sprintf("%s_%s/table/deseq2/funTabforFishtaco",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out_fishtaco)) dir.create(out_fishtaco)
  
  
  
  
  path <- file_path
  print(path)
  filenames <- list.files(path, pattern="Forfishtaco.csv")
  print(filenames)
  
  for(f in filenames){
    top.deseq2 <- read.csv(sprintf("%s/%s",path,f),row.names=1,check.names=FALSE)
    print(head(top.deseq2))
    ranks <-c("KO", "KO.des","Genus","Species","Path","Path.des")
    # otu table
    otu.filt <- as.data.frame(otu_table(psIN)) 
    
    otu.filt$ko <- getTaxonomy(otus=rownames(otu.filt), tax_tab=tax_table(psIN), taxRanks=ranks,  level="KO")
    agg <- aggregate(. ~ ko, otu.filt, sum)
    
    
    agg.clean <- agg[agg$ko %in% top.deseq2$x, ]# agg[!agg$ko %in% top.deseq2$x, ]  top.deseq2$x만 제거
    
    write.csv(agg.clean, quote = FALSE,col.names = NA,file=sprintf("%s/%s",out_fishtaco,f,sep="/"))
  }
}
