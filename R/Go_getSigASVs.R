
#' Extract Significant ASVs (Amplicon Sequence Variants) and Generate FASTA File
#'
#' @param project String representing the project name.
#' @param file_path Path to the directory containing the files.
#' @param files File extension filter (default is ".csv").
#' @param SigASVs File path of the significant ASVs data.
#' @param sampledata Data frame containing sample information.
#' @param target.column The target column to filter (e.g., "deseq2", "ancom").
#' @param target.group Groups to consider (e.g., "up", "down").
#' @param name Optional name for the output file.
#'
#' @details
#' This function filters ASVs based on significance from DESeq2 or ANCOM analysis. It creates FASTA files for each
#' significant ASV group (e.g., upregulated or downregulated). This is useful for further sequence-based analysis or
#' phylogenetic studies.
#'
#' @return
#' A data frame of significant ASVs with corresponding FASTA files in the specified output directory.
#'
#' @examples
#' Go_getSigASVs(project = "ProjectName",
#'              file_path = "/path/to/files",
#'              SigASVs = "/path/to/SigASVs.csv",
#'              sampledata = sample_dataframe,
#'              target.column = "deseq2",
#'              target.group = c("up", "down"))
#'
#' @export

Go_getSigASVs <- function(project, file_path, files= ".csv", 
                          SigASVs,
                          sampledata,
                          target.column,
                          target.group,
                          name=NULL){
  
  #===== output files
  out <- file.path(sprintf("%s_%s",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out)) dir.create(out)
  
  out_path <- file.path(sprintf("%s_%s/sigTab",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out_path)) dir.create(out_path)
  
  #out_csv <- file.path(sprintf("%s_%s/sigTab/csv",project, format(Sys.Date(), "%y%m%d"))) 
  #if(!file_test("-d", out_csv)) dir.create(out_csv)
  
  out_fasta <- file.path(sprintf("%s_%s/sigTab/fasta",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out_fasta)) dir.create(out_fasta)

  #===== check input
  if(!target.column %in% c("deseq2", "ancom")) {
    stop(paste("Invalid column name:", target.column, ". Expected 'deseq2' or 'ancom'."))
  }
    
  if (!all(target.group %in% c("up", "down")) || length(target.group) > 2) {
    stop(paste("Invalid group. Expected 'up', 'down', or both ('up', 'down')."))
  }
  

  #===== make the table
  sigTaxa <- Go_mergeTab(pattern=files, file_path = path);sigTaxa
  tabSig <- as.data.frame(subset(sigTaxa, sigTaxa[,target.column] %in% target.group)) 
  tabSig$ShortName <- gsub(" ","_", tabSig$ShortName)
  
  #===== get fasta header names
  
  for (i in 1:dim(tabSig)[1]) {
    if(tabSig[,target.column] [i] == "up"){
      tabSig$headers[i] <- sprintf(">%s|%s|%s|%s",tabSig[,target.column] [i],tabSig$mvar[i],tabSig$smvar[i],tabSig$ShortName[i])
    }else if(tabSig[,target.column] [i] == "down"){
      tabSig$headers[i] <- sprintf(">%s|%s|%s|%s",tabSig[,target.column] [i],tabSig$mvar[i],tabSig$basline[i],tabSig$ShortName[i])
    }
  }
  
  #===== subset by names
  bugs <- unique(tabSig$ShortName);bugs
  
  for(bug in bugs){
    tabSig.sel <- subset(tabSig, ShortName == bug)
    
    # Simplify all names
    tabSig.sel$headers <- gsub("(\\w)\\w+_(\\w+)", "\\1.\\2", tabSig.sel$headers)

    headers <- vector(dim(tabSig.sel)[1], mode="character")
    for (i in 1:dim(tabSig.sel)[1]) {
      tabSig.sel$headers[i] <- paste(tabSig.sel$headers[i], i, sep="")
    }
    
    
    df <- tabSig.sel %>% select(all_of(c("headers", target.column)))
    
    
    df$headers <- gsub(">","",df$headers)
    rownames(df) <- df$headers
    
    seqs.tab <- tabSig.sel$taxa
    seqs <- getSequences(seqs.tab)
    fasta <- c(rbind(tabSig.sel$headers, seqs))
    
    # write.csv(df,quote = FALSE,col.names = NA, file=sprintf("%s/%s.%s.map.csv",out_csv,bug,format(Sys.Date(), "%y%m%d"),sep="/"))
    write(fasta, file=sprintf("%s/%s.%s.seqs.fna", out_fasta,bug,format(Sys.Date(), "%y%m%d"),sep="/"))
  }

  
  tab <- read.csv(SigASVs,row.names=1,check.names=FALSE);dim(tab)

  tab.sel <- tab[unique(tabSig$taxa),] 

  write.csv(tab.sel,quote = FALSE,col.names = NA, file=sprintf("%s/%s.%s.%s%s.csv",out_path,project, "sigASVs",
                                                               ifelse(is.null(name), "", paste(name, ".", sep = "")), 
                                                               format(Sys.Date(), "%y%m%d"),sep="/"))

  return(tabSig)
}
