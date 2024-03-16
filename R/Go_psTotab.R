
#' Convert Phyloseq Object to Table and Export Sequences
#'
#' @param psIN Input phyloseq object.
#' @param project Name of the project for naming output files.
#'
#' @details
#' This function is designed to extract data from a phyloseq object and convert it into a more accessible table format. It also extracts sequence data and saves it as a FASTA file. This is particularly useful for further analysis or data sharing.
#'
#' @return
#' The function generates two files:
#' 1. A CSV file containing the abundance data and taxonomy information from the phyloseq object.
#' 2. A FASTA file containing sequence data.
#'
#' Both files are named using the project name and are saved in a specified output directory.
#'
#' @examples
#' Go_psTotab(psIN = my_phyloseq_object,
#'            project = "MyMicrobiomeProject")
#'
#' @export

Go_psTotab <- function(psIN, project){
  out <- file.path("1_out")
  if(!file_test("-d", out)) dir.create(out)

  #====== step 1 read ps object
  seqtab.nochim <- as.matrix(otu_table(psIN))
  tax <- as.matrix(tax_table(psIN))


  #====== step 2 extract fna
  seqs <- getSequences(t(seqtab.nochim))
  headers <- paste(">", seqs, sep="")
  fasta <- c(rbind(headers, seqs))


  if (nchar(headers[1]) < 100){
    seqs <- getSequences(seqtab.nochim)
    headers <- paste(">", seqs, sep="")
    fasta <- c(rbind(headers, seqs))
    write(fasta, file=sprintf("%s/%s.%s.psTotab.seqs.fna",out, project, format(Sys.Date(), "%y%m%d"),sep="/"))
  }else{
    write(fasta, file=sprintf("%s/%s.%s.psTotab.seqs.fna",out, project, format(Sys.Date(), "%y%m%d"),sep="/"))
  }



  #====== step 3 get the table
  otu <- as.data.frame(t(otu_table(psIN)));dim(otu)
  tax <- tax_table(psIN);dim(tax)


  tt <- try( otuTable <- cbind(otu,tax),T)
  if(class(tt) == "try-error"){
    otu <- as.data.frame(otu_table(psIN));dim(otu)
    tax <- tax_table(psIN);dim(tax)
    otuTable <- cbind(otu,tax)

    seqs <- getSequences(t(seqtab.nochim))
    headers <- paste(">", seqs, sep="")
    fasta <- c(rbind(headers, seqs))

    write(fasta, file=sprintf("%s/%s.%s.psTotab.seqs.fna",out, project, format(Sys.Date(), "%y%m%d"),sep="/"))

  }else{
    otuTable <- cbind(otu,tax)
  }

  write.csv(tax, quote = FALSE,col.names = NA,#row.names = FALSE,
            file=sprintf("%s/%s.%s.psTotab.tax.csv",out,project,format(Sys.Date(), "%y%m%d"), sep="/"))

  write.csv(otu, quote = FALSE,col.names = NA,#row.names = FALSE,
            file=sprintf("%s/%s.%s.psTotab.asv.csv",out,project,format(Sys.Date(), "%y%m%d"), sep="/"))

  write.csv(otuTable, quote = FALSE,col.names = NA,#row.names = FALSE,
            file=sprintf("%s/%s.%s.psTotab.asvTable.csv",out,project,format(Sys.Date(), "%y%m%d"), sep="/"))

}



