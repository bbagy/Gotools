
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
  if(!dir.exists(out)) dir.create(out)
  date_tag <- format(Sys.Date(), "%y%m%d")
  seq_path <- file.path(out, sprintf("%s.%s.psTotab.seqs.fna", project, date_tag))
  tax_path <- file.path(out, sprintf("%s.%s.psTotab.tax.csv", project, date_tag))
  otu_path <- file.path(out, sprintf("%s.%s.psTotab.asv.csv", project, date_tag))
  otu_table_path <- file.path(out, sprintf("%s.%s.psTotab.asvTable.csv", project, date_tag))

  #====== step 1 read ps object
  seqtab.nochim <- as.matrix(otu_table(psIN))
  tax <- as.matrix(tax_table(psIN))


  #====== step 2 extract fna
  fasta_written <- FALSE
  taxa_ids <- taxa_names(psIN)
  looks_like_dna <- length(taxa_ids) > 0 && all(grepl("^[ACGTN]+$", taxa_ids))

  if (looks_like_dna) {
    seqs <- dada2::getSequences(t(seqtab.nochim))
    headers <- paste(">", seqs, sep="")
    fasta <- c(rbind(headers, seqs))

    if (length(headers) > 0 && nchar(headers[1]) < 100){
      seqs <- dada2::getSequences(seqtab.nochim)
      headers <- paste(">", seqs, sep="")
      fasta <- c(rbind(headers, seqs))
    }

    write(fasta, file = seq_path)
    fasta_written <- TRUE
  } else {
    message("Go_psTotab(): taxa names do not look like DNA sequences; skipping FASTA export.")
  }



  #====== step 3 get the table
  otu <- as.data.frame(t(otu_table(psIN)));dim(otu)
  tax <- tax_table(psIN);dim(tax)


  tt <- try( otuTable <- cbind(otu,tax),T)
  if(inherits(tt, "try-error")){
    otu <- as.data.frame(otu_table(psIN));dim(otu)
    tax <- tax_table(psIN);dim(tax)
    otuTable <- cbind(otu,tax)
    if (looks_like_dna && !fasta_written) {
      seqs <- dada2::getSequences(t(seqtab.nochim))
      headers <- paste(">", seqs, sep="")
      fasta <- c(rbind(headers, seqs))
      write(fasta, file = seq_path)
    }

  }else{
    otuTable <- cbind(otu,tax)
  }

  write.csv(tax, quote = TRUE, file = tax_path)

  write.csv(otu, quote = TRUE, file = otu_path)

  write.csv(otuTable, quote = TRUE, file = otu_table_path)

}
