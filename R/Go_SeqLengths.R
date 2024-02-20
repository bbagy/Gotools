#' Calculate Sequence Lengths for Phyloseq Data
#'
#' This function calculates the sequence lengths from ASV sequences in a Phyloseq object.
#' It can operate without any filters or apply a minimum sequence length filter before
#' calculating and plotting the distribution of sequence lengths.
#'
#' @param psIN A Phyloseq object containing OTU/ASV counts and associated metadata.
#' @param filter An optional numeric value specifying the minimum sequence length to include
#'               in the analysis. If NULL (default), no filtering is applied.
#'
#' @return If no filter is applied, it prints and plots the distribution of sequence lengths.
#'         If a filter is applied, it returns a new Phyloseq object containing only sequences
#'         with lengths greater than or equal to the filter value.
#'
#' @examples
#' # Assuming 'ps' is a Phyloseq object
#' Go_SeqLengths(ps)
#' Go_SeqLengths(ps, filter = 100)
#'
#' @export
#' @importFrom phyloseq otu_table sample_data merge_phyloseq
#' @importFrom stats nchar table
#' @importFrom graphics hist
Go_SeqLengths <- function(psIN,
                          filter=NULL){

  # Extract ASV sequences from the phyloseq object
  seqtab <- otu_table(psIN)

  asv_sequences <- colnames(seqtab)

  if(all(strsplit(asv_sequences, "")[[1]] %in% c("A", "T", "C", "G", "N"))) {
  } else {
    seqtab <- t(otu_table(psIN))
    asv_sequences <- colnames(seqtab)
  }




  if(is.null(filter)){
    # Calculate the length of each sequence
    sequence_lengths <- nchar(asv_sequences)
    # Create a table of sequence length distribution
    length_distribution <- table(sequence_lengths)

    # table(nchar(getSequences(seqtab)))

    # Print the distribution
    print(length_distribution)
    # Plot a histogram of sequence lengths
    hist(sequence_lengths, breaks=100, main="Distribution of Sequence Lengths", xlab="Sequence Length", ylab="Frequency")
  }else{
    seqtab2 <- seqtab[,nchar(colnames(seqtab)) >= filter]
                      # Calculate the length of each sequence
    sequence_lengths <- nchar(colnames(seqtab2))
    # Create a table of sequence length distribution
    length_distribution <- table(sequence_lengths)
    hist(sequence_lengths, breaks=100, main="Distribution of Sequence Lengths (filtered)", xlab="Sequence Length", ylab="Frequency")


    new_ps1 <- phyloseq(seqtab2, tax_table(psIN))

    if (!is.null(sample_data(psIN))) {
      new_ps1 <- merge_phyloseq(new_ps1, sample_data(psIN))
    }

    if (!is.null(phy_tree(psIN))) {
      new_ps1 <- merge_phyloseq(new_ps1, phy_tree(psIN))
    }

  return(new_ps1)
  }


}
