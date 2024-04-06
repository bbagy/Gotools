#' Calculate Sequence Lengths for Phyloseq Data
#'
#' This function calculates the sequence lengths from ASV sequences in a Phyloseq object.
#' It can operate without any filters or can filter sequences based on their length before
#' calculating and plotting the distribution of sequence lengths.
#'
#' @param psIN A Phyloseq object containing OTU/ASV counts and associated metadata.
#' @param from An optional numeric value specifying the minimum sequence length to include
#'             in the analysis. If NULL (default), no minimum length filtering is applied.
#' @param to An optional numeric value specifying the maximum sequence length to include
#'           in the analysis. If NULL (default), no maximum length filtering is applied.
#'
#' @return If no filter is applied, it prints and plots the distribution of sequence lengths.
#'         If filters are applied, it returns a new Phyloseq object containing only sequences
#'         within the specified length range.
#'
#' @examples
#' # Assuming 'ps' is a Phyloseq object
#' Go_SeqLengths(ps)
#' Go_SeqLengths(ps, from = 100, to = 200)
#'
#' @export
#' @importFrom phyloseq otu_table sample_data merge_phyloseq
#' @importFrom stats nchar table
#' @importFrom graphics hist


Go_SeqLengths <- function(psIN, from=NULL, to=NULL){
  # Ensure the phyloseq package is loaded
  if (!requireNamespace("phyloseq", quietly = TRUE)) {
    stop("The 'phyloseq' package is not installed. Please install it first.")
  }

  # Extract ASV sequences from the phyloseq object
  seqtab <- otu_table(psIN)

  # Ensure seqtab is a matrix and get sequences
  if (!is(seqtab, "matrix")) {
    stop("OTU table in the phyloseq object is not a matrix.")
  }
  asv_sequences <- colnames(seqtab)

  # Check if sequences contain only valid characters
  if(!all(sapply(strsplit(asv_sequences, ""), function(x) all(x %in% c("A", "T", "C", "G", "N"))))) {
    #stop("Sequences contain invalid characters.")
    asv_sequences <- rownames(seqtab)
  }
  # Define sequence lengths
  sequence_lengths <- nchar(asv_sequences)

  if(!is.null(from) && !is.null(to)) {
    # Filter sequences based on the provided range
    valid_indices <- which(sequence_lengths >= from & sequence_lengths <= to)
    tt <- try(seqtab <- seqtab[, valid_indices],T)
    if (class(tt) =="try-error"){
      seqtab.t <- t(seqtab)
      seqtab <- seqtab.t[, valid_indices]
    }

    asv_sequences <- asv_sequences[valid_indices]
    sequence_lengths <- sequence_lengths[valid_indices] # Update sequence lengths after filtering
  }

  # Create a table of sequence length distribution and plot it
  length_distribution <- table(sequence_lengths)
  print(length_distribution)

  # Adjust the main title of the histogram based on filtering
  main_title <- if (!is.null(from) && !is.null(to)) {
    paste("Distribution of Sequence Lengths (Filtered: ", from, "-", to, ")", sep="")
  } else {
    "Distribution of Sequence Lengths"
  }

  hist(sequence_lengths, breaks=100, main=main_title, xlab="Sequence Length", ylab="Frequency")

  # Create a new phyloseq object if filtered
  if(!is.null(from) && !is.null(to)) {
    new_ps <- phyloseq::phyloseq(phyloseq::otu_table(seqtab, taxa_are_rows = is_taxa_are_rows(psIN)), phyloseq::tax_table(tax_table(psIN)))

    # Attempt to merge sample_data if present
    tryCatch({
      if (!is.null(sample_data(psIN))) {
        new_ps <- phyloseq::merge_phyloseq(new_ps, sample_data(psIN))
      }
    }, error = function(e) {
      cat("No sample_data present in the phyloseq object.\n")
    })


    # Attempt to merge tree if present
    tryCatch({
      if (!is.null(phy_tree(psIN))) {
        new_ps <- phyloseq::merge_phyloseq(new_ps, phy_tree(psIN))
      }
    }, error = function(e) {
      cat("No phy_tree present in the phyloseq object.\n")
    })
    print(new_ps)

    return(new_ps)
  }
}
