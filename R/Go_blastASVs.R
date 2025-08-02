#' Run BLAST on ASV Table Sequences
#'
#' This function takes a CSV file containing ASV sequences, converts them to FASTA format,
#' runs BLAST searches against a specified database, parses the results, and returns
#' a CSV file with added BLAST hit information.
#'
#' #===== download DB in local
#' cd DB
#' Download the nt db (300~500GB)
#' update_blastdb.pl nt
#' Download the RefSeq bacterial database
#' update_blastdb.pl --decompress refseq_rna
#' Or download the 16S microbial database
#' update_blastdb.pl --decompress 16S_ribosomal_RNA
#'
#' @param asvsTable The path to the CSV file containing ASV sequences with sequence IDs as row names.
#' @param blastDB The path to the BLAST database against which the sequences will be queried.
#'                 For example, "/Users/username/DB/blastDB/16S_ribosomal_RNA".
#'
#' @return Writes a CSV file to "1_out/updated_sequences_with_blast_results.csv"
#'         with the original ASV data and added columns for BLAST hits.
#'         The function also prints a message upon successful completion.
#'
#' @details The function processes each sequence individually, creates a temporary FASTA file,
#'          runs a BLAST search, and parses the output to extract the hit number and hit definition.
#'          It assumes that the BLAST tool is properly installed and accessible in the system's PATH.
#'          The results include the genus and species extracted from the hit definitions.
#'
#' @examples
#' Go_blastASVs(asvsTable="path/to/your/asvTable.csv",
#'              blastDB="/path/to/your/blastDB/16S_ribosomal_RNA")
#'
#' @export

Go_blastASVs <- function(project,
                         asvsTable,
                         blastDB) {
  # Load necessary libraries
  if (!requireNamespace("Biostrings", quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager")
    }
    BiocManager::install("Biostrings", ask = FALSE)
  }
  library(Biostrings, quietly = TRUE)

  # Set up the PATH environment variable for BLAST
  current_path <- Sys.getenv("PATH")
  new_path <- paste(current_path, "/Users/heekukpark/miniconda3/envs/blast/bin/", sep=":")# which blastn
  Sys.setenv(PATH = new_path)

  # Verify BLAST installation
  if (system("blastn -version", intern = F) == 0) {
    cat("BLAST is correctly installed and found in PATH.\n")
  } else {
    stop("BLAST not found. Please check your BLAST installation and PATH.")
  }

  # Read sequences
  sequences_df <- read.csv(asvsTable, row.names = 1, check.names = FALSE)
  dna_strings <- DNAStringSet(rownames(sequences_df))
  names(dna_strings) <- paste("Sequence", seq_along(dna_strings), sep="_")

  # Setup output directory
  output_dir <- "1_out/blast"
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # Run BLAST for each sequence
  blast_results <- sapply(names(dna_strings), function(seq_name) {
    fasta_file <- file.path(output_dir, paste0(seq_name, ".fasta"))
    output_file <- file.path(output_dir, paste0(seq_name, "_blast.xml"))
    writeXStringSet(dna_strings[seq_name], fasta_file)
    system2("blastn", args = c("-query", fasta_file, "-db", blastDB, "-out", output_file, "-outfmt", "5"))
    file.remove(fasta_file) # Clean up FASTA file immediately after use
    return(output_file)
  })

  sapply(names(dna_strings), function(seq_name) {
    fasta_file_path <- file.path(output_dir, paste0(seq_name, ".fasta"))
    file.remove(fasta_file_path)
  })

  # Check results
  print(blast_results)


  library(XML)

  # Parse BLAST results
  parse_blast_results <- function(file) {
    # Load and parse the XML file
    doc <- try(xmlParse(file), silent = TRUE)
    if (inherits(doc, "try-error")) {
      cat("Failed to parse XML file:", file, "\n")
      return(list(Hit_Number = NA, Hit_Definition = NA))
    }

    # Get iterations
    iterations <- getNodeSet(doc, "//BlastOutput_iterations/Iteration")
    if (length(iterations) == 0) {
      cat("No iterations found in file:", file, "\n")
      return(list(Hit_Number = NA, Hit_Definition = NA))
    }

    # Check for hits
    hits <- getNodeSet(iterations[[1]], "Iteration_hits/Hit")
    if (length(hits) == 0) {
      cat("No hits found in iteration of file:", file, "\n")
      return(list(Hit_Number = NA, Hit_Definition = NA))
    }

    # Extract the first hit's number and definition
    first_hit <- hits[[1]]
    hit_num <- xmlValue(getNodeSet(first_hit, "Hit_num")[[1]])
    hit_def <- xmlValue(getNodeSet(first_hit, "Hit_def")[[1]])

    return(list(Hit_Number = hit_num, Hit_Definition = hit_def))
  }


  # Apply the function to each result file
  parsed_results <- lapply(blast_results, parse_blast_results)

  # Convert results to a data frame
  results_df <- as.data.frame(do.call(rbind, parsed_results))

  # If sequences_df has the same order and number of sequences as blast_results
  #sequences_df$Hit_Number <- results_df$Hit_Number
  sequences_df$Hit_Definition <- results_df$Hit_Definition

  # Assuming Hit_Number and Hit_Definition might still be lists, convert them to character vectors
  #sequences_df$Hit_Number <- sapply(sequences_df$Hit_Number, function(x) ifelse(is.null(x), NA, x))
  sequences_df$Hit_Definition <- sapply(sequences_df$Hit_Definition, function(x) ifelse(is.null(x), NA, x))
  # Split each string by spaces and select the first two words
  sequences_df$Hit_Definition <- sapply(strsplit(sequences_df$Hit_Definition, " "), function(words) {
    paste(words[1:2], collapse=" ")
  })

  # Print the modified column to verify the changes
  print(sequences_df$Hit_Definition)


  # Write results to CSV
  write.csv(sequences_df, sprintf("1_out/%s.updated_sequences_with_blast_results.%s.csv",project, format(Sys.Date(), "%y%m%d")), row.names = TRUE)
  cat("Results have been successfully saved to 1_out/updated_sequences_with_blast_results.csv\n")
}




