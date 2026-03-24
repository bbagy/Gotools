#' Run BLAST on ASV Table Sequences
#'
#' This function takes a CSV file containing ASV sequences, converts them to FASTA format,
#' runs BLAST searches against a specified database, parses the results, and writes:
#' 1. the original BLAST-annotated table
#' 2. a final ASV table with blast-guided species correction
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
#' @return Writes two CSV files to `1_out/`:
#' 1. `project.updated_sequences_with_blast_results.YYMMDD.csv`
#' 2. `project.final_asvTable.YYMMDD.csv`
#' If a matching `project.final_asvTable.*.csv` already exists in `1_out/`,
#' the function stops without running BLAST again.
#'
#' @details The final ASV table keeps existing DADA2 species names when they are already
#'          non-missing genus+species assignments. When species is effectively `Genus NA`,
#'          the function uses the BLAST species only if the BLAST genus matches the existing genus.
#'          Missing taxonomic ranks from phylum through genus are filled with the most recent
#'          available higher rank so the taxonomy remains report-friendly.
#'
#' @examples
#' Go_blastASVs(project = "Example",
#'              asvsTable = "path/to/your/asvTable.csv",
#'              blastDB = "/path/to/your/blastDB/16S_ribosomal_RNA")
#'
#' @export

Go_blastASVs <- function(project,
                         asvsTable,
                         blastDB) {
  existing_final_outputs <- Sys.glob(sprintf("1_out/%s.final_asvTable.*.csv", project))
  if (length(existing_final_outputs) > 0) {
    message(
      sprintf(
        "A final ASV table already exists for project '%s'. Reusing existing file(s): %s",
        project,
        paste(existing_final_outputs, collapse = ", ")
      )
    )
    return(invisible(existing_final_outputs))
  }

  if (!requireNamespace("Biostrings", quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager")
    }
    BiocManager::install("Biostrings", ask = FALSE)
  }
  if (!requireNamespace("XML", quietly = TRUE)) {
    install.packages("XML")
  }

  library(Biostrings, quietly = TRUE)
  library(XML, quietly = TRUE)

  current_path <- Sys.getenv("PATH")
  new_path <- paste(current_path, "/Users/heekukpark/miniconda3/envs/blast/bin/", sep = ":")
  Sys.setenv(PATH = new_path)

  if (system("blastn -version", intern = FALSE) == 0) {
    cat("BLAST is correctly installed and found in PATH.\n")
  } else {
    stop("BLAST not found. Please check your BLAST installation and PATH.")
  }

  sequences_df <- read.csv(asvsTable, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
  dna_strings <- DNAStringSet(rownames(sequences_df))
  names(dna_strings) <- paste("Sequence", seq_along(dna_strings), sep = "_")

  output_dir <- "1_out/blast"
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  is_valid_blast_xml <- function(file) {
    if (!file.exists(file)) {
      return(FALSE)
    }
    file_info <- file.info(file)
    if (is.na(file_info$size) || file_info$size <= 0) {
      return(FALSE)
    }
    doc <- try(XML::xmlParse(file), silent = TRUE)
    if (inherits(doc, "try-error")) {
      return(FALSE)
    }
    iterations <- XML::getNodeSet(doc, "//BlastOutput_iterations/Iteration")
    length(iterations) > 0
  }

  blast_results <- sapply(names(dna_strings), function(seq_name) {
    fasta_file <- file.path(output_dir, paste0(seq_name, ".fasta"))
    output_file <- file.path(output_dir, paste0(seq_name, "_blast.xml"))

    if (is_valid_blast_xml(output_file)) {
      return(output_file)
    }

    writeXStringSet(dna_strings[seq_name], fasta_file)
    system2("blastn", args = c("-query", fasta_file, "-db", blastDB, "-out", output_file, "-outfmt", "5"))
    file.remove(fasta_file)

    if (!is_valid_blast_xml(output_file)) {
      file.remove(output_file)
      writeXStringSet(dna_strings[seq_name], fasta_file)
      system2("blastn", args = c("-query", fasta_file, "-db", blastDB, "-out", output_file, "-outfmt", "5"))
      file.remove(fasta_file)
    }

    output_file
  })

  cat(length(blast_results), "BLAST XML files prepared in", output_dir, "\n")

  parse_blast_results <- function(file) {
    doc <- try(xmlParse(file), silent = TRUE)
    if (inherits(doc, "try-error")) {
      cat("Failed to parse XML file:", file, "\n")
      return(list(Hit_Number = NA_character_, Hit_Definition = NA_character_))
    }

    iterations <- getNodeSet(doc, "//BlastOutput_iterations/Iteration")
    if (length(iterations) == 0) {
      cat("No iterations found in file:", file, "\n")
      return(list(Hit_Number = NA_character_, Hit_Definition = NA_character_))
    }

    hits <- getNodeSet(iterations[[1]], "Iteration_hits/Hit")
    if (length(hits) == 0) {
      cat("No hits found in iteration of file:", file, "\n")
      return(list(Hit_Number = NA_character_, Hit_Definition = NA_character_))
    }

    first_hit <- hits[[1]]
    hit_num <- xmlValue(getNodeSet(first_hit, "Hit_num")[[1]])
    hit_def <- xmlValue(getNodeSet(first_hit, "Hit_def")[[1]])

    list(Hit_Number = hit_num, Hit_Definition = hit_def)
  }

  parsed_results <- lapply(blast_results, parse_blast_results)
  results_df <- as.data.frame(do.call(rbind, parsed_results), stringsAsFactors = FALSE)

  sequences_df$Hit_Number <- as.character(results_df$Hit_Number)
  sequences_df$Hit_Definition <- as.character(results_df$Hit_Definition)
  sequences_df$Hit_Definition <- sapply(base::strsplit(as.character(sequences_df$Hit_Definition), " "), function(words) {
    words <- words[!is.na(words) & nzchar(words)]
    if (length(words) >= 2) {
      paste(words[1:2], collapse = " ")
    } else if (length(words) == 1) {
      words[1]
    } else {
      NA_character_
    }
  })

  is_missing_taxon <- function(x) {
    x <- trimws(as.character(x))
    is.na(x) | x == "" | toupper(x) %in% c("NA", "N/A", "NULL")
  }

  find_tax_col <- function(nms, target) {
    hits <- nms[tolower(nms) == tolower(target)]
    if (length(hits) == 0) NA_character_ else hits[1]
  }

  normalize_spaces <- function(x) {
    x <- gsub("_", " ", as.character(x))
    x <- gsub("\\s+", " ", x)
    trimws(x)
  }

  final_df <- sequences_df
  blast_parts <- base::strsplit(as.character(normalize_spaces(final_df$Hit_Definition)), " ")
  blast_genus <- vapply(blast_parts, function(words) {
    if (length(words) >= 1) words[1] else NA_character_
  }, character(1))
  blast_species_word <- vapply(blast_parts, function(words) {
    if (length(words) >= 2) words[2] else NA_character_
  }, character(1))
  blast_species_full <- ifelse(
    is_missing_taxon(blast_genus) | is_missing_taxon(blast_species_word),
    NA_character_,
    paste(blast_genus, blast_species_word)
  )

  tax_cols <- names(final_df)
  phylum_col <- find_tax_col(tax_cols, "Phylum")
  class_col <- find_tax_col(tax_cols, "Class")
  order_col <- find_tax_col(tax_cols, "Order")
  family_col <- find_tax_col(tax_cols, "Family")
  genus_col <- find_tax_col(tax_cols, "Genus")
  species_col <- find_tax_col(tax_cols, "Species")

  fill_cols <- c(phylum_col, class_col, order_col, family_col, genus_col)
  fill_cols <- fill_cols[!is.na(fill_cols)]

  if (length(fill_cols) > 0) {
    prev_values <- NULL
    for (col_name in fill_cols) {
      current_values <- as.character(final_df[[col_name]])
      current_values <- normalize_spaces(current_values)
      if (!is.null(prev_values)) {
        missing_idx <- is_missing_taxon(current_values)
        current_values[missing_idx] <- prev_values[missing_idx]
      }
      final_df[[col_name]] <- current_values
      prev_values <- current_values
    }
  }

  if (!is.na(species_col) && !is.na(genus_col)) {
    genus_values <- normalize_spaces(final_df[[genus_col]])
    species_values <- normalize_spaces(final_df[[species_col]])

    species_is_genus_na <- (!is_missing_taxon(genus_values)) & (
      is_missing_taxon(species_values) |
        tolower(species_values) == tolower(paste(genus_values, "NA")) |
        tolower(species_values) == tolower(paste(genus_values, "sp."))
    )

    blast_genus_match <- (!is_missing_taxon(genus_values)) &
      (!is_missing_taxon(blast_genus)) &
      (tolower(genus_values) == tolower(blast_genus))

    use_blast_species <- species_is_genus_na & blast_genus_match & !is_missing_taxon(blast_species_word)
    species_source <- rep("dada2", nrow(final_df))
    species_source[species_is_genus_na] <- "genus_only"
    species_source[use_blast_species] <- "blast_corrected"

    species_values[use_blast_species] <- paste(genus_values[use_blast_species], blast_species_word[use_blast_species])

    unresolved_species <- is_missing_taxon(species_values) & !is_missing_taxon(genus_values)
    species_values[unresolved_species] <- paste(genus_values[unresolved_species], "NA")

    final_df[[species_col]] <- species_values
    final_df$Blast_Genus <- blast_genus
    final_df$Blast_Species <- blast_species_full
    final_df$Species_Source <- species_source
  }

  original_output <- sprintf("1_out/%s.updated_sequences_with_blast_results.%s.csv", project, format(Sys.Date(), "%y%m%d"))
  final_output <- sprintf("1_out/%s.final_asvTable.%s.csv", project, format(Sys.Date(), "%y%m%d"))

  write.csv(sequences_df, original_output, row.names = TRUE)
  write.csv(final_df, final_output, row.names = TRUE)

  cat("Results have been successfully saved to:\n")
  cat("- ", original_output, "\n", sep = "")
  cat("- ", final_output, "\n", sep = "")

  invisible(list(raw_output = original_output, final_output = final_output))
}
