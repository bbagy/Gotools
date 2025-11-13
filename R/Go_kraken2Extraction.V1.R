#' Go_kraken2Extraction
#'
#' Parse Kraken2 (and optional Bracken) MPA output into a phyloseq object
#' with unified taxonomy for Bacteria, Virus, and Fungi.
#'
#' @description
#' This function extracts taxonomy and abundance tables from Kraken2 and
#' optionally Bracken outputs, then converts them into a standardized
#' phyloseq object with harmonized taxonomic ranks.
#'
#' Supported kingdoms:
#' - "Bacteria": Requires Bracken; merges taxonomy with refined abundance.
#' - "Virus": Kraken-only; missing taxonomic ranks are filled with "Unknown".
#' - "Fungi": Kraken-only; missing ranks also filled with "Unknown".
#'
#' The final taxonomy table always uses the rank order:
#' Rank1, Phylum, Class, Order, Family, Genus, Species.
#'
#' @param kraken2 Path to Kraken2 MPA file.
#' @param bracken Optional Bracken MPA file (required when kingdom = "Bacteria").
#' @param kingdom One of "Bacteria", "Virus", or "Fungi".
#'
#' @return A phyloseq object containing OTU + taxonomy tables.
#'
#' @author Heekuk Park (hp2523@@cumc.columbia.edu)
#'
#' @examples
#' \dontrun{
#'
#' # -----------------------
#' # Example 1: Bacteria
#' # Requires both Kraken2 and Bracken MPA files
#' # -----------------------
#' ps_bac <- Go_kraken2Extraction(
#'   kraken2 = "1_out/20251111_kraken2_mpa.txt",
#'   bracken = "1_out/20251111_bracken_mpa.txt",
#'   kingdom = "Bacteria"
#' )
#' print(ps_bac)
#'
#'
#' # -----------------------
#' # Example 2: Virus
#' # Kraken2 only (Bracken not needed)
#' # -----------------------
#' ps_vir <- Go_kraken2Extraction(
#'   kraken2 = "1_out/20251111_kraken2_mpa.txt",
#'   kingdom = "Virus"
#' )
#' print(ps_vir)
#'
#'
#' # -----------------------
#' # Example 3: Fungi
#' # Kraken2 only (Bracken not needed)
#' # -----------------------
#' ps_fun <- Go_kraken2Extraction(
#'   kraken2 = "1_out/20251111_kraken2_mpa.txt",
#'   kingdom = "Fungi"
#' )
#' print(ps_fun)
#'
#' }
#'
#' @export


Go_kraken2Extraction <- function(kraken2, bracken=NULL, kingdom=c("Bacteria","Virus","Fungi")) {

  library(dplyr)
  library(phyloseq)

  kingdom <- match.arg(kingdom)

  # ----------------------------
  # 1. Kraken2 MPA load
  # ----------------------------
  kraken <- read.delim(kraken2, sep="\t", header=TRUE, check.names=FALSE)
  kraken$ID <- as.character(kraken$ID)


  # ==========================================================
  #  CASE 1: VIRUS  (NO Bracken)
  # ==========================================================
  if (kingdom == "Virus") {

    # Only keep viral rows
    virus <- kraken %>% filter(grepl("^d__Viruses", ID))
    virus$ID <- as.character(virus$ID)

    # Keep ONLY species-level rows
    virus <- virus[grepl("s__", virus$ID), ]

    # --------------------
    # Viral taxonomy parser
    # --------------------
    parse_viral_tax <- function(x) {
      parts <- unlist(strsplit(x, "\\|"))
      get <- function(prefix) {
        m <- grep(paste0("^", prefix), parts, value=TRUE)
        if (length(m)==0) return(NA)
        return(sub("^[a-z]__+", "", m))
      }
      return(c(
        Phylum = get("p__"),
        Class  = get("c__"),
        Order  = get("o__"),
        Family = get("f__"),
        Genus  = get("g__"),
        Species= get("s__")
      ))
    }

    taxonomy <- as.data.frame(do.call(rbind, lapply(virus$ID, parse_viral_tax)))

    # Replace NA → "Unknown"
    taxonomy[is.na(taxonomy)] <- "Unknown"

    taxonomy$Rank1 <- "Viruses"
    taxonomy <- taxonomy[, c("Rank1","Phylum","Class","Order","Family","Genus","Species")]
    rownames(taxonomy) <- virus$ID

    # OTU matrix
    otu_mat <- as.matrix(virus[, setdiff(colnames(virus),"ID")])
    rownames(otu_mat) <- virus$ID

    otu <- otu_table(otu_mat, taxa_are_rows=TRUE)
    tax <- tax_table(as.matrix(taxonomy))

    return(phyloseq(otu,tax))
  }



  # ==========================================================
  #  CASE 2: BACTERIA  (Bracken REQUIRED)
  # ==========================================================
  if (kingdom == "Bacteria") {

    if (is.null(bracken)) stop("Bracken file required for Bacteria.")

    br <- read.delim(bracken, sep="\t", header=TRUE, check.names=FALSE)

    # Extract species id
    kraken$species_id <- sub(".*\\|s__", "s__", kraken$ID)
    kraken$species_id[!grepl("^s__", kraken$species_id)] <- NA

    parse_taxonomy <- function(x) {
      parts <- unlist(strsplit(x,"\\|"))
      get <- function(prefix){
        m <- grep(paste0("^",prefix), parts, value=TRUE)
        if (length(m)==0) return(NA)
        return(sub("^[a-z]__+","",m))
      }
      return(c(
        phylum  = get("p__"),
        class   = get("c__"),
        order   = get("o__"),
        family  = get("f__"),
        genus   = get("g__"),
        species = get("s__")
      ))
    }

    taxonomy <- as.data.frame(
      do.call(rbind, lapply(kraken$ID, parse_taxonomy)),
      stringsAsFactors = FALSE
    )

    taxonomy_df <- cbind(species_id = kraken$species_id, taxonomy)
    taxonomy_df$species_id <- gsub(" ","_", taxonomy_df$species_id)

    merged <- merge(br, taxonomy_df, by.x="ID", by.y="species_id", all.x=TRUE)

    ranks <- c("phylum","class","order","family","genus","species")
    colnames(merged) <- sapply(colnames(merged), function(x){
      if (x %in% ranks) paste0(toupper(substr(x,1,1)), substr(x,2,nchar(x)))
      else x
    })

    # safe rank selection
    safe_ranks <- intersect(ranks, colnames(merged))

    # only fill NA in the ranks that actually exist
    merged[, safe_ranks] <- lapply(merged[, safe_ranks, drop=FALSE], function(x) {
      ifelse(is.na(x), "Unknown", x)
    })

    otu_cols <- setdiff(colnames(merged), c("ID","Phylum","Class","Order","Family","Genus","Species"))
    otu_mat <- as.matrix(merged[, otu_cols])
    rownames(otu_mat) <- merged$ID

    taxonomy <- merged[, c("Phylum","Class","Order","Family","Genus","Species")]
    taxonomy$Rank1 <- "Bacteria"
    taxonomy <- taxonomy[, c("Rank1","Phylum","Class","Order","Family","Genus","Species")]
    rownames(taxonomy) <- rownames(otu_mat)

    otu <- otu_table(otu_mat, taxa_are_rows=TRUE)
    tax <- tax_table(as.matrix(taxonomy))

    return(phyloseq(otu,tax))
  }



  # ==========================================================
  #  CASE 3: FUNGI (NO Bracken)
  # ==========================================================
  if (kingdom == "Fungi") {

    fungi <- kraken %>% filter(grepl("k__Fungi", ID))

    fungi <- fungi[grepl("s__", fungi$ID), ]  # species only

    parse_fun_tax <- function(x){
      parts <- unlist(strsplit(x,"\\|"))
      get <- function(prefix){
        m <- grep(paste0("^",prefix), parts, value=TRUE)
        if (length(m)==0) return(NA)
        return(sub("^[a-z]__+","",m))
      }
      return(c(
        Phylum = get("p__"),
        Class  = get("c__"),
        Order  = get("o__"),
        Family = get("f__"),
        Genus  = get("g__"),
        Species= get("s__")
      ))
    }

    taxonomy <- as.data.frame(do.call(rbind, lapply(fungi$ID, parse_fun_tax)))
    taxonomy[is.na(taxonomy)] <- "Unknown"

    taxonomy$Rank1 <- "Fungi"
    taxonomy <- taxonomy[, c("Rank1","Phylum","Class","Order","Family","Genus","Species")]
    rownames(taxonomy) <- fungi$ID

    otu_mat <- as.matrix(fungi[, setdiff(colnames(fungi),"ID")])
    rownames(otu_mat) <- fungi$ID

    otu <- otu_table(otu_mat, taxa_are_rows=TRUE)
    tax <- tax_table(as.matrix(taxonomy))

    return(phyloseq(otu,tax))
  }

}
