#' Go_kraken2Extraction
#'
#' Parse Kraken2 (and optional Bracken) MPA output into a phyloseq object
#' with unified taxonomy for Bacteria, Virus, and Fungi.
#'
#' @description
#' This function extracts taxonomy and abundance tables from Kraken2/Bracken
#' outputs and converts them into a standardized phyloseq object.
#'
#' @param kraken2 Path to Kraken2 MPA file
#' @param bracken Optional Bracken file (required for Bacteria)
#' @param kingdom One of "Bacteria", "Virus", "Fungi"
#'
#' @return A phyloseq object
#'
#' @author Heekuk Park <hp2523@cumc.columbia.edu>
#' @date 2025-11-01
#'
#' @export
#'
#' ----
#' Personal signature:
#'   Heekuk Park, Columbia University Medical Center
#'   Email: hp2523@cumc.columbia.edu
#'   Date: 2025-11-01


Go_kraken2Extraction <- function(kraken2, bracken=NULL, kingdom=c("Bacteria","Virus","Fungi")) {
  
  library(dplyr)
  library(phyloseq)
  
  kingdom <- match.arg(kingdom)
  
  # ---------------------
  # 0. Helper: taxonomy order standardizer
  # ---------------------
  FixTaxOrder <- function(tax_df, kingdom_name){
    tax_df$Rank1 <- kingdom_name
    tax_df <- tax_df[, c("Rank1","Phylum","Class","Order","Family","Genus","Species")]
    return(tax_df)
  }
  
  # ---------------------
  # 1. Load Kraken2 MPA
  # ---------------------
  kraken <- read.delim(kraken2, sep="\t", header=TRUE, check.names=FALSE)
  
  # ---------------------
  # 2. Case 1: Virus (Bracken 없음)
  # ---------------------
  if (kingdom == "Virus") {
    
    # ---- Viral taxonomy formatter ----
    FixViralTaxonomy <- function(tax_vector) {
      viral_ranks <- c("d__", "k__", "p__", "c__", "o__", "f__", "g__", "s__")
      out_list <- list()
      
      for (i in seq_along(tax_vector)) {
        parts <- unlist(strsplit(tax_vector[i], "\\|"))
        filled <- setNames(rep(NA, length(viral_ranks)), viral_ranks)
        for (p in parts) {
          pref <- substr(p,1,3)
          if (pref %in% viral_ranks) filled[pref] <- p
        }
        out_list[[i]] <- filled
      }
      
      out <- do.call(rbind, out_list)
      colnames(out) <- c("domain","kingdom","phylum","class","order","family","genus","species")
      return(as.data.frame(out))
    }
    
    # ---- Virus만 추출 ----
    virus <- kraken %>% filter(grepl("^d__Viruses", ID))
    virus$ID <- as.character(virus$ID)
    
    # ---- raw taxonomy parsing ----
    tax_raw <- FixViralTaxonomy(virus$ID)
    
    # ---- standard taxonomy columns ----
    tax_std <- tax_raw[, c("phylum","class","order","family","genus","species")]
    colnames(tax_std) <- c("Phylum","Class","Order","Family","Genus","Species")
    
    # ---- standard order + Rank1 ----
    tax_std <- FixTaxOrder(tax_std, "Viruses")
    rownames(tax_std) <- virus$ID
    
    # ---- abundance ----
    otu_mat <- as.matrix(virus[, !(colnames(virus) %in% "ID")])
    rownames(otu_mat) <- virus$ID
    
    otu <- otu_table(otu_mat, taxa_are_rows = TRUE)
    tax <- tax_table(as.matrix(tax_std))
    
    return(phyloseq(otu, tax))
  }
  
  # ---------------------
  # 3. Case 2: Bacteria (Bracken 사용)
  # ---------------------
  if (kingdom == "Bacteria") {
    
    if (is.null(bracken)) stop("Bracken file required for kingdom='Bacteria'")
    
    br <- read.delim(bracken, sep="\t", header=TRUE, check.names=FALSE)
    
    # ---- extract species ----
    kraken$species_id <- sub(".*\\|s__", "s__", kraken$ID)
    kraken$species_id[!grepl("^s__", kraken$species_id)] <- NA
    
    # ---- taxonomy parser ----
    parse_taxonomy <- function(id_string) {
      parts <- unlist(strsplit(id_string, "\\|"))
      get_level <- function(prefix) {
        m <- grep(paste0("^", prefix), parts, value=TRUE)
        if (length(m)==0) return(NA)
        return(sub("^[a-z]__+", "", m))
      }
      return(c(
        phylum  = get_level("p__"),
        class   = get_level("c__"),
        order   = get_level("o__"),
        family  = get_level("f__"),
        genus   = get_level("g__"),
        species = get_level("s__")
      ))
    }
    
    taxonomy_parsed <- as.data.frame(
      do.call(rbind, lapply(kraken$ID, parse_taxonomy)),
      stringsAsFactors = FALSE
    )
    
    taxonomy_df <- cbind(species_id = kraken$species_id, taxonomy_parsed)
    taxonomy_df$species_id <- gsub(" ", "_", taxonomy_df$species_id)
    
    # ---- merge Bracken + taxonomy ----
    merged <- merge(br, taxonomy_df, by.x="ID", by.y="species_id", all.x=TRUE)
    
    # ---- format taxonomy ----
    tax <- merged[, c("phylum","class","order","family","genus","species")]
    colnames(tax) <- c("Phylum","Class","Order","Family","Genus","Species")
    
    tax <- FixTaxOrder(tax, "Bacteria")
    rownames(tax) <- merged$ID
    
    # ---- abundance matrix ----
    otu_cols <- setdiff(colnames(merged), c("ID","phylum","class","order","family","genus","species"))
    otu_mat <- as.matrix(merged[, otu_cols])
    rownames(otu_mat) <- merged$ID
    
    otu <- otu_table(otu_mat, taxa_are_rows = TRUE)
    tax <- tax_table(as.matrix(tax))
    
    return(phyloseq(otu, tax))
  }
  
  # ---------------------
  # 4. Case 3: Fungi (Bracken 없음)
  # ---------------------
  if (kingdom == "Fungi") {
    
    fungi <- kraken %>% filter(grepl("k__Fungi", ID))
    fungi$ID <- as.character(fungi$ID)
    
    # ---- taxonomy parser ----
    parse_taxonomy <- function(id_string) {
      parts <- unlist(strsplit(id_string, "\\|"))
      get_level <- function(prefix) {
        m <- grep(paste0("^", prefix), parts, value=TRUE)
        if (length(m)==0) return(NA)
        return(sub("^[a-z]__+", "", m))
      }
      return(c(
        phylum  = get_level("p__"),
        class   = get_level("c__"),
        order   = get_level("o__"),
        family  = get_level("f__"),
        genus   = get_level("g__"),
        species = get_level("s__")
      ))
    }
    
    tax_raw <- as.data.frame(do.call(rbind, lapply(fungi$ID, parse_taxonomy)))
    colnames(tax_raw) <- c("Phylum","Class","Order","Family","Genus","Species")
    
    tax_std <- FixTaxOrder(tax_raw, "Fungi")
    rownames(tax_std) <- fungi$ID
    
    otu_mat <- as.matrix(fungi[, !(colnames(fungi) %in% "ID")])
    rownames(otu_mat) <- fungi$ID
    
    otu <- otu_table(otu_mat, taxa_are_rows = TRUE)
    tax <- tax_table(as.matrix(tax_std))
    
    
    return(phyloseq(otu, tax))
  }
}

