
#' Alpha Diversity Analysis for Microbiome Data
#'
#' This function computes various alpha diversity metrics for microbiome data.
#'
#' @param psIN A phyloseq object containing microbiome data.
#' @param project A string representing the project name, used for file naming and directory creation.
#' @param alpha_metrics A vector of alpha diversity metrics to be calculated.
#'
#' @details
#' The function calculates several alpha diversity indices, including Observed, Chao1, ACE, Shannon, Simpson, InvSimpson, Fisher, and PD (Phylogenetic Diversity: Faith’s PD). It allows for a comprehensive analysis of species richness and evenness within individual samples in a microbiome dataset.
#'
#' @return
#' Returns a data frame containing the calculated alpha diversity indices for each sample in the dataset. The result includes both the diversity indices and the corresponding sample metadata.
#'
#' @examples
#' # Example usage:
#' adiv_results <- Go_adiv(psIN = my_phyloseq_object,
#'                         project = "MyMicrobiomeStudy",
#'                         alpha_metrics = c("Observed", "Shannon", "Simpson"))
#'
#' @export

Go_adiv <- function(psIN, project, alpha_metrics, name=NULL){
  # out dir
  out <- file.path(sprintf("%s_%s",project, format(Sys.Date(), "%y%m%d")))
  if(!file_test("-d", out)) dir.create(out)
  out_path <- file.path(sprintf("%s_%s/table",project, format(Sys.Date(), "%y%m%d")))
  if(!file_test("-d", out_path)) dir.create(out_path)
  out_adiv <- file.path(sprintf("%s_%s/table/adiv",project, format(Sys.Date(), "%y%m%d")))
  if(!file_test("-d", out_adiv)) dir.create(out_adiv)

  print("You can measure Observed, Chao1, ACE, Shannon, Simpson, InvSimpson, Fisher, and PD (Phylogenetic diversity:Faith’s PD).")

  # adiv table
  adiv <- estimate_richness(psIN, measures=alpha_metrics) # se.chao1 stand error
  mapping.sel <- data.frame(sample_data(psIN))
  rownames(adiv) <- gsub("^X", "", rownames(adiv))
  adiv$SampleID <- rownames(adiv)
  rownames(adiv) <- rownames(mapping.sel)


  # add pd (Phylogenetic diversity:Faith’s PD) to adiv table
  pd <- grep("PD", alpha_metrics)

  if (length(pd) > 0){
    psIN.tab <- as.data.frame(psIN@otu_table)
    psIN.tree <- psIN@phy_tree

    tt <- try(df.pd <- pd(t(psIN.tab), psIN.tree, include.root=T),T)

    if(class(tt) =="try-error"){
      df.pd <- pd(psIN.tab, psIN.tree, include.root=T)
    }else{
      df.pd <- pd(t(psIN.tab), psIN.tree, include.root=T)
    }

    adiv <- merge(df.pd, adiv, by="row.names")
    rownames(adiv) <- adiv$Row.names
    adiv$Row.names <- NULL
  }

  # add pd mapping data to adiv table
  adiv <- merge(adiv, mapping.sel, by="row.names")
  rownames(adiv) <- adiv$Row.names
  adiv$Row.names <- NULL





  cat(sprintf("adiv table is saved in %s.\n",out_path))
  cat("                                                       \n")
  write.csv(adiv, quote = FALSE, col.names = NA,
            file=sprintf("%s/adiv.%s.%s%s.csv",out_adiv,
                         project,
                         ifelse(is.null(name), "", paste(name, ".", sep = "")),
                         format(Sys.Date(), "%y%m%d"),sep="/"))
  return(adiv)
}

