#' Calculate Alpha Diversity Metrics for Phyloseq Data
#'
#' This function calculates various alpha diversity metrics for microbiome data
#' provided in a Phyloseq object. It supports multiple diversity measures and
#' can include phylogenetic diversity (Faith's PD).
#'
#' @param psIN A Phyloseq object containing OTU/ASV counts and associated sample data.
#' @param project Name or identifier for the project or dataset.
#' @param alpha_metrics A character vector specifying the alpha diversity measures to be calculated.
#'                      Supported measures include standard diversity indices and Faith's PD.
#'
#' @return A data frame containing calculated alpha diversity metrics for each sample
#'         in the Phyloseq object, merged with the sample data.
#'
#' @examples
#' # Assuming 'ps' is a Phyloseq object and 'project_name' is the project identifier
#' Go_adiv_markdown(ps, project_name, c("Shannon", "Simpson", "PD"))
#'
#' @export
#' @importFrom phyloseq estimate_richness sample_data
#' @importFrom picante pd

Go_adiv_markdown <- function(psIN, project, alpha_metrics){
  
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
  
  return(adiv)
} 

