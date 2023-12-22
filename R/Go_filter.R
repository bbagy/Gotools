
#' Filter Low Abundance Taxa from Phyloseq Object
#'
#' @param psIN Phyloseq object to be filtered.
#' @param cutoff Abundance threshold for filtering taxa.
#'
#' @details
#' This function filters out low-abundance taxa from a Phyloseq object based on a specified abundance threshold. It first removes samples with zero or very low total counts, then applies a relative abundance filter to exclude taxa below the cutoff threshold.
#'
#' @return
#' A filtered Phyloseq object.
#'
#' @examples
#' ps_filtered <- Go_filter(psIN = my_phyloseq_object, cutoff = 0.00005)
#'
#' @export

Go_filter <- function(psIN, cutoff){ #project


  # remove 0 ASVs
  
  tt = try(psIN.prune <- prune_samples(sample_sums(psIN) > 1, psIN),T)
  if (class(tt) == "try-error"){
    psIN.prune = prune_samples(sample_sums(psIN) > 0, psIN);psIN.prune
  }else{
    psIN.prune <- prune_samples(sample_sums(psIN) > 1, psIN);psIN.prune
  }
  
  
  x <- sample_sums(psIN.prune) > 1
  cat("#--  Removing 0 sum column   --#\n")
  cat(sprintf("#--  %s column(s) are removed   --#\n", length(x[x== F])))
  cat("\n")

  phylo_relabun <- transform_sample_counts(psIN.prune, function(x) x / sum(x))
  phylo_filter <- filter_taxa(phylo_relabun, function(x) mean(x) < cutoff, TRUE) #.00005
  rmtaxa <- taxa_names(phylo_filter)
  alltaxa <- taxa_names(phylo_relabun)
  myTaxa <- alltaxa[!alltaxa %in% rmtaxa]
  phylo_relabun_filtered <- prune_taxa(myTaxa,phylo_relabun)
  ps_filtered <- prune_taxa(myTaxa,psIN.prune)

  cat("#--  Before filter  --#\n")
  print(psIN)
  cat("\n")
  cat("#--  After filter   --#\n")
  print(ps_filtered)

  prune_taxa(myTaxa,psIN)
  
  # out dir
  # out <- file.path("2_rds") 
  # if(!file_test("-d", out)) dir.create(out)
  #saveRDS(ps_filtered, sprintf("%s/ps_filtered.%s.(%s).%s.rds", out, project, cutoff,format(Sys.Date(), "%y%m%d")))
  return(ps_filtered)
  #cat("\n")
  #print(sprintf("ps_filtered is saved as 2_rds/ps_filtered.%s.(%s).%s.rds",  project, cutoff,format(Sys.Date(), "%y%m%d")))

}
