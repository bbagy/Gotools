#' Calculate Distances for Microbiome Data
#'
#' This function calculates distances for microbiome data in a Phyloseq object using specified distance metrics.
#' It supports multiple categorical variables and distance metrics.
#'
#' @param psIN A Phyloseq object containing the OTU/ASV counts and associated sample data.
#' @param project Name or identifier for the project or dataset being analyzed.
#' @param cate.vars Vector of column names in the sample data of 'psIN' representing categorical variables for distance calculation.
#' @param name Optional name for the analysis output.
#' @param distance_metrics Vector of distance metrics to be used for the calculations, such as "bray", "unifrac", etc.
#'
#' @return A list of distance matrices, one for each specified distance metric.
#'
#' @examples
#' # Assuming 'ps' is a phyloseq object with appropriate data
#' Go_dist_markdown(ps, "my_project", c("Group"), distance_metrics = c("bray", "unifrac"))
#'
#' @export
#' @importFrom phyloseq phyloseq sample_data distance
#' @importFrom reshape2 melt

Go_dist_markdown <- function(psIN, project, cate.vars, name=NULL, distance_metrics){
  # run distance
  dm <- list()
  for (distance_metric in distance_metrics) {
    dm[[length(dm)+1]] <- phyloseq::distance(psIN, method=distance_metric)
    
    names(dm) <- distance_metric
    
    for (mvar in cate.vars) {
      map <- data.frame(sample_data(psIN))
      
      if (length(unique(map[,mvar])) < 2){
        next
      }
      
      dm1 <- as.dist(dm[[distance_metric]])
      dm1.mat  <- as.matrix(dm1)
      dm1.df <- data.frame(dm1.mat)
      
      df <- melt(as.matrix(dm1), varnames = c("row", "col"))
      
      
      sub_dist <- list()
      groups_all <- map[,mvar]
      
      for (group in unique(groups_all)) { 
        row_group <- which(groups_all == group)
        sample_group <- sample_names(psIN)[row_group]
        sub_dist[[group]] <- dm1.mat[sample_group, sample_group]
        sub_dist[[group]][!lower.tri(sub_dist[[group]])] <- NA
      }
      
      dmgroups<- melt(sub_dist)
      dmgroups.df <- dmgroups[complete.cases(dmgroups), ]
      dmgroups.df[,mvar] <- factor(dmgroups.df$L1, levels=names(sub_dist))
      dmgroups.df$L1 <- NULL
      
      colnames(dmgroups.df) <- c("Var1",	"Var2",	distance_metric, mvar)
      
      # for across samples
      # Initialize a list to store the pairwise distances
      group_dists <- list()
      
      # Get unique group labels
      unique_groups <- unique(groups_all)
      
      # Loop over all pairs of groups
      for (i in 1:(length(unique_groups)-1)) {
        for (j in (i+1):length(unique_groups)) {
          
          # Get the indices of samples in each group
          row_group1 <- which(groups_all == unique_groups[i])
          row_group2 <- which(groups_all == unique_groups[j])
          
          # Get the sample names for each group
          sample_group1 <- sample_names(psIN)[row_group1]
          sample_group2 <- sample_names(psIN)[row_group2]
          
          # Get the pairwise distances between groups
          group_dist <- dm1.mat[sample_group1, sample_group2]
          
          # Add the distances to the list, with a name indicating the pair of groups
          group_dists[[paste(unique_groups[i], unique_groups[j], sep=".VS.")]] <- group_dist
        }
      }
      
      
      across_dist.groups<- melt(group_dists);head(across_dist.groups)
      colnames(across_dist.groups) <- c("Var1",	"Var2",	distance_metric, mvar);head(across_dist.groups)
      
      merge.dm.df <- rbind(dmgroups.df,across_dist.groups)
    }
  }
  return(dm)
}

