
#' Calculate Distance Matrices for Phyloseq Object
#'
#' This function calculates distance matrices for a given Phyloseq object using various distance metrics.
#'
#' @param psIN A Phyloseq object to calculate the distance matrices for.
#' @param project A string indicating the project name for file naming.
#' @param cate.vars An array of categorical variables to consider in the distance matrix calculations.
#' @param name (Optional) A string to add as a suffix to the file name.
#' @param distance_metrics An array of distance metrics to use for calculations.
#'
#' @details
#' The function calculates distance matrices using specified distance metrics and saves them in specified directories.
#' Supported distance metrics include unifrac, wunifrac, dpcoa, jsd, manhattan, euclidean, canberra, bray, kulczynski, and jaccard.
#' It also handles subsetting for different categories and calculates distances both within and across different groups.
#'
#' @return
#' A list of distance matrices for each specified metric.
#'
#' @examples
#' Go_dist(psIN = my_phyloseq_object,
#'         project = "MyProject",
#'         cate.vars = c("Variable1", "Variable2"),
#'         name = "Analysis1",
#'         distance_metrics = c("bray", "unifrac"))
#'
#' @export

Go_dist <- function(psIN, project, cate.vars, name=NULL, distance_metrics){
  # out dir
  out <- file.path(sprintf("%s_%s",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out)) dir.create(out)
  out_path <- file.path(sprintf("%s_%s/table",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out_path)) dir.create(out_path)
  out_dist <- file.path(sprintf("%s_%s/table/dist",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out_dist)) dir.create(out_dist)
  
  print("You can measure unifrac, wunifrac, dpcoa, jsd, manhattan, euclidean, canberra, bray, kulczynski, and jaccard.")
  
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
      
      
      write.csv(dm1.df, quote = FALSE,col.names = NA, sprintf("%s/distance1.%s.%s.%s%s.csv", out_dist, 
                                                              project, 
                                                              distance_metric, 
                                                              ifelse(is.null(name), "", paste(name, ".", sep = "")), 
                                                              format(Sys.Date(), "%y%m%d")),sep="/")
      
      write.csv(merge.dm.df, quote = FALSE,col.names = NA, sprintf("%s/distance2.%s.%s.%s%s%s.csv", out_dist, 
                                                                   project, 
                                                                   distance_metric, 
                                                                   ifelse(is.null(mvar), "", paste(mvar, ".", sep = "")), 
                                                                   ifelse(is.null(name), "", paste(name, ".", sep = "")), 
                                                                   format(Sys.Date(), "%y%m%d")),sep="/")
      
    

    }
  }
  return(dm)
}

