
#' Global Permutation Analysis for Microbiome Data
#'
#' This function performs global permutation (PERMANOVA) analysis on microbiome data.
#'
#' @param psIN A phyloseq object containing microbiome data.
#' @param cate.vars A vector of categorical variables for which PERMANOVA will be performed.
#' @param project A string representing the project name, used for file naming and directory creation.
#' @param distance_metrics A vector of distance metrics to be used in the analysis.
#' @param mul.vars A boolean value indicating whether to use multiple variables in the PERMANOVA model.
#' @param name An optional string for naming the output files.
#'
#' @details
#' The function conducts PERMANOVA analysis to understand the impact of various categorical variables on microbial community composition. It supports multiple distance metrics and can handle both single and multiple variables in the model. The function outputs the results in a CSV file and returns them as a data frame.
#'
#' @return
#' Returns a data frame containing the results of the PERMANOVA analysis, including degrees of freedom, sums of squares, R-squared values, F-model statistics, p-values, and adjusted p-values.
#'
#' @examples
#' # Example usage:
#' permanova_results <- Go_perm(psIN = my_phyloseq_object,
#'                              cate.vars = c("TreatmentGroup", "AgeGroup"),
#'                              project = "MyMicrobiomeStudy",
#'                              distance_metrics = c("bray", "unifrac"),
#'                              mul.vars = FALSE,
#'                              name = "MyAnalysis")
#'
#' @export

Go_perm <- function(psIN, cate.vars, project, distance_metrics, mul.vars=FALSE, name=NULL){
  # out dir
  out <- file.path(sprintf("%s_%s",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out)) dir.create(out)
  out_path <- file.path(sprintf("%s/table",out)) 
  if(!file_test("-d", out_path)) dir.create(out_path)
  out_perm <- file.path(sprintf("%s/perm",out_path)) 
  if(!file_test("-d", out_perm)) dir.create(out_perm)
  
  out_distance <- file.path(sprintf("%s/distance",out_path)) 
  if(!file_test("-d", out_distance)) dir.create(out_distance)
  

  # Run
  set.seed(1)
  mapping.sel <-data.frame(sample_data(psIN))
  
  res <-{}
  # Run
  for (mvar in cate.vars) {
    mapping.sel.na <- mapping.sel[!is.na(mapping.sel[,mvar]), ]
    
    
    if (length(unique(mapping.sel.na[,mvar])) == 1){
      cat(sprintf("there is no group campare to %s\n",unique(mapping.sel[,mvar])))
      next
    }
    for (distance_metric in distance_metrics) {
      
      psIN.sel <- prune_samples(rownames(mapping.sel[!is.na(mapping.sel.na[,mvar]), ]), psIN)
      
      
      
      ## fix factor  and  numeric
      mapping.sel.na[,mvar] <- factor(mapping.sel.na[,mvar])
      
      distance <- Go_dist(psIN = psIN.sel, project = project, distance_metrics = distance_metric)
      
      x <- as.dist(distance[[distance_metric]])
      
      
      
      # write.csv(as.data.frame(x), quote = FALSE,col.names = NA, sprintf("%s/distance.%s.%s.%s%s.csv", out_distance, 
      #                                                     project, 
      #                                                     distance_metric, 
      #                                                     ifelse(is.null(name), "", paste(name, ".", sep = "")), 
      #                                                     format(Sys.Date(), "%y%m%d")),sep="/")
      
      
      
      factors <-  mapping.sel.na[,mvar]
      
      R2 <- c()
      p.value <- c()
      F.Model <- c()
      pairs <- c()
      SumsOfSqs <- c()
      Df <- c()
      
      x1=as.matrix(x)[factors %in% unique(factors), factors %in% unique(factors)]
      
      # run
      map.pair <- subset(mapping.sel.na, mapping.sel.na[,mvar] %in% unique(factors))
      
      # count to table
      
      if (mul.vars == T) {

        form <- as.formula(sprintf("x1 ~ %s", paste(setdiff(cate.vars, "SampleType"), collapse="+")))
        print(form)
      }else{
        form <- as.formula(sprintf("x1 ~ %s", mvar))
        print(form)
      }
      
      ad <- adonis2(form, data = map.pair, permutations=999, by="terms")# "terms"  "margin" NULL
      
      plot(ad)
      
      
      Df <- c(Df,ad[1,1])
      SumsOfSqs <- c(SumsOfSqs, ad[1,2])
      R2 <- round(c(R2,ad[1,3]), digits=3)
      F.Model <- c(F.Model,ad[1,4]);
      p.value <- c(p.value,ad[1,5])
      
      pairw.res <- data.frame(Df,SumsOfSqs,R2,F.Model,p.value)
      
      class(pairw.res) <- c("pwadonis", "data.frame")
      # end adonis end
      tmp <- as.data.frame(pairw.res)
      tmp$distance_metric <- distance_metric

      
      # end adonis end
      tmp <- as.data.frame(pairw.res)
      tmp$distance_metric <- distance_metric
      
      if(mul.vars == TRUE){
        mul.form <- sprintf("x1 ~ %s", paste(setdiff(cate.vars, "SampleType"), collapse="+"))
        tmp$formula <- mul.form
        type <- "mulit"
      } else {
        tmp$mvar <- mvar
        type <- "uni"
      }
      
      res <- rbind(res, tmp)
      

    }
    if(mul.vars == TRUE){
      break
    }
  }
  
  res$padj <- p.adjust(res$p.value, method="bonferroni")
  
 
  if(mul.vars == TRUE){
    res <- res[,c("Df","SumsOfSqs","R2","F.Model", "p.value", "padj", "distance_metric","formula")]
  } else {
    res <- res[,c("Df","SumsOfSqs","R2","F.Model", "p.value", "padj", "distance_metric","mvar")]
  }
  
  # output
    write.csv(res, quote = FALSE,col.names = NA, sprintf("%s/global_permanova.%s.%s.%s%s.csv", out_perm, 
              project, 
              type, 
              ifelse(is.null(name), "", paste(name, ".", sep = "")), 
              format(Sys.Date(), "%y%m%d")),sep="/")
  

  return(res)
}
