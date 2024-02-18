
#' Perform Paired Permutational Multivariate Analysis of Variance (PERMANOVA)
#'
#' This function conducts a paired PERMANOVA to assess differences in microbial community composition
#' between groups, taking into account paired samples or repeated measures.
#'
#' @param psIN A `phyloseq` object containing the microbiome data.
#' @param cate.vars A vector of categorical variables to be tested.
#' @param project A string indicating the project name for directory and file naming.
#' @param distance A `phyloseq::distance` object or the method to calculate distances.
#' @param distance_metrics A vector of distance metrics to use for the PERMANOVA.
#' @param cate.conf (Optional) A vector of confounding variables to adjust for in the model.
#' @param des A description or label for the PERMANOVA run.
#' @param name (Optional) A name for the analysis for labeling and file naming.
#'
#' @details
#' The function performs a PERMANOVA analysis on each pair of levels within each categorical variable.
#' It can adjust for confounding variables and allows for multiple distance metrics.
#'
#' @return
#' Outputs the results as CSV files in a specified directory. The function also returns the results as
#' a data frame.
#'
#' @examples
#' Go_pairedperm(psIN = ps,
#'               cate.vars = c("TreatmentGroup"),
#'               project = "MicrobiomeStudy",
#'               distance = "bray",
#'               distance_metrics = c("bray", "unifrac"),
#'               cate.conf = c("Age", "BMI"),
#'               des = "TreatmentComparison",
#'               name = "TreatmentAnalysis")
#'
#' @export

Go_pairedperm <- function(psIN, cate.vars, project, distance_metrics, cate.conf=NULL, des, name=NULL){
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
  if (!is.null(des)) {
    # Uni
    print(sprintf("#--- Running Paired-PERMANOVA (%s) ---#", des))
  }  else {
    print("#--- Running Paired-PERMANOVA  ---#")
  }
  set.seed(1)
  mapping.sel <-data.frame(sample_data(psIN))

  res.pair <-{}
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

      distance <- Go_dist(psIN = psIN.sel, project = project, cate.vars=mvar, distance_metrics = distance_metric)


      # pairwise.adonis2
      # pair.ado <- pairwise.adonis2(x=as.dist(distance[[distance_metric]]), factors = mapping.sel.na[,mvar], map=mapping.sel.na, cate.conf=adjust, mvar=mvar)

      x <- as.dist(distance[[distance_metric]])
      factors <-  mapping.sel.na[,mvar]
      map <- mapping.sel.na

      co <- combn(unique(as.character(map[,mvar])),2)
      R2 <- c()
      p.value <- c()
      F.Model <- c()
      pairs <- c()
      SumsOfSqs <- c()
      Df <- c()





      for(elem in 1:ncol(co)){
        x1=as.matrix(x)[factors %in% c(as.character(co[1,elem]),as.character(co[2,elem])),
                        factors %in% c(as.character(co[1,elem]),as.character(co[2,elem]))]

        # run
        map.pair <- subset(map, map[,mvar] %in% c(co[1,elem],co[2,elem]))

# count to table
        if (table(map.pair[,mvar])[co[1,elem]] <=2 | table(map.pair[,mvar])[co[2,elem]] <=2){
          next
        }

        if (!is.null(cate.conf)) {
          form <- as.formula(sprintf("x1 ~ %s + %s", mvar, paste(setdiff(cate.conf, "SampleType"), collapse="+")))
          print(form)
        }else{
          form <- as.formula(sprintf("x1 ~ %s", mvar))
          print(form)
        }

        ad <- adonis2(form, data = map.pair, permutations=999, by="terms")# "terms"  "margin" NULL

        pairs <- c(pairs,paste(co[1,elem],'vs',co[2,elem]));
        Df <- c(Df,ad[1,1])
        SumsOfSqs <- c(SumsOfSqs, ad[1,2])
        R2 <- c(R2,ad[1,3])
        F.Model <- c(F.Model,ad[1,4]);
        p.value <- c(p.value,ad[1,5])
      }

      pairw.res <- data.frame(pairs,Df,SumsOfSqs,R2,F.Model,p.value)

      class(pairw.res) <- c("pwadonis", "data.frame")
      # end adonis end
      tmp <- as.data.frame(pairw.res)
      tmp$distance_metric <- distance_metric
      tmp$mvar <- mvar
      tmp$adjusted <- paste(setdiff(cate.conf, "SampleType"), collapse="+")
      res.pair <- rbind(res.pair, tmp)
    }
  }

  res.pair$padj <- p.adjust(res.pair$p.value, method="bonferroni")

  res.pair <- res.pair[,c("pairs", "Df","SumsOfSqs","R2","F.Model", "p.value", "padj", "distance_metric","mvar", "adjusted")]


  # output
    write.csv(res.pair, quote = FALSE,col.names = NA, sprintf("%s/pair_permanova.%s.%s%s%s%s.csv", out_perm,
              project,
              ifelse(is.null(cate.conf), "", paste(cate.conf, "adjusted.", sep = "")),
              ifelse(is.null(des), "", paste(des, ".", sep = "")),
              ifelse(is.null(name), "", paste(name, ".", sep = "")),
              format(Sys.Date(), "%y%m%d")),sep="/")


  return(res.pair)
}
