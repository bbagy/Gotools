

#' Generate Grouped Boxplots from Phyloseq Data
#'
#' This function creates grouped boxplots for Phyloseq data, useful for comparing distributions
#' across different groups or conditions.
#'
#' @param psIN Phyloseq object containing the data for analysis.
#' @param mainGroup The main grouping variable for boxplot comparison.
#' @param project Name of the project or analysis.
#' @param orders Optional ordering of levels in the main grouping variable.
#' @param top Number of top taxa to include in the analysis (if NULL, all taxa are considered).
#' @param name Optional name for the analysis.
#' @param rank Taxonomic rank to consider for boxplot generation.
#' @param pval P-value threshold for significance in the analysis.
#' @param standardsize Logical value indicating whether to standardize plot size.
#' @param mycols Color specifications for the plot.
#' @param ylim Y-axis limits for the plot.
#' @param flip Logical value to flip the axis of the plot.
#' @param height Height of the output plot.
#' @param width Width of the output plot.
#'
#' @return Generates grouped boxplot as a PDF file and saves it in a specified directory.
#'
#' @details
#' The function performs statistical tests (Kruskal-Wallis or Wilcoxon) to determine significant differences
#' between groups. It then generates boxplots for taxa abundances across the specified groups.
#'
#' @examples
#' # psIN is a Phyloseq object
#' # Example usage:
#' Go_groupBox(psIN = psIN,
#'             mainGroup = "Treatment",
#'             project = "MyProject",
#'             orders = c("Control", "Treatment1", "Treatment2"),
#'             top = 10,
#'             name = "Analysis1",
#'             rank = "Genus",
#'             pval = 0.05,
#'             height = 8,
#'             width = 10)
#'
#' @export

Go_groupBox <- function(psIN, mainGroup, project,
                        orders=NULL,
                        top=NULL,
                        name =NULL,
                        rank,
                        pval,
                        standardsize=TRUE,
                        mycols=NULL,
                        ylim=NULL,
                        flip=T,
                        height,
                        width){

  if (!requireNamespace("compositions", quietly = TRUE))
    install.packages("compositions")
  library(compositions)



  if(!is.null(dev.list())) dev.off()
  # out dir
  out <- file.path(sprintf("%s_%s",project, format(Sys.Date(), "%y%m%d")))
  if(!file_test("-d", out)) dir.create(out)
  out_path <- file.path(sprintf("%s_%s/pdf",project, format(Sys.Date(), "%y%m%d")))
  if(!file_test("-d", out_path)) dir.create(out_path)
  set.seed(151)

  # out file
  # "name" definition
  if (class(name) == "function"){
    name <- NULL
  }

  tt <- try(mycols,T)
  if(class(tt) == "try-error"){
    print("mycols is not defined.")
    mycols <- NULL
  }

  tt <- try(orders,T)
  if(class(tt) == "try-error"){
    print("orders is not defined.")
    orders <- NULL
  }


  if(!is.null(top)){
    Top = names(sort(taxa_sums(psIN), TRUE)[1:top])
    ps.top = prune_taxa(Top, psIN);ps.top
  }else{
    ps.top = psIN
  }


  ### Centered Log Ratio (CLR)
  otu_data <- otu_table(ps.top)
  # Replace zeros if necessary, since CLR cannot handle zeros
  otu_data[otu_data == 0] <- 1e-5

  # Apply the CLR transformation
  clr_transformed <- clr(otu_data)

  # Convert back to a matrix if needed (clr returns an array)
  clr_transformed <- matrix(clr_transformed, nrow = nrow(otu_data), ncol = ncol(otu_data), dimnames = dimnames(otu_data))

  # Update the OTU table in the phyloseq object
  otu_table(ps.top) <- otu_table(clr_transformed, taxa_are_rows = taxa_are_rows(ps.top))


  tab <- data.frame(otu_table(ps.top));head(tab)

  nsamps_threshold <- 0.01 # fraction of relabund to call a sample positive
  filt_threshold <- 0.1 # fraction of samples that need to be positive to keep an OTU for association testing
  nperm <- 100000

  ### aggregation by rank
  otu.filt <- as.data.frame((otu_table(ps.top))) # for dada2  t(otu_table(ps.relative)
  otu.filt$func <- getTaxonomy(otus=rownames(otu.filt), tax_tab=tax_table(ps.top), taxRanks =colnames(tax_table(psIN)), level= rank)
  agg <- aggregate(. ~ func, otu.filt, sum, na.action=na.pass);dim(agg)


  funcNames <- agg$func
  rownames(agg) <- funcNames
  #ftk <- names(which(unlist(apply(agg, 1, function(x) length(which(x>=nsamps_threshold)))) > ceiling(filt_threshold*ncol(agg))))
  # agg <- agg[intersect(ftk,ftk),]
  # agg_t <- t(agg)



  aggDf<-as.data.frame(agg, row.names = agg$func)
  aggDf$func <- NULL
  agg_t <- t(aggDf)

  colnames(agg_t)
  rownames(agg_t)

  # Add grouping information
  map <- sample_data(ps.top);dim(map)
  df <- data.frame(agg_t, Group = map[,mainGroup]) #, name = map$StudyID,  NoOfFMT= map$NoOfFMT );head(df)

  df[,mainGroup] <- as.character(df[,mainGroup]);df[,mainGroup]
  df[,mainGroup][df[,mainGroup]==""] <- "NA";df[,mainGroup]
  df.na <- subset(df, df[,mainGroup] != "NA");df.na[,mainGroup]  # subset 를 사용한 NA 삭제
  df.na[,mainGroup] <- as.factor(df.na[,mainGroup]);df.na[,mainGroup]



  #=========== KW or WX test for function screening =========#
  df$Group <- NULL

  set.seed(151)
  # Ensure that 'maingroup' is the last column in your data frame
  group_1 <- as.factor(df.na[,mainGroup]); group_1

  # Exclude 'maingroup' column from the dataframe for the test
  df_test <- df.na[, -which(names(df.na) == mainGroup)]

  result.table <- data.frame()
  for (i in 1:dim(df_test)[2]) {
    # Detect number of groups
    num.groups <- length(unique(group_1))

    if (num.groups > 2) {
      test.result <- kruskal.test(df_test[,i], g=group_1)
      sig.test <- "KW"
      # Report number of values tested
      cat(paste("Kruskal-Wallis test for ",names(df_test)[i]," ", i, "/",
                dim(df_test)[2], "; p-value=", test.result$p.value,"\n", sep=""))
    } else if (num.groups == 2) {
      test.result <- wilcox.test(df_test[,i], g=group_1)
      sig.test <- "WX"
      # Report number of values tested
      cat(paste("Wilcoxon test for ",names(df_test)[i]," ", i, "/",
                dim(df_test)[2], "; p-value=", test.result$p.value,"\n", sep=""))
    } else {
      next
    }
    # Store the result in the data frame
    result.table <- rbind(result.table,
                          data.frame(id=names(df_test)[i],
                                     p.value=test.result$p.value
                          ))
  }




  result.table <- result.table[order(result.table$p.value, decreasing = FALSE), ] # Order by p-value

  result.table <- result.table %>%
    mutate(padj_BH = p.adjust(p.value, method = "BH"))

  sig.result <- result.table[which(result.table$p.value < pval),]

  # Reporting the number of significant results
  cat(paste("There are ", nrow(sig.result), " significant results at p < ", pval, "\n", sep=""))

  sig.mat <- as.matrix(sig.result)
  funcNames.sig <- sig.mat[,1]

  df.sel <- df.na[funcNames.sig]
  df.sel <- data.frame(df.sel, Group = map[,mainGroup])




  df.sel.melt <- melt(df.sel, id.vars = mainGroup, measure.vars = funcNames.sig)
  df.sel.melt$value <- as.numeric(df.sel.melt$value)
  df.sel.melt.clean <- subset(df.sel.melt, variable != "Group" &  value > 0)




  if (!is.null(orders)) {
    df.sel.melt.clean[,mainGroup] <- factor(df.sel.melt.clean[,mainGroup], levels = rev(orders))
  } else {
    df.sel.melt.clean[,mainGroup] <- factor(df.sel.melt.clean[,mainGroup])
  }

  df.sel.melt.clean$variable <- as.character(df.sel.melt.clean$variable)

  df.sel.melt.clean <- df.sel.melt.clean[order(df.sel.melt.clean$variable ,  decreasing = F), ]


  df.sel.melt.clean <- subset(df.sel.melt.clean, variable != mainGroup)

  print(unique(df.sel.melt.clean$variable))
  p <- ggplot(df.sel.melt.clean, aes_string(x="variable", y="value", fill=mainGroup)) +  geom_boxplot(outlier.shape = NA,lwd=0.3) +
    theme_bw() + theme(strip.background = element_blank()) +
    labs(y="Centered Log Ratio (CLR)", x= NULL) + ggtitle(sprintf("%s p < %s", sig.test,pval))

  # + stat_compare_means(aes_string(group = mainGroup),label = "p.format") +

  #+ scale_x_discrete(limits = rev)


  if(!is.null(mycols)){
    p <- p + scale_fill_manual(values = mycols)
  }else{
    p <- p
  }

  if(flip == T){
    p <- p+ coord_flip()
  }else{
    p <- p + theme(text=element_text(size=9), axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
  }



  # tt <- try(mycols, T)
  # if(class(tt) == "try-error"){
  #  p <- p
  # }else{
  #   p <- p + scale_fill_manual(values = mycols)
  # }

  if(!is.null(ylim)){
    p = p + ylim(ylim[1] , ylim[2])
  }else{
    p=p
  }

  #=== image size ===#
  #height <- 0.4*length(unique(df.sel.melt.clean[,mainGroup])) + 0.4*dim(kw.sig)[1];height
  #width <- log((max(nchar(funcNames.sig)))*max(nchar(as.character(unique(df.sel.melt.clean[,mainGroup])))));width

  # plot size ratio
  if (length(unique(df.sel.melt.clean$variable)) < 6){
    if(standardsize==TRUE){
      num.subgroup <- length(unique(df.sel.melt.clean$variable))*0.08
    }else{
      num.subgroup <- 1
    }
  }else{
    num.subgroup <- length(unique(df.sel.melt.clean$variable))*0.08
  }

  p  <- p  + theme(panel.background = element_rect(fill = "white", colour = "grey50"),aspect.ratio = num.subgroup/1)

  #plotlist[[length(plotlist)+1]] <- p
  pdf(sprintf("%s/groupBox.%s.%s.%s%s%s.pdf", out_path,
              project,
              mainGroup,
              ifelse(is.null(rank), "", paste(rank, ".", sep = "")),
              ifelse(is.null(name), "", paste(name, ".", sep = "")),
              format(Sys.Date(), "%y%m%d")), height = height, width = width)

  print(p)
  dev.off()
  return(result.table)
}








