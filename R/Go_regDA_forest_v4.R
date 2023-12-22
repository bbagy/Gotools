
#' Generate Forest Plots from Regression Differential Abundance Analysis
#'
#' This function generates forest plots from regression differential abundance (regDA) analysis.
#'
#' @param project Name of the project or analysis.
#' @param file_path Directory path containing input files for analysis.
#' @param mycols Optional color palette for the plot.
#' @param fdr False Discovery Rate threshold for significance.
#' @param estimate Minimum estimate threshold for including in the plot.
#' @param files Pattern to match filenames for analysis.
#' @param name Optional name for the analysis.
#' @param orders Optional order of levels in the categorical variables.
#' @param height Height of the output plot.
#' @param width Width of the output plot.
#'
#' @return Generates a series of forest plots based on regDA analysis and saves them as PDF files.
#'
#' @details
#' The function processes input files from regDA analysis, applies thresholds for significance and estimate,
#' and creates forest plots visualizing the results. It supports custom color schemes and handles multiple comparisons.
#'
#' @examples
#' # Example usage:
#' Go_regDA_fore(project = "MyProject",
#'               file_path = "path/to/data",
#'               files = "result_pattern.*\\.csv",
#'               fdr = 0.05,
#'               estimate = 1)
#'
#' @export

Go_regDA_fore <- function(project,
                         file_path, 
                         mycols=NULL, 
                         fdr=0.05, 
                         estimate=1, 
                         files, 
                         name, 
                         orders, 
                         height, 
                         width){
  if(!is.null(dev.list())) dev.off()
  # out dir
  out <- file.path(sprintf("%s_%s",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out)) dir.create(out)
  out_path <- file.path(sprintf("%s_%s/pdf",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out_path)) dir.create(out_path)

  # add input files
  path <- file_path
  filenames <- list.files(path, pattern=files);filenames
  
  print(path)
  print(filenames)
  
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
  
  # out file
  pdf(sprintf("%s/RegDA.forest.%s.%s(fdr=%s.estimate=%s).%s.pdf", out_path, 
              project, 
              ifelse(is.null(name), "", paste(name, ".", sep = "")), 
              fdr,
              estimate,
              format(Sys.Date(), "%y%m%d")), height = height, width = width)
  
  for (fn in 1:length(filenames)) {
    df <- read.csv(sprintf("%s/%s",path, filenames[fn]), row.names=NULL ,check.names=FALSE,quote = "")

    print(filenames[fn])
    parts <- strsplit(filenames[fn], "\\.")[[1]]  # Split the string by dot
    last_three_parts <- tail(parts, 3)  # Get the last two parts
    filename1 <- paste(last_three_parts, collapse = ".") # Combine them together with a dot
    filename2 <- gsub(paste0(project, "\\."), "", gsub(".csv","",filename1))
    print(filename2)
    
    
    df.sel <- df
    
    df.sel.sel <- as.data.frame(subset(df.sel, abs(df.sel$Estimate) > estimate))
    
    resSig <- as.data.frame(subset(df.sel.sel, pvalue < fdr )); resSig <- resSig[order(resSig$Estimate),]
    # resSig$smvar <- factor(resSig$smvar)
    if (dim(resSig)[1] == 0)
      next
    resSig$dir <- ifelse(resSig$pvalue < 0.05, ifelse(sign(resSig$Estimate)== 1, "up", "down"), "NS")
    resSig$dirPadj <- ifelse(resSig$padj < 0.05, TRUE, FALSE)
    
    padj_shape <- c(19,1); names(padj_shape) <- c(TRUE, FALSE)
    
    
    
    
    # 중복 이름 처리 하기
    headers <- vector(dim(resSig)[1], mode="character")
    
    for (i in 1:dim(resSig)[1]) {
      headers[i] <- paste("ASV", i, sep="_")
    }
    resSig$ASV <- headers
    
    
    for (plot in unique(resSig$mvar)){
      resSig.sel <- subset(resSig, mvar == plot)

      if (unique(resSig.sel$bas.count) > 1){
        baseline <- unique(resSig.sel$baseline)
        compare <- unique(resSig.sel$V3)
        compare <- gsub(unique(resSig.sel$mvar), "", compare);compare
      }else{
        baseline <-"Negative"
        compare <- "Positive"
      }

      
      
      # color
      resSig.sel$dir <- gsub('down',baseline, gsub('up',compare, resSig.sel$dir))
      resSig.sel$dir <- factor(resSig.sel$dir, levels = c(as.character(baseline), "NS", as.character(compare)))
      
      # color end
      
      if(!is.null(mycols)){
        dircolors <- c(mycols[1], "grey",mycols[2]); names(dircolors) <- c(as.character(baseline), "NS", as.character(compare))
      }else{
        dircolors <- c("#f8766d", "grey","#7cae00"); names(dircolors) <- c(as.character(baseline), "NS", as.character(compare))
      }
      
      if (unique(resSig.sel$bas.count) > 1){
        legend.labs <- 
          c(paste(baseline, " (n=", unique(resSig.sel$bas.count),")",sep=""),
            paste("NS"),
            paste(compare, " (n=", unique(resSig.sel$coef.count), ")",sep=""))
      }
      
      
      lims <- max(abs(resSig.sel$Estimate) + abs(resSig.sel$SE))*1.0
      p1 <- ggplot(resSig.sel, aes(x=reorder(ASV,Estimate), y=Estimate, color=dir)) + geom_point(aes(shape = dirPadj)) +  # Differentiate by shape
        geom_errorbar(aes(x=ASV, ymin=Estimate-SE, max=Estimate+SE), width=0.2) + 
        geom_hline(yintercept=0) + theme_classic()  + coord_flip() +  
        ylim(c(-lims, lims)) +
        scale_x_discrete(breaks = as.character(resSig.sel$ASV), labels = resSig.sel$Species)+
        theme(plot.title = element_text(size=8, hjust = .5)) +
        ggtitle(sprintf("regDA-%s ( pvalue < %s, cutoff=%s (closed circle fdr < 0.05)) \n %s", plot, fdr, estimate, filename2)) + labs(y = "Estimate") +labs(x = NULL)
      
      if(unique(resSig.sel$bas.count) > 1){
        p1 <- p1 + scale_color_manual(values=dircolors,labels=legend.labs,drop = FALSE) # drop = FALSE keep legend even there is no value
        p1 <- p1+ scale_shape_manual(values = padj_shape)   # Choose shapes for the points
      } else {
        p1 <- p1 + scale_color_manual(values=dircolors, drop = FALSE) # drop = FALSE keep legend even there is no value
        p1 <- p1 + scale_shape_manual(values = padj_shape)  # Choose shapes for the points
      }
      
      
      print(p1)
      
    }
  }
  dev.off()
}





