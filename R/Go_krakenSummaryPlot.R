
#' Create Kraken Summary Plot
#'
#' This function creates a bar plot summarizing the results of a Kraken analysis.
#'
#' @param project A string indicating the project name for file naming.
#' @param logTab The path to the log file containing Kraken analysis results.
#' @param name (Optional) A string to add as a suffix to the file name.
#' @param with_human A logical value indicating whether to include human reads in the plot (default is TRUE).
#' @param legend A string indicating the position of the legend ("right" or other valid ggplot2 positions).
#' @param height The height of the output plot.
#' @param width The width of the output plot.
#'
#' @details
#' The function reads Kraken analysis results and prepares a bar plot showing the relative abundance of various
#' categories like Archaea, Bacteria, Fungi, Viruses, and optionally human reads. The plot is saved as a PDF.
#'
#' @return
#' A bar plot saved as a PDF file in the specified output directory.
#'
#' @examples
#' Go_krakenSummaryPlot(project = "MyProject",
#'                      logTab = "path/to/kraken_log.csv",
#'                      name = "Analysis1",
#'                      with_human = TRUE,
#'                      legend = "right",
#'                      height = 6,
#'                      width = 8)
#'
#' @export


Go_krakenSummaryPlot <- function(project, logTab,
                                 name=NULL,
                                 with_human = T,
                                 legend="right",
                                 height, width){

  if(!is.null(dev.list())) dev.off()
  logTab <- read.csv(logTab,row.names=1,check.names=F);head(logTab)
  
  logTab.t <-  t(logTab)
  if(with_human == T){
    get.names <- c("Human_reads", "ST_unclassified", "d__Archaea", "d__Bacteria", "d__Eukaryota|k__Fungi", "d__Viruses")
  }else{
    get.names <- c("ST_unclassified", "d__Archaea", "d__Bacteria", "d__Eukaryota|k__Fungi", "d__Viruses")
  }
  
  
  logTab.t.sel <- subset(logTab.t, rownames(logTab.t) %in% get.names);logTab.t.sel
  
  
  
  # Melt the data
  df_melted <- melt(logTab.t.sel)
  
  # Rename the variables
  colnames(df_melted) <- c("Category", "Sample", "Count")
  
  df_melted <- df_melted %>%
    group_by(Sample) %>%
    mutate(Percentage = Count / sum(Count) * 100)
  
  
  logTab.t <-  t(logTab)
  if(with_human == T){
    mycols <- c("#5FB233FF", "#6A7F93FF", "#F57206FF", "#EB0F13FF", "#8F2F8BFF", "#1396DBFF")
  }else{
    mycols <- c("#6A7F93FF", "#F57206FF", "#EB0F13FF", "#8F2F8BFF", "#1396DBFF")
  }
  
  
  
  
  getPalette = colorRampPalette(mycols)
  colourCount = length(unique(df_melted$Category));colourCount
  
  # Bar chart
  p <- ggplot(df_melted, aes(x = Sample, y = Percentage, fill = Category)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = paste0(round(Percentage, 1), "%")), 
              position = position_stack(vjust = 0.5), size = 2.5) +
    theme_minimal() + scale_fill_manual(values = getPalette(colourCount))+
    labs(x = "Sample", y = "Percentage (%)", fill = "Category") +
    theme(legend.position=legend, 
          legend.title = element_blank(), 
          axis.text.x = element_text(angle = 45, hjust = 1)) + ggtitle(sprintf("Summary of %s project",project))
  

  
  out <- file.path(sprintf("%s_%s",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out)) dir.create(out)
  out_path <- file.path(sprintf("%s_%s/pdf",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out_path)) dir.create(out_path)
  
  
  
  pdf(sprintf("%s/kraken2.summary.%s.%s%s%s.pdf", out_path, 
              project, 
              ifelse(with_human, "withHuman", "withoutHuman"),
              ifelse(is.null(name), "", paste(name, ".", sep = "")), 
              format(Sys.Date(), "%y%m%d")), height = height, width = width)
  
  
  
  print(p)
  
  dev.off()
  
} 

