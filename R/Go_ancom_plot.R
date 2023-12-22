
#' Generate ANCOM Log Fold Change Plots
#'
#' This function creates bar plots showing log fold changes from ANCOM results for each metadata variable comparison.
#'
#' @param project A string representing the name of the project.
#' @param file_path A string indicating the path to the input files.
#' @param files A pattern string to identify the specific CSV files to be processed.
#' @param type The type of analysis data, default is "taxonomy".
#' @param mycols Optional, a vector of colors for plot aesthetics.
#' @param name An optional string to define a specific naming convention for output files.
#' @param overlaps Integer, the maximum number of labels to allow overlapping.
#' @param font The font size for plot text.
#' @param height The height of the output plot.
#' @param width The width of the output plot.
#'
#' @return
#' Generates PDF files containing ANCOM log fold change plots for each metadata variable comparison.
#'
#' @details
#' The function reads in CSV files containing ANCOM results, processes the data to highlight significant differences, and then generates bar plots showing the log fold changes for the comparisons of interest. These plots are saved as PDF files.
#'
#' @examples
#' # Example usage:
#' Go_ancom_plot(project = "MyProject", file_path = "data/", files = "ancom_results*.csv",
#'               type = "taxonomy", mycols = c("red", "blue"), name = "experiment1",
#'               overlaps = 10, font = 12, height = 6, width = 8)
#'
#' @export


Go_ancom_plot <- function(project, file_path,files, type="taxonomy", mycols=NULL, name, overlaps=10, font, height, width){
  
  if(!is.null(dev.list())) dev.off()
  
  # out dir
  out <- file.path(sprintf("%s_%s",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out)) dir.create(out)
  out_path <- file.path(sprintf("%s_%s/pdf",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out_path)) dir.create(out_path)
  out_AC <- file.path(sprintf("%s_%s/pdf/AC_plot",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out_AC)) dir.create(out_AC)
  
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
  
  
  for (fn in 1:length(filenames)) {
    
    df <- read.csv(sprintf("%s/%s",path, filenames[fn]), row.names=NULL ,check.names=FALSE,quote = "")
    # remove NA
    df[df==""] <- "NA"
    #df$deseq2 <- ifelse(abs(df$log2FoldChange) > fc, ifelse(sign(df$log2FoldChange)==1, "up", "down"), "NS")
    
    
    df$ancom <- ifelse(df$diff_abn == T, ifelse(sign(df$lfc_ancombc)==1, "up", "down"), "NS")
    
    
    df.na <- df[!is.na(df$ancom), ]
    df.na$ancom
    basline <- unique(df$basline)
    smvar <- unique(df$smvar)
    mvar <- unique(df$mvar)
    
    
    # colors and names
    df.na$ancom <- gsub('down',basline, gsub('up',smvar, df.na$ancom))
    df.na$ancom <- factor(df.na$ancom, levels = c(as.character(basline), "NS", as.character(smvar)))
    
    
    if(!is.null(mycols)){
      dircolors <- c(mycols[1], mycols[2]); names(dircolors) <- c(as.character(basline), as.character(smvar))
    }else{
      dircolors <- c("#f8766d", "#7cae00"); names(dircolors) <- c(as.character(basline), as.character(smvar))
    }
    
    legend.labs <- 
      c(paste(basline, " (n=", unique(df.na$bas.count),")",sep=""),
        paste(smvar, " (n=", unique(df.na$smvar.count), ")",sep=""))
    
    names(legend.labs) <-  c(as.character(basline), as.character(smvar))
    
    df.sel <- subset(df.na, diff_abn == T);df.sel
    
    df.sel$ancom <- as.factor(df.sel$ancom)
    
    
    levels(df.sel$ancom)
     
    p = df.sel %>%
      ggplot(aes(x=reorder(taxa,lfc_ancombc), y = lfc_ancombc, fill = ancom)) + 
      geom_bar(stat = "identity", width = 0.7, color = "black", 
               position = position_dodge(width = 0.4)) +
      geom_errorbar(aes(ymin = lfc_ancombc - se_ancombc, ymax = lfc_ancombc + se_ancombc), 
                    width = 0.2, position = position_dodge(0.05), color = "black") + 
      labs(x = NULL, y = "Log fold change", 
           title = "Log fold changes as one unit increase of age") + 
      scale_fill_manual(values=dircolors, labels=legend.labs, drop = FALSE) +
      scale_color_discrete(name = NULL)+
    theme_bw() + 
      theme(text = element_text(size=font), 
            plot.title = element_text(size=font, hjust = 0.5),  #hjust = 1
            legend.key = element_rect(fill = "transparent"),
            panel.grid.minor.y = element_blank(),
            axis.text.x = element_text(angle = 60, hjust = 1),
            axis.text.y = element_text(angle=0, vjust=0.5, hjust=1, size=font,face = "italic")) 
    
    
    
    p1<- p + scale_x_discrete(breaks = as.character(df.sel$taxa), labels = as.character(df.sel$ShortName)) + theme(legend.position = c(0.85, 0.15))

    
    pdf(sprintf("%s/%s%s.(%s.vs.%s).%s.%s.%s.pdf", out_AC, 
                ifelse(is.null(plot), "", paste("DA_ancom", ".", sep = "")), 
                mvar,
                basline, 
                smvar,
                project, 
                ifelse(is.null(name), "", paste(name, ".", sep = "")), 
                format(Sys.Date(), "%y%m%d")), height = height, width = width)
    print(p1)
    dev.off()
  }
}





