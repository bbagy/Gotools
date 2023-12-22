
#' Generate Alluvial Plots for Microbial Community Analysis
#'
#' This function creates alluvial plots to visualize the associations between microbial species and
#' various experimental groups or conditions.
#'
#' @param project A string indicating the project name, used for directory and file naming.
#' @param SigASVs Either a file path to a CSV containing significant ASVs or a data frame.
#' @param map A file path to a CSV containing sample metadata or a data frame.
#' @param targets.bac A vector of target bacterial species for the plot.
#' @param Addcolumn A vector of additional columns from the metadata to be included in the plot.
#' @param outcome The primary outcome or variable of interest.
#' @param column1 (Optional) An additional column for stratifying the data.
#' @param column2 (Optional) Another column for further stratification.
#' @param orders The order of factors to arrange the plot.
#' @param name (Optional) A name for the plot for labeling purposes.
#' @param height Height of the output plot.
#' @param width Width of the output plot.
#' @param plotCols Number of columns in the plot layout.
#' @param plotRows Number of rows in the plot layout.
#'
#' @details
#' The function merges microbial data with sample metadata to create informative alluvial plots. It
#' allows for customizing the stratification and visualization of data based on specified microbial
#' targets and metadata columns.
#'
#' @return
#' Saves alluvial plots as PDF files in the specified project directory.
#'
#' @examples
#' # Example usage:
#' Go_alluvialplot(project = "MyProject",
#'                 SigASVs = "path_to_sig_ASVs.csv",
#'                 map = "path_to_metadata.csv",
#'                 targets.bac = c("Bacteroides_fragilis", "Escherichia_coli"),
#'                 Addcolumn = c("Age", "Gender"),
#'                 outcome = "TreatmentResponse",
#'                 column1 = "Diet",
#'                 column2 = "BMI",
#'                 orders = c("Control", "Treatment"),
#'                 name = "AlluvialAnalysis",
#'                 height = 2, width = 3,
#'                 plotCols = 2, plotRows = 1)
#'
#' @export


Go_alluvialplot <- function(project, 
                            SigASVs, 
                            map,
                            targets.bac,
                            Addcolumn,
                            outcome,
                            column1 = NULL,
                            column2 = NULL,
                            orders,
                            name = NULL,
                            height = 2, width=3,
                            plotCols=2, plotRows=1){
  
  if(!is.null(dev.list())) dev.off()
  
  #===== out dir
  out <- file.path(sprintf("%s_%s",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out)) dir.create(out)
  out_path <- file.path(sprintf("%s_%s/pdf",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out_path)) dir.create(out_path)


  #===== input
  if (class(SigASVs) == "character"){
    tab <- read.csv(SigASVs,row.names=1,check.names=FALSE)
  }else{
    tab <- SigASVs
  }
  
  
  if (class(map) == "character"){
    sampledata <- read.csv(map,row.names=1,check.names=FALSE)
  }else{
    sampledata <- map
  }
  

  
  
  #===== check input
  if(!"names" %in% colnames(tab)) {
    stop(paste("Column", "name", "does not exist in the dataframe. (name_read counts)"))
  }
  
  if(!outcome %in% Addcolumn) {
    stop(paste("outcome", outcome, "does not exist in the data. (name_read counts)"))
  }
  
  
  
  rownames(tab) <- tab$names
  rownames(tab) <- gsub(" ","_", rownames(tab))
  tab$Species <- gsub(" ","_", tab$Species)
  

  
  plotlist <- list()
  for(target in targets.bac){

    print(target)
    
    #===== Colors
    if (target == "Gardnerella_vaginalis" ){
      mycols <- c("#CBD588", "#5F7FC7", "orange",  "#DA5724", "#508578", "#CD9BCD", "#AD6F3B", "#673770", "#D14285", "#652926", "#C84248", "#8569D5", "#5E738F", "#D1A33D", "#8A7C64", "#599861")
    }else{
      mycols <- c("#EB5291FF", "#1794CEFF", "#972C8DFF", "#FBBB68FF", "#F5BACFFF", "#9DDAF5FF", "#6351A0FF", "#ECF1F4FF", "#FEF79EFF" )#pony
    }
    
    #===== Merge tab + sampledata
    tab.sel <- subset(tab, Species == target)
    
    taxaTab <- merge(sampledata, t(tab.sel), by="row.names");head(taxaTab)
    
    rownames(taxaTab) <- taxaTab$Row.names
    
    
    #===== remove unused ranks
    for (rank in c("Kingdom","Phylum","Class","Order","Family","Genus","Species","Sum","names")){
      tab.sel[,rank] <- NULL
    }
    
  
    
    taxaTab.log <- taxaTab
    for(species in rownames(tab.sel)){
      taxaTab[,species] <- as.numeric(taxaTab[,species])
      taxaTab.log[,species] <- log(taxaTab[,species])
    }
    
    
    # melting tab
    bacteria.taxa <- rownames(tab.sel)
    
    taxaTab.melt <- melt(taxaTab.log, id.vars =  Addcolumn, measure.vars = bacteria.taxa)
    

    # Check for missing values
    
    if (any(is.na(taxaTab.melt))) {
      print("Data contains missing values.")
    } else {
      print("Data doesn't contain any missing values.")
    } #is_alluvia_form(as.data.frame(taxaTab.melt), axes = 1:3, silent = TRUE)
    
    
    
    
    taxaTab.melt <- arrange(taxaTab.melt, taxaTab.melt$variable)
    
    # Simplify all names
    taxaTab.melt$variable <- gsub("(\\w)\\w+_(\\w+)", "\\1.\\2", taxaTab.melt$variable)
    
    
    for(vari in Addcolumn){
      taxaTab.melt[,vari] <- as.character(taxaTab.melt[,vari] )
    }
    

    
    #===== Define column1 and column2 variables

    if (!is.null(column1) && !is.null(column2)) {
      # both column1 and column2 are provided
      if (!(column1 %in% colnames(taxaTab.melt) & column2 %in% colnames(taxaTab.melt))) {
        # neither column1 nor column2 is in the dataframe
        stop("Neither 'column1' nor 'column2' exists in the dataframe.")
      }else{
        print(paste("column1 is", ifelse(is.null(column1), "NULL.", paste("'", column1, "'", "and exists in the dataframe:", column1 %in% colnames(taxaTab.melt)))))
        print(paste("column2 is", ifelse(is.null(column2), "NULL.", paste("'", column2, "'", "and exists in the dataframe:", column2 %in% colnames(taxaTab.melt)))))
      }
    } else{
      print(paste("column1 is", ifelse(is.null(column1), "NULL.", paste("'", column1, "'", "and exists in the dataframe:", column1 %in% colnames(taxaTab.melt)))))
      print(paste("column2 is", ifelse(is.null(column2), "NULL.", paste("'", column2, "'", "and exists in the dataframe:", column2 %in% colnames(taxaTab.melt)))))
      
    }
    

    
    # conntol is.infinite
    df <- taxaTab.melt
    df$value[is.infinite(df$value)] <- NA
    df <- df[!is.na(df$value), ]
  

    df <- df %>% 
      group_by(variable) %>%
      mutate(count = n())
    
    df <- df %>% 
      mutate(label = ifelse(variable %in% unique(variable), paste(variable, " (n=", count, ")", sep=""), variable))
    
    
    
    #===== Logic for the vatiation
    if (!is.null(column1) & !is.null(column2)){
      p <- ggplot(data = df, aes(y = value, axis1 = label, axis2 = !!sym(column1), axis3 = !!sym(column2), axis4 = !!sym(outcome)))
    } else if (!is.null(column1) | !is.null(column2)){
      column = ifelse(!is.null(column1), column1, column2)
      p <- ggplot(data = df, aes(y = value, axis1 = label, axis2 = !!sym(column), axis3 = !!sym(outcome)))
    } else {
      p <- ggplot(data = df, aes(y = value, axis1 = label, axis2 = !!sym(outcome)))
    }
  
    
    # Modify the geom_text call in the ggplot function to include the counts
    p1 <- p +
      theme_void() +
      theme(legend.position="none",
            axis.title.y = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            plot.title = element_text(hjust = 0.1, size = 10),
            text = element_text(size = 8)) +
      geom_alluvium(aes(fill = variable), width = 5/12) +
      geom_stratum(width = 5/12, fill = "aliceblue", color = "black") +
      geom_text(stat = "stratum", aes(label = after_stat(stratum)), size=3) +
      scale_x_discrete(limits = c("CST", "Virus_CINB2"), expand = c(.2, .2)) +
      scale_fill_manual(values = mycols) +
      ggtitle(paste("Contribution of ",target))
    
    
    plotlist[[length(plotlist)+1]] <- p1 
  }
  pdf(sprintf("%s/%s.Alluvial.plots.%s.(%s%s)%s%s.pdf", out_path, 
              project,
              outcome, 
              ifelse(is.null(column1), "", paste(column1, ".", sep = "")), 
              ifelse(is.null(column2), "", paste(column2, ".", sep = "")), 
              ifelse(is.null(name), "", paste(name, ".", sep = "")), 
              format(Sys.Date(), "%y%m%d")), height = height, width = width)
  
  multiplot(plotlist=plotlist, cols=plotCols, rows=plotRows)
  dev.off()
}
