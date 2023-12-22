#' Generate Bar Charts for Microbiome Data
#'
#' This function creates bar charts for microbiome data analysis from a Phyloseq object. It supports various options including taxonomic aggregation, relative abundance, customization of colors, and facetting.
#'
#' @param psIN Phyloseq object containing the OTU/ASV counts and associated sample data.
#' @param cate.vars Vector of column names in the sample data of 'psIN' representing categorical variables for analysis.
#' @param project Name or identifier for the project or dataset being analyzed.
#' @param taxanames Vector of taxonomic ranks to be included in the analysis.
#' @param orders Optional ordering for the levels of the categorical variables.
#' @param simple Boolean indicating if a simplified version of the bar chart should be generated.
#' @param mycols Color palette for the bar chart.
#' @param relative Boolean indicating whether relative abundance should be used (default TRUE).
#' @param y_axis Label for the y-axis.
#' @param x_label Label for the x-axis, default is "SampleIDfactor".
#' @param facet Optional parameter for facetting the plot.
#' @param legend Position of the legend in the plot.
#' @param cutoff Cutoff value for filtering low-abundance taxa.
#' @param name Optional name for the analysis output.
#' @param ncol Number of columns for facetting.
#'
#' @return An invisible ggplot object representing the bar chart.
#'
#' @examples
#' # Assuming 'ps' is a phyloseq object with appropriate data
#' Go_barchart_markdown(ps, cate.vars = c("Group"), project = "my_project",
#'                      taxanames = c("Phylum", "Genus"), relative = TRUE)
#'
#' @export
#' @importFrom phyloseq phyloseq sample_data otu_table tax_table
#' @importFrom ggplot2 ggplot aes_string geom_bar theme_classic labs theme guides scale_fill_manual
#' @importFrom reshape2 melt
#' @importFrom gridExtra grobTree textGrob
#' @importFrom grid gpar

Go_barchart_markdown <- function(psIN, cate.vars, project, taxanames, orders=NULL,
                                 simple = FALSE,  
                                 mycols=NULL, 
                                 relative = T,
                                 y_axis=NULL,
                                 x_label=NULL, 
                                 facet=NULL, 
                                 legend="bottom", 
                                 cutoff=0.005, 
                                 name=NULL, 
                                 ncol=NULL){
  
  taxRanks <- taxanames
  
  if(!is.null(x_label)){
    x_label = x_label
  }else{
    x_label="SampleIDfactor"
  }
  
  
  
  # out file
  # "name" definition
  if (class(name) == "function"){
    name <- NULL
  }
  
  tt <- try(mycols,T)
  if(class(tt) == "try-error"){
    mycols <- NULL
  }
  
  tt <- try(orders,T)
  if(class(tt) == "try-error"){
    orders <- NULL
  }
  
  # order by bdiv
  ordi.tt <- try(ordi <- ordinate(psIN , method = "PCoA", distance = "bray"),T)
  
  if (class(ordi.tt) == "try-error"){
    map <- data.frame(sample_data(psIN))
    ordering.pc1 <- unique(map$SampleID)
  }else{
    ordi <- ordinate(psIN , method = "PCoA", distance = "bray")
    ordering.pc1 <- names(sort(ordi$vectors[,"Axis.1"]))
  }
  
  
  
  mapping.sel <- data.frame(sample_data(psIN))
  
  plotlist <- list()
  for(i in 1:length(taxanames)){
    
    # try table type
    otu.filt <- as.data.frame(otu_table(psIN)) 
    tt <- try(otu.filt[,taxanames[i]] <- getTaxonomy(otus=rownames(otu.filt), tax_tab=tax_table(psIN), taxRanks=colnames(tax_table(psIN)),level=taxanames[i]),T)
    
    if(class(tt) == "try-error"){
      otu.filt <- as.data.frame(t(otu_table(psIN))) 
      otu.filt[,taxanames[i]] <- getTaxonomy(otus=rownames(otu.filt), tax_tab=tax_table(psIN), taxRanks=colnames(tax_table(psIN)),level=taxanames[i])
    }else{
      otu.filt <- as.data.frame(otu_table(psIN)) 
      otu.filt[,taxanames[i]] <- getTaxonomy(otus=rownames(otu.filt), tax_tab=tax_table(psIN), taxRanks=colnames(tax_table(psIN)),level=taxanames[i])
    }
    
    
    
    
    #if (dim(otu.filt)[2] == 2){
    #  next
    #}
    
    agg <- aggregate(as.formula(sprintf(". ~ %s" , taxanames[i])), otu.filt, sum, na.action=na.pass)
    
    if (taxanames[i] == "Species"){
      agg <- agg[grepl("NA NA", agg$Species)==F,]
    }
    
    
    if (relative == TRUE){
      genera <- agg[,taxanames[i]]
      agg[,taxanames[i]] <- NULL
      #agg <- agg[,-1]
      agg <- normalizeByCols(agg)
      inds_to_grey <- which(rowMeans(agg) < cutoff)
      genera[inds_to_grey] <- "[1_#Other]"
      agg[,taxanames[i]] <- genera
      #saving table
      agg_other_out <- subset(agg, agg[,taxanames[i]] != "[1_#Other]")
      
      
      df <- melt(agg, variable="SampleID")
    }else if(relative == FALSE){
      genera <- agg[,taxanames[i]]
      agg[,taxanames[i]] <- NULL
      #agg <- agg[,-1]
      agg.rel <- normalizeByCols(agg)
      inds_to_grey <- which(rowMeans(agg.rel) < cutoff)
      genera[inds_to_grey] <- "[1_#Other]"
      agg[,taxanames[i]] <- genera
      #saving table
      agg_other_out <- subset(agg, agg[,taxanames[i]] != "[1_#Other]")
      df <- melt(agg, variable="SampleID")
    }
    
    
    
    
    # add StduyID
    
    
    df2 <- aggregate(as.formula(sprintf("value ~ %s + SampleID" , taxanames[i])), df, sum)
    df2$SampleID <- as.character(df2$SampleID)
    df2$SampleIDfactor <- factor(df2$SampleID, levels=ordering.pc1)
    df.SampleIDstr <- unique(df2[,c("SampleID", "SampleIDfactor")]);head(df.SampleIDstr)
    
    #mapping.sel[df2$SampleID, "StudyID"]
    
    # add groups
    for (mvar in cate.vars) {
      df.SampleIDstr$Group <- as.character(mapping.sel[df.SampleIDstr$SampleID, mvar])
      df2[,mvar] <- mapping.sel[df2$SampleID, mvar]
      
      # order
      if (length(orders) >= 1) {
        df2[,mvar] <- factor(df2[,mvar], levels = orders)
      }
      else {
        df2[,mvar] <- factor(df2[,mvar])
      }
    }
    
    # adding facet to groups
    if (!is.null(facet)) {
      for (fa in facet){
        rownames(mapping.sel) <- as.character(rownames(mapping.sel))
        df.SampleIDstr$Group <- as.character(mapping.sel[df.SampleIDstr$SampleID, fa])
        df2[,fa] <- mapping.sel[df2$SampleID, fa]
        df2[,fa] <- factor(df2[,fa], levels = orders)
      }
    }
    
    
    if (x_label == "SampleID"| x_label == "SampleIDfactor"){
      df2 <- df2
    } else if (length(x_label) >= 1){
      df2[,x_label] <- mapping.sel[df2$SampleID, x_label]
      df2[,x_label] <- factor(df2[,x_label], levels = orders)
    } 
    
    
    # color
    colourCount = length(unique(df2[,taxanames[i]]));colourCount
    
    if(!is.null(mycols)){
      getPalette = colorRampPalette(mycols)
    }else{
      p=p
    }
    
    
    
    
    
    
    # pdf size height = 5, width=9
    
    if (legend == "bottom"){
      if (colourCount < 30) {
        coln <- 4
      }else if (colourCount > 30) {
        coln <- 5
      }
    } else if (legend == "right") {
      if (colourCount <= 18) {
        coln <- 1
      } else if (colourCount > 19 & colourCount  < 35) {
        coln <- 2
      } else if (colourCount > 36) {
        coln <- 3
      }
    }
    
    # plot
    # df2 <- df2[order(df2$value, decreasing=T),]
    
    p <- ggplot(df2, aes_string(x= x_label, y="value", fill=taxanames[i], order=taxanames[i])) + 
      geom_bar(stat="identity", position="stack") + theme_classic()  + labs(fill=NULL)+
      theme(legend.position=legend, # legend.text=element_text(size=8), 
            legend.text = element_text(face = c(rep("italic", 5), rep("plain", 5))),
            axis.title.x = element_blank(), axis.text.x = element_text(angle=90, vjust=0.5, hjust=1, size=8)) + 
      guides(fill=guide_legend(ncol= coln))   #guides(col = guide_legend(ncol = coln)) + 
    
    
    if(!is.null(y_axis)){
      p <- p + labs(y = y_axis)
    }else{
      if (relative == TRUE){
        p <- p + labs(y = "Relative abundance") + ylim(c(-.1, 1.01))
      }else if(relative == FALSE){
        p <- p + labs(y = "Absolute abundance")
      }
    }
    
    
    
    
    
    
    if(!is.null(mycols)){
      p=p + scale_fill_manual(values = getPalette(colourCount)) 
    }else{
      p=p
    }
    
    
    
    if (!is.null(facet)) {
      for (mvar in cate.vars) {        
        if (facet == mvar) {
          next
        }
        
        df2[,facet] <- factor(df2[,facet], levels = orders)
        
        if (!is.null(ncol)) {
          p <- p+ facet_wrap(as.formula(sprintf("~ %s + %s", paste(setdiff(facet, "SampleType"), collapse="+"), mvar)), scales = "free_x", ncol = ncol) 
        }else{
          p <- p+ facet_grid(as.formula(sprintf("~ %s + %s", paste(setdiff(facet, "SampleType"), collapse="+"), mvar)), scales = "free_x", space = "free") 
        }
        
        
        
        
        
        if (!is.null(name)) {
          p = p+ ggtitle(sprintf("Taxa barplots overall of %s-%s (cut off < %s)",mvar,name, cutoff))
        }
        else {
          p= p+ ggtitle(sprintf("Taxa barplots overall of %s (cut off < %s)",mvar, cutoff))
        }
        
        plot(p)
      }
      
    }     else if (is.null(facet) & simple == FALSE) {
      for (mvar in cate.vars) {
        
        if (!is.null(ncol)) {
          p <- p + facet_wrap(as.formula(sprintf("~ %s"  ,mvar)), scales = "free_x", ncol = ncol)  
        }else{
          p <- p + facet_grid(as.formula(sprintf("~ %s"  ,mvar)), scales = "free_x", space = "free_x") 
        }
        
        
        if (length(name) == 1) {
          p= p+ ggtitle(sprintf("%s barplots overall of %s-%s (cut off < %s)",taxanames[i],mvar,name, cutoff))
        }
        else {
          p= p+ ggtitle(sprintf("%s barplots overall of %s (cut off < %s)",taxanames[i], mvar, cutoff))
        }
        #plotlist[[length(plotlist)+1]] <- p
        plot(p)
      }
    } else if (is.null(facet) & simple == TRUE) {
      for (mvar in cate.vars) {
        
        
        
        p = p
        
        if (!is.null(name)) {
          p= p+ ggtitle(sprintf("%s barplots overall of %s-%s (cut off < %s)",taxanames[i],mvar,name, cutoff))
        }
        else {
          p= p+ ggtitle(sprintf("%s barplots overall of %s (cut off < %s)",taxanames[i],mvar, cutoff))
        }
        #plotlist[[length(plotlist)+1]] <- p
        plot(p)
      }
    }
  }
  invisible(p)
}
