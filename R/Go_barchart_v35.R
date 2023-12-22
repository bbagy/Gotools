#' Create Bar Charts from Phyloseq Data
#'
#' This function generates bar charts from a Phyloseq object, allowing for analysis
#' of taxonomic abundance at various levels. It supports both relative and absolute
#' abundance measures and offers customization options like color, faceting, and ordering.
#'
#' @param psIN Phyloseq object containing the data for analysis.
#' @param cate.vars Categorical variables for plotting.
#' @param project Name of the project or analysis.
#' @param taxanames Taxonomic levels for aggregation and analysis.
#' @param orders Custom order of factors in the plots, if applicable.
#' @param simple Logical value to choose between simple and detailed plotting.
#' @param mycols Color palette for the bars in the chart.
#' @param relative Logical value to select between relative (TRUE) or absolute (FALSE) abundance.
#' @param y_axis Label for the Y-axis.
#' @param x_label Label for the X-axis.
#' @param facet Faceting variable for the plot.
#' @param legend Position of the legend ("bottom" or "right").
#' @param cutoff Threshold for filtering low-abundance taxa.
#' @param name Optional name for the analysis.
#' @param ncol Number of columns for facet wrapping in the plot.
#' @param height Height of the output plot.
#' @param width Width of the output plot.
#'
#' @return Saves bar chart plots as PDF files in a specified directory and outputs tables
#' containing taxa abundance data. The function supports various customization options for the plots.
#'
#' @details
#' The function preprocesses the Phyloseq object, potentially aggregating data at a specified taxonomic level,
#' and then generates bar charts for each taxonomic level specified. It handles relative and absolute abundance
#' data and allows the user to specify various aesthetic and layout options for the charts.
#'
#' @examples
#' # psIN is a Phyloseq object
#' # Example usage:
#' Go_barchart(psIN = psIN,
#'             cate.vars = c("Treatment", "Condition"),
#'             project = "MyProject",
#'             taxanames = c("Phylum", "Genus"),
#'             orders = c("Control", "Treatment"),
#'             relative = TRUE,
#'             y_axis = "Abundance",
#'             x_label = "SampleID",
#'             facet = "Treatment",
#'             legend = "bottom",
#'             cutoff = 0.005,
#'             name = "Analysis1",
#'             ncol = 2,
#'             height = 7,
#'             width = 10)
#'
#' @export

Go_barchart <- function(psIN, cate.vars, project, taxanames, orders=NULL,
                        simple = FALSE,  
                        mycols=NULL, 
                        relative = T,
                        y_axis=NULL,
                        x_label=NULL, 
                        facet=NULL, 
                        legend="bottom", 
                        cutoff=0.005, 
                        name=NULL, 
                        ncol=NULL, 
                        height, width){
    
  if(!is.null(dev.list())) dev.off()
  
  
  taxRanks <- taxanames
  
  # out dir
  out <- file.path(sprintf("%s_%s",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out)) dir.create(out)
  out_path <- file.path(sprintf("%s_%s/pdf",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out_path)) dir.create(out_path)
  out_tab <- file.path(sprintf("%s_%s/table",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out_tab)) dir.create(out_tab)
  out_taxa <- file.path(sprintf("%s_%s/table/taxa",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out_taxa)) dir.create(out_taxa)

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
    print("mycols is not defined.")
    mycols <- NULL
  }

  tt <- try(orders,T)
  if(class(tt) == "try-error"){
    print("orders is not defined.")
    orders <- NULL
  }
  

if(relative == T){
  pdf(sprintf("%s/barchart.relative.%s.%s%s(%s).%s.pdf", out_path, 
              project, 
              ifelse(is.null(facet), "", paste(facet, ".", sep = "")), 
              ifelse(is.null(name), "", paste(name, ".", sep = "")), 
              cutoff,
              format(Sys.Date(), "%y%m%d")), height = height, width = width)
}else{
  pdf(sprintf("%s/barchart.absolute.%s.%s%s(%s).%s.pdf", out_path, 
              project, 
              ifelse(is.null(facet), "", paste(facet, ".", sep = "")), 
              ifelse(is.null(name), "", paste(name, ".", sep = "")), 
              cutoff,
              format(Sys.Date(), "%y%m%d")), height = height, width = width)
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
      print("DADA2 table")
      otu.filt <- as.data.frame(t(otu_table(psIN))) 
      otu.filt[,taxanames[i]] <- getTaxonomy(otus=rownames(otu.filt), tax_tab=tax_table(psIN), taxRanks=colnames(tax_table(psIN)),level=taxanames[i])
    }else{
      otu.filt <- as.data.frame(otu_table(psIN)) 
      print("other table")
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
      write.csv(agg_other_out, quote = FALSE, col.names = NA, file=sprintf("%s/%s.taxa_relative_abundance.(%s).%s.%s%s.csv", out_taxa,
                                                                           project,cutoff,taxanames[i],
                                                                           ifelse(is.null(name), "", paste(name, ".", sep = "")), 
                                                                           format(Sys.Date(),"%y%m%d"))) #,sep="/"


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
      write.csv(agg_other_out, quote = FALSE, col.names = NA, file=sprintf("%s/%s.taxa_absolute_abundance.(%s).%s.%s%s.csv", out_taxa,
                                                                           project,cutoff,taxanames[i],
                                                                           ifelse(is.null(name), "", paste(name, ".", sep = "")), 
                                                                           format(Sys.Date(),"%y%m%d"))) #,sep="/"
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
    


    print(1)
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
    print(2)

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
        
        print(sprintf("Facet by %s-%s",mvar, facet))

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

        print(p)
      }

    }     else if (is.null(facet) & simple == FALSE) {
      for (mvar in cate.vars) {
        print("B")
        print(sprintf("Facet by %s",mvar))

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
        print(p)
      }
    } else if (is.null(facet) & simple == TRUE) {
      for (mvar in cate.vars) {

        print("C")
        print("Simpe plot")
        
        p = p
        
        if (!is.null(name)) {
          p= p+ ggtitle(sprintf("%s barplots overall of %s-%s (cut off < %s)",taxanames[i],mvar,name, cutoff))
        }
        else {
          p= p+ ggtitle(sprintf("%s barplots overall of %s (cut off < %s)",taxanames[i],mvar, cutoff))
        }
        #plotlist[[length(plotlist)+1]] <- p
        print(p)
      }
    }
  }
  dev.off()
}

