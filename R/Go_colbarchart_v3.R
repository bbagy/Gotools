
#' Generate Color Bar Charts for Microbiome Data
#'
#' Produces color bar charts representing the relative or absolute abundance of microbial taxa, categorized by specified variables.
#'
#' @param psIN A phyloseq object containing microbiome data.
#' @param cate.vars A vector of categorical variables to be used for grouping the data in the bar charts.
#' @param project A string representing the project name, used for file naming.
#' @param taxanames A vector of taxonomic ranks to be visualized in the charts.
#' @param data_type A string indicating the type of data, affects taxonomy naming.
#' @param orders A vector specifying the order of categories within the plot facets.
#' @param simple A logical value, if TRUE, produces simpler plots without facets.
#' @param mycols Optional, custom color palette for the plot.
#' @param relative A logical value indicating whether to plot relative abundances.
#' @param x_label A string for the x-axis label, default is "SampleIDfactor".
#' @param facet A string or vector of variables to create facets in the plot.
#' @param legend Position of the legend, either "bottom" or "right".
#' @param cutoff A numeric value for filtering out low-abundance taxa.
#' @param name Optional, a string for additional naming in the plot title.
#' @param ncol Optional, number of columns in the facet layout.
#' @param height Plot height.
#' @param width Plot width.
#'
#' @details
#' The function creates bar charts showing the abundance of microbial taxa, grouped by specified variables. It can handle both relative and absolute abundance data, and allows for flexible customization of plot appearance.
#'
#' @return
#' Creates and saves bar chart plots in PDF format in the specified output directory.
#'
#' @examples
#' # Example usage:
#' Go_colbarchart(psIN = my_phyloseq_object,
#'                cate.vars = c("Treatment", "TimePoint"),
#'                project = "MyMicrobiomeStudy",
#'                taxanames = c("Phylum", "Genus"),
#'                data_type = "16S",
#'                orders = c("Control", "Treatment"),
#'                relative = TRUE,
#'                height = 6,
#'                width = 8)
#'
#' @export


Go_colbarchart <- function(psIN, cate.vars, project, taxanames, data_type, orders, 
                        simple = FALSE,  
                        mycols=NULL, 
                        relative = T,
                        x_label=NULL, 
                        facet=NULL, 
                        legend="bottom", 
                        cutoff=0.005, 
                        name=NULL, 
                        ncol=NULL,
                        height, width){
    if(!is.null(dev.list())) dev.off()
    
  # out dir
  out <- file.path(sprintf("%s_%s",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out)) dir.create(out)
  out_path <- file.path(sprintf("%s_%s/pdf",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out_path)) dir.create(out_path)

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
    print("orders is not defined.")
    mycols <- NULL
  }
  
  # logic for out file
if(relative == T){
  pdf(sprintf("%s/colbarchart.relative.%s.%s%s(%s).%s.pdf", out_path, 
              project, 
              ifelse(is.null(facet), "", paste(facet, ".", sep = "")), 
              ifelse(is.null(name), "", paste(name, ".", sep = "")), 
              cutoff,
              format(Sys.Date(), "%y%m%d")), height = height, width = width)
}else{
  pdf(sprintf("%s/colbarchart.absolute.%s.%s%s(%s).%s.pdf", out_path, 
              project, 
              ifelse(is.null(facet), "", paste(facet, ".", sep = "")), 
              ifelse(is.null(name), "", paste(name, ".", sep = "")), 
              cutoff,
              format(Sys.Date(), "%y%m%d")), height = height, width = width)
}

  taxRanks <- taxanames
  
  # order by bdiv
  ordi <- ordinate(psIN , method = "PCoA", distance = "bray")
  ordering.pc1 <- names(sort(ordi$vectors[,"Axis.1"]))
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
    

    # continue
    otu.filt[,taxanames[i]] <- getTaxonomy(otus=rownames(otu.filt), tax_tab=tax_table(psIN), taxRanks=taxRanks,level=taxanames[i])
    otu.filt$PhylumCol <- getTaxonomy(otus=rownames(otu.filt), tax_tab=tax_table(psIN), taxRanks=taxRanks, level="Phylum")

    if (dim(otu.filt)[2] == 2){
      next
    }

    agg <- aggregate(as.formula(sprintf(". ~ %s + PhylumCol" , taxanames[i])), otu.filt, sum, na.action=na.pass)
    genera <- agg[,taxanames[i]]
    PhylumCol <- agg$PhylumCol
    agg[,taxanames[i]] <- NULL
    agg$PhylumCol <- NULL

    agg <- normalizeByCols(agg)
    inds_to_grey <- which(rowMeans(agg) < cutoff)
    genera[inds_to_grey] <- "[1_#Other]"
    agg[,taxanames[i]] <- genera
    agg$PhylumCol <- PhylumCol 
    
    
    
    if (taxanames[i] == "Phylum"){
      agg$Phylum <- genera
    }
    
    df <- melt(agg, variable.name="SampleID")


    # add StduyID

    df2 <- aggregate(as.formula(sprintf("value ~ %s + PhylumCol + SampleID" , taxanames[i])), df, sum)
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
    if (length(facet) == 1) {
      for (fa in facet){
        df.SampleIDstr$Group <- as.character(mapping.sel[df.SampleIDstr$SampleID, fa])
        df2[,fa] <- mapping.sel[df2$SampleID, fa]
      }
    }
    
    


    if (x_label == "SampleID"| x_label == "SampleIDfactor"){
      df2 <- df2
    } else if (length(x_label) >= 1){
      df2[,x_label] <- mapping.sel[df2$SampleID, x_label]
      df2[,x_label] <- factor(df2[,x_label], levels = orders)
    }  


    print(1)
    #------------------------#
    # ---  Color table   --- #
    #------------------------#
    agg$PhylumCol <- PhylumCol 
    agg[,taxanames[i]] <- genera
    

    #-------- remove other from taxa table --------#
    TaxaTab <- agg[order(agg[,taxanames[i]] ,  decreasing = TRUE), ]
    cdf <- data.frame(subset(TaxaTab, select=c("PhylumCol", taxanames[i])))
    cdf.sel <- subset(cdf, cdf[,taxanames[i]] != "[1_#Other]");dim(cdf.sel)[1]
    
    # 몇개인지 결정후 Phylum 으로 정리
    N <- dim(cdf.sel)[1]
    cdf.sel <- cdf.sel[order(cdf.sel$PhylumCol ,  decreasing = FALSE), ]
    cdf.sel <- data.frame(as.character(cdf.sel$PhylumCol[1:N]), as.character(cdf.sel[,taxanames[i]][1:N]))
    colnames(cdf.sel) <- c("PhylumCol", taxanames[i])
    #cdf.sel[ ,c("Kingdom","Class", "Order", "Family","Genus")] <- list(NULL)
    
    cdf.sel[,taxanames[i]] <-  gsub("p__", "", gsub("c__", "", gsub("o__", "", gsub("f__", "", gsub("g__", "", gsub("s__", "", cdf.sel[,taxanames[i]]))))))
    
    # save species name
    taxaName <- cdf.sel[,taxanames[i]]
    cdf.sel[,taxanames[i]] <- NULL
    
    # -----  create color table   ---- #
    coltab <- Go_color(cdf=cdf.sel, taxaName=taxaName)
    
    # hsv code
    #print(coltab$color_table)
    #coltab$color_table$Phylum
    # color code
    #print(coltab$coloring)
    
    print(2)
    # pdf size height = 5, width=9
    if (legend == "bottom"){
      if (N < 30) {
        col <- 5
      }
    } else if (legend == "right") {
      if (N < 18) {
        col <- 1
      }
      else if (N > 19 & N  < 35) {
        col <- 2
      }
      else if (N > 36) {
        col <- 3
      }
    }

    # plot
    # df2 <- df2[order(df2$value, decreasing=T),]
    print(3)
    level <- unique(df2[,taxanames[i]])
    #facet <- "SampleType"
    #mvar <- "TreatmentGroup"
    df2[,facet] <- factor(df2[,facet], levels = orders)
    if (length(facet) == 1) {
      for (mvar in cate.vars) {
        if (class(ncol) == "numeric") {
          ncol <- ncol
        }else if(length(unique(df2[,mvar])) >= 1){
          ncol <- length(unique(df2[,mvar]))*length(unique(df2[,facet]))
        }
        if (facet == mvar) {
          next
        }
        
        df2[,facet] <- factor(df2[,facet], levels = orders)
        print(4)
        
        p <- ggplot(df2, aes_string(x= x_label, y="value", fill=factor(df2[,taxanames[i]], levels=level), order=taxanames[i])) + geom_bar(stat="identity", position="stack") + theme_classic()  + theme(legend.position=legend, legend.text=element_text(size=8), axis.title.x = element_blank(), axis.text.x = element_text(angle=90, vjust=0.5, hjust=1, size=8)) + guides(fill=guide_legend(ncol= col))  + guides(col = guide_legend(ncol = col)) + ylim(c(-.1, 1.01)) + scale_fill_manual(values=coltab$coloring) + facet_wrap(as.formula(sprintf("~ %s + %s", paste(setdiff(facet, "SampleType"), collapse="+"), mvar)), scales="free_x", ncol = ncol) + labs(y = "Relative abundance") + labs(fill = taxanames[i])

        
        if (length(name) == 1) {
          p = p+ ggtitle(sprintf("Taxa barplots overall of %s-%s (cut off < %s)",mvar,name, cutoff))
        }
        else {
          p= p+ ggtitle(sprintf("Taxa barplots overall of %s (cut off < %s)",mvar, cutoff))
        }

        print(p)
        #plotlist[[length(plotlist)+1]] <- p
      }

    } else if (length(facet) != "NULL") {
      for (mvar in cate.vars) {
        if (class(ncol) == "numeric") {
          ncol <- ncol
        }else if(length(unique(df2[,mvar])) >= 1){
          ncol <- length(unique(df2[,mvar]))
        }

        p <- ggplot(df2, aes_string(x= x_label, y="value", fill=factor(df2[,taxanames[i]], levels=level), order=taxanames[i])) + 
        geom_bar(stat="identity", position="stack") + theme_classic()  + 
        theme(legend.position= legend, legend.text=element_text(size=8), axis.title.x = element_blank(), 
        axis.text.x = element_text(angle=90, vjust=0.5, hjust=1, size=8)) + guides(fill=guide_legend(ncol= col)) + 
        guides(col = guide_legend(ncol = col)) + ylim(c(-.1, 1.01)) + scale_fill_manual(values=coltab$coloring) + 
        facet_wrap(as.formula(sprintf("~ %s"  ,mvar)), scales="free_x", ncol = ncol) + labs(y = "Relative abundance")+ labs(fill = taxanames[i])
        if (length(name) == 1) {
          p= p+ ggtitle(sprintf("Taxa barplots overall of %s-%s (cut off < %s)",mvar,name, cutoff))
        }
        else {
          p= p+ ggtitle(sprintf("Taxa barplots overall of %s (cut off < %s)",mvar, cutoff))
        }
        #plotlist[[length(plotlist)+1]] <- p
        print(p)
      }
    }
  }
  #multiplot(plotlist=plotlist, cols=plotCols, rows=plotRows)
  dev.off()
}
