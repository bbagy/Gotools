
#' Create a Differential Abundance Heatmap
#'
#' @param df Data frame containing differential abundance results.
#' @param project Name of the project.
#' @param data_type Type of data used, e.g., 'dada2' or 'other'.
#' @param facet Name of the facet column in the data frame.
#' @param groupby Name of the grouping column in the data frame.
#' @param font Font size for the plot text.
#' @param addnumber Logical indicating whether to add numbers to the plot.
#' @param fdr False Discovery Rate threshold for filtering significant features.
#' @param fc Fold Change threshold for filtering significant features.
#' @param orders Vector of ordered factor levels for arranging groups.
#' @param name Optional name for the saved plot.
#' @param height Height of the output plot.
#' @param width Width of the output plot.
#'
#' @details
#' This function generates a heatmap to visualize differential abundance results. It highlights significant features based on FDR and fold change thresholds and arranges them according to specified orders. The function is designed for comparative analysis in microbiome studies.
#'
#' @return
#' A heatmap plot saved as a PDF file, showing differentially abundant features.
#'
#' @examples
#' Go_DA_heat(df = diff_abundance_results,
#'            project = "MyProject",
#'            data_type = "dada2",
#'            facet = "Condition",
#'            groupby = "Group",
#'            font = 12,
#'            addnumber = TRUE,
#'            fdr = 0.05,
#'            fc = 2,
#'            orders = c("Group1", "Group2"),
#'            name = "HeatmapPlot",
#'            height = 10,
#'            width = 15)
#'
#' @export

Go_DA_heat <- function(df, project, data_type, facet,groupby,font,
                       addnumber=TRUE,
                       fdr,fc, orders, name, height, width){
    
  if(!is.null(dev.list())) dev.off()
   
  # out dir
  out <- file.path(sprintf("%s_%s",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out)) dir.create(out)
  out_path <- file.path(sprintf("%s_%s/pdf",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out_path)) dir.create(out_path)
  
  # out file
  # "name" definition
  if (class(name) == "function"){
    name <- NULL
  }
  pdf(sprintf("%s/DA.heatmap.%s.%s%s(%s.%s).%s.pdf", out_path, 
              project, 
              ifelse(is.null(facet), "", paste(facet, ".", sep = "")), 
              ifelse(is.null(name), "", paste(name, ".", sep = "")), 
              fdr, 
              fc, 
              format(Sys.Date(), "%y%m%d")), height = height, width = width)
  
  
  resSig <- as.data.frame(subset(df, padj < fdr)); resSig <- resSig[order(resSig$log2FoldChange),]
  
  
  if (length(subset(resSig, ancom == TRUE)) > 1){
    print(sprintf("Combination Deseq2(%s) and Ancom(%s)",length(resSig$deseq2),length(subset(resSig, ancom == TRUE))))
    resSig.top <- as.data.frame(subset(resSig, abs(resSig$log2FoldChange) > fc))
    resSig.top <- subset(resSig.top, ancom == TRUE)
    
  } else{
    print(sprintf("Use only Deseq2(%s)",length(resSig.top$deseq2)))
    resSig.top <- as.data.frame(subset(resSig, abs(resSig$log2FoldChange) > fc))
  }
  
  #print("c")
  #if (length(unique(resSig$smvar)) >=2 ){
  
  if (dim(resSig)[1] >= 1) {
    # re-order
    if (!is.null(orders)) {
      resSig.top[,groupby] <- factor(resSig.top[,groupby], levels = intersect(orders, resSig.top[,groupby]))
      resSig.top[,facet] <- factor(resSig.top[,facet], levels = intersect(orders, resSig.top[,facet]))
      print("Re-ordered")
    } else {
      resSig.top[,groupby] <- factor(resSig.top[,groupby])
      resSig.top[,facet] <- factor(resSig.top[,facet])
      print(2)
    }
    


    resSig.top$basline <- paste(resSig.top$basline," (n=",resSig.top$bas.count, ")",sep="")
    resSig.top$smvar <- paste(resSig.top$smvar," (n=", resSig.top$smvar.count, ")",sep="")
    
    
    # re-order using number
    new.orders <- c()
    for(i in orders){
      if(length(order <- grep(i, unique(resSig.top$smvar)))){
        order <- c(unique(resSig.top$smvar)[order])
      }
      new.orders <- c(new.orders, order)
    }
    
    tt <- try(resSig.top$smvar  <- factor(resSig.top[,facet], levels = intersect(new.orders, resSig.top$smvar)),T)
    
    if (class(tt) =="try-error"){
      resSig.top$smvar  <- factor(resSig.top[,groupby], levels = intersect(new.orders, resSig.top$smvar))
    } else{
      resSig.top$smvar  <- factor(resSig.top[,facet], levels = intersect(new.orders, resSig.top$smvar))
    }
    


    
    

    print(1)
    if (groupby == "smvar"){
      p <- ggplot(resSig.top, aes(x=reorder(taxa,log2FoldChange), y=smvar, color=smvar)) + theme_classic()+ coord_flip() #x=reorder(taxa,Estimate); 원래 x=factor(taxa). 값에 따라 정열 하기 위해x=reorder(taxa,Estimate)를 사용함
 
    }  else {
      p <- ggplot(resSig.top, aes(x=reorder(taxa,log2FoldChange), y=mvar, color=mvar)) + theme_classic()+ coord_flip()#x=reorder(taxa,Estimate); 원래 x=factor(taxa). 값에 따라 정열 하기 위해x=reorder(taxa,Estimate)를 사용함
    }

    
    p = p + geom_tile(aes(fill = log2FoldChange), colour = "white") + 
      labs(y = "Comparison Group") +labs(x = NULL) +
      scale_fill_gradient2(low = "#1170aa", mid = "white", high = "#fc7d0b")+
      ggtitle(sprintf("%s baseline %s vs %s (padj < %s, cutoff=%s) ", unique(resSig$mvar), unique(resSig$basline), "All groups",  fdr,fc))  + 
      theme(plot.title = element_text(hjust = 0.5),legend.position= "right")+ #0.5
      theme(axis.text.x = element_text(angle=0, vjust=0.5, hjust=1, size=8),
             axis.text.y = element_text(angle=0, vjust=0.5, hjust=1, size=8,face = "italic"),
            plot.title=element_text(size=9,face="bold")) 
    
    
    print(2)
    if (data_type == "dada2" | data_type == "DADA2") {
      p1 = p + scale_x_discrete(breaks = as.character(resSig$taxa), labels = as.character(paste(resSig$Phylum, resSig$ShortName)))
    } else if (data_type == "Other" | data_type == "other") {
      p1 = p + scale_x_discrete(breaks = as.character(resSig$taxa), labels = as.character(paste(resSig$KOName)))
    }
    
    print(3)
    if (groupby == "smvar"){
      if (length(facet) == 1) {
        ncol <- length(unique(resSig.top[,facet]))*length(unique(resSig.top[,"smvar"]))
        p2 = p1 + facet_wrap(as.formula(sprintf("~ %s+%s", "smvar", facet)), scales="free_x", ncol = ncol)
      } else {
        p2 = p1 + facet_wrap(~  smvar, scales="free_x", ncol = 10)
      }
    }else if (groupby == "des"){
      if (length(facet) == 1) {
        ncol <- length(unique(resSig.top[,facet]))*length(unique(resSig.top[,"des"]))
        p2 = p1 + facet_wrap(as.formula(sprintf("~ %s+%s", "des", facet)), scales="free_x", ncol = ncol)
      } else {
        p2 = p1 + facet_wrap(~  des, scales="free_x", ncol = 10)
      }
    }
    #print(4)
    #plotlist[[length(plotlist)+1]] <- p
    p3 = p2 + theme(axis.text.x = element_blank(), axis.ticks = element_blank()) + theme(text = element_text(size=font), plot.title = element_text(hjust=1))
   # print(p3)
  }else{
    next
  }

  
  p4 <- ggplotGrob(p3)
  id <- which(p4$layout$name == "title")
  p4$layout[id, c("l","r")] <- c(1, ncol(p4))
  #grid.newpage()
  grid.draw(p4)
  dev.off()
}

