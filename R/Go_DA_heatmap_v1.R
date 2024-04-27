
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

Go_DA_heat <- function(df, project, data_type, font,
                       addnumber=TRUE, facet=NULL,
                       pval,fc, orders, name, height, width){

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

  # Assuming df is your dataframe
  if("deseq2.P" %in% colnames(df)){
    tool <- "deseq2"
  } else if("aldex2.P" %in% colnames(df)){
    tool <- "aldex2"
  } else if("ancom2.P" %in% colnames(df)){
    tool <- "ancom2"
  } else {
    tool <- NA  # or some default value or error handling
  }

  print(tool)
  if(tool == "deseq2"){
    asvs <- "taxa"
    p.value <- "pvalue"
    fdr <- "padj"
    lfc <- "log2FoldChange"
  }else if(tool == "ancom2"){
    asvs <- "ASV"
    p.value <- "pvalue_ancombc"
    fdr <- "qvalue_ancombc"
    lfc <- "lfc_ancombc"
  }else if(tool == "aldex2"){
    asvs <- "ASV"
    p.value <- "pvalue"
    fdr <- "padj"
    lfc <- "Est"
  }




  pdf(sprintf("%s/DA.heatmap.%s.%s%s(%s.%s).%s.%s.pdf", out_path,
              project,
              ifelse(is.null(facet), "", paste(facet, ".", sep = "")),
              ifelse(is.null(name), "", paste(name, ".", sep = "")),
              pval,
              fc,
              paste("(",tool,"_heatmap",")",sep=""),
              format(Sys.Date(), "%y%m%d")), height = height, width = width)


  resSig <- subset(df, df[,p.value] < pval)
  resSig <- resSig[order(resSig[[lfc]]),]
  resSig <- as.data.frame(resSig)

  print(sprintf("p < %s (n = %s) and padj < 0.05 (n = %s)",pval, length(resSig[,p.value]),dim(subset(resSig, resSig[,fdr] < 0.05))[1]))
  resSig.top <- as.data.frame(subset(resSig, abs(resSig[,lfc]) > fc))
  resSig.top <- subset(resSig.top, resSig.top[,p.value] < pval)
  resSig.top$Significance <-cut(resSig.top[,fdr], breaks=c(-Inf, 0.01, 0.05, 0.08, Inf), label=c("**", "*", ".", ""))
  resSig.top$smvar


  if (dim(resSig)[1] >= 1) {
    # re-order
    if (!is.null(orders)) {
      resSig.top$smvar <- factor(resSig.top$smvar, levels = intersect(orders, resSig.top$smvar))
      resSig.top[,facet] <- factor(resSig.top[,facet], levels = intersect(orders, resSig.top[,facet]))
    } else {
      resSig.top$smvar <- factor(resSig.top$smvar)
      resSig.top[,facet] <- factor(resSig.top[,facet])
    }

    resSig.top$basline <- paste(resSig.top$basline," (n=",resSig.top$bas.count, ")",sep="")

    # Create formatted labels
    formatted_labels <- paste(resSig.top$smvar, " (n=", resSig.top$smvar.count, ")", sep="")
    formatted_orders <- paste(orders, " (n=", resSig.top$smvar.count[match(orders, resSig.top$smvar)], ")", sep="")



    # resSig.top$smvar <- paste(resSig.top$smvar," (n=", resSig.top$smvar.count, ")",sep="")
    # re-order using number
    # new.orders <- c()
    # for(i in orders){
    #   if(length(order <- grep(i, unique(resSig.top$smvar)))){
    #     order <- c(unique(resSig.top$smvar)[order])
    #  }
    #  new.orders <- c(new.orders, order)
    # }

    new.orders <- c(formatted_orders, orders)
    # resSig.top$smvar  <- factor(resSig.top$smvar, levels = intersect(new.orders, resSig.top$smvar))
    # Convert formatted labels back into a factor with the specified order
    resSig.top$smvar <- factor(formatted_labels, levels = new.orders)





    if(tool == "deseq2"){
      p <- ggplot(resSig.top, aes(x=reorder(taxa, log2FoldChange), y=smvar, color=smvar))
    }else if(tool == "ancom2"){
      p <- ggplot(resSig.top, aes(x=reorder(ASV, lfc_ancombc), y=smvar, color=smvar))
    }else if(tool == "aldex2"){
      p <- ggplot(resSig.top, aes(x=reorder(ASV, Est), y=smvar, color=smvar))
    }





    p <- p + labs(y = "Comparison Group") + theme_classic() + coord_flip() +
      geom_tile(aes_string(fill = lfc), colour = "white") +
      scale_fill_gradient2(low = "#149BEDFF", mid = "white", high = "#FA6B09FF") +
      ggtitle(sprintf("%s baseline %s vs All group (pvalue < %s, cutoff=%s) ", unique(resSig$mvar), unique(resSig.top$basline), pval, fc)) +
      theme(plot.title = element_text(hjust = 0.5), legend.position = "right") +
      facet_wrap(~ smvar, scales="free_x", ncol = 10) +
      theme(axis.text.x = element_blank(), axis.ticks = element_blank(), text = element_text(size=12), plot.title = element_text(hjust=1),
            axis.text.y = element_text(angle=0, vjust=0.5, hjust=1, size=font,face = "italic"))


    if (data_type == "dada2" | data_type == "DADA2") {
      p1 = p + scale_x_discrete(breaks = as.character(resSig[,asvs]), labels = as.character(paste(resSig$Phylum, resSig$Species)))
    } else if (data_type == "Other" | data_type == "other") {
      p1 = p + scale_x_discrete(breaks = as.character(resSig$taxa), labels = as.character(paste(resSig$KOName)))
    }

    # Assuming groupby is always 'smvar' and there is no facet variable


    if (!is.null(facet)) {
      ncol <- length(unique(resSig.top[,facet]))*length(unique(resSig.top$smvar))
      p2 = p1 + facet_wrap(as.formula(sprintf("~ %s+%s", facet,  "smvar")), scales="free_x", ncol = ncol) +
        geom_text(aes_string(label="Significance"), color="black", size=3)
    } else {

      p2 = p1 + facet_wrap(~ smvar, scales="free_x", ncol = 10) +
        geom_text(aes_string(label="Significance"), color="black", size=3)
    }
    p3 <- p2 + theme(
      axis.text.x = element_blank(),
      axis.ticks = element_blank(),
      text = element_text(size=font),  # Ensure 'font' is defined, e.g., font <- 12
      plot.title = element_text(hjust=1),
      strip.background = element_blank()
    )



  } else{
    next
  }

  p4 <- ggplotGrob(p3)
  id <- which(p4$layout$name == "title")
  p4$layout[id, c("l","r")] <- c(1, ncol(p4))
  #grid.newpage()
  grid.draw(p4)
  dev.off()
}




