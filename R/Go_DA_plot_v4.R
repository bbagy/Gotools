
#' Generate Differential Abundance (DA) Plots
#'
#' This function generates various plots (volcano, MA, forest) for visualizing differential abundance (DA) data.
#'
#' @param project Name of the project or analysis.
#' @param file_path Directory path containing input files for analysis.
#' @param files Pattern to match filenames for analysis.
#' @param plot Type of plot to generate: "volcano", "maplot", or "forest".
#' @param fdr False Discovery Rate threshold for significance.
#' @param fc Fold change threshold for including in the plot.
#' @param mycols Optional color palette for the plot.
#' @param name Optional name for the analysis.
#' @param overlaps Maximum number of overlaps allowed in text repelling.
#' @param font Font size for plot text.
#' @param height Height of the output plot.
#' @param width Width of the output plot.
#'
#' @return Generates DA plots based on specified parameters and saves them as PDF files.
#'
#' @details
#' The function processes input files containing DA analysis results and generates specified types of plots. It supports custom color schemes and handles multiple comparisons.
#'
#' @examples
#' # Example usage:
#' Go_DA_plot(project = "MyProject",
#'            file_path = "path/to/data",
#'            files = "result_pattern.*\\.csv",
#'            plot = "volcano",
#'            fdr = 0.05,
#'            fc = 1)
#'
#' @export

Go_DA_plot <- function(project, file_path,files, plot = "volcano", fdr, fc, mycols=NULL, name, overlaps=10, font, height, width){
    
  if(!is.null(dev.list())) dev.off()
   
  # out dir
  out <- file.path(sprintf("%s_%s",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out)) dir.create(out)
  out_path <- file.path(sprintf("%s_%s/pdf",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out_path)) dir.create(out_path)
  out_DA <- file.path(sprintf("%s_%s/pdf/DA_plot",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out_DA)) dir.create(out_DA)

  
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
    
    # get data tyep
    print("Check the data type")
    taxtab.col <- colnames(df)
    
    if (any(grepl("Species", taxtab.col))){
      type <- "taxonomy"
      print(type)
    }else if(any(grepl("KO", taxtab.col))){
      type <- "function"
      print(type)
    }else if(any(grepl("pathway", taxtab.col))){
      type <- "function"
      print(type)
    }else if(any(grepl("symbol", taxtab.col))){
      type <- "RNAseq"
      print(type)
    }
    
    
    
    # remove NA
    df[df==""] <- "NA"
    #df$deseq2 <- ifelse(abs(df$log2FoldChange) > fc, ifelse(sign(df$log2FoldChange)==1, "up", "down"), "NS")
    df$deseq2 <- ifelse(df$padj < fdr & abs(df$log2FoldChange) > fc, ifelse(sign(df$log2FoldChange)==1, "up", "down"), "NS")
    df.na <- df[!is.na(df$deseq2), ]

    basline <- unique(df.na$basline)
    smvar <- unique(df.na$smvar)
    mvar <- unique(df.na$mvar)


    # colors and names
    df.na$deseq2<- gsub('down',basline, gsub('up',smvar, df.na$deseq2))
    
    df.na$deseq2 <- factor(df.na$deseq2, levels = c(as.character(basline), "NS", as.character(smvar)))
    

   if(!is.null(mycols)){
    dircolors <- c(mycols[1], "grey",mycols[2]); names(dircolors) <- c(as.character(basline), "NS", as.character(smvar))
    }else{
     dircolors <- c("#f8766d", "grey","#7cae00"); names(dircolors) <- c(as.character(basline), "NS", as.character(smvar))
    }
    
    legend.labs <- 
      c(paste(basline, " (n=", unique(df.na$bas.count),")",sep=""),
        paste("NS"),
        paste(smvar, " (n=", unique(df.na$smvar.count), ")",sep=""))

    names(legend.labs) <-  c(as.character(basline), "NS", as.character(smvar))

    
    
    if (type == "taxonomy") {
      df.na$diff_abn[is.na(df.na$diff_abn)] <- FALSE
      diff_abnshape <- c(16,1); names(diff_abnshape) <- c(TRUE, FALSE)# 16,1,1
    }
    

    

    #-----------------------------#
    #   plot style for Deseq2     #
    #-----------------------------#
    if(plot == "volcano"){
      print("Generating Volcano plots.")
      p1 <- ggplot(data=df.na, aes(x=log2FoldChange, y=-log10(pvalue),colour=deseq2)) + 
        xlab("log2 fold change") + ylab("-log10 (p-value)")+ 
        geom_vline(xintercept = c(-log2(fc), 0,log2(fc)),col = dircolors, linetype = "dotted", size = 1)
      
    } else if(plot == "maplot"){
      print("Generating M (log ratio) A (mean average)  plots.")
      p1 <-  ggplot(df.na, aes(x=log2(baseMean +1), y=log2FoldChange, colour=deseq2)) +
        xlab("Log2 mean of normalized counts") + ylab("Log2 fold change")+ 
        geom_hline(yintercept = c(-log2(fc), 0,log2(fc)),col = dircolors, linetype = "dotted", size = 1)
      
    } else if(plot == "forest"){
      print("Generating forest plots.")
      resSig <- as.data.frame(subset(df.na, padj < fdr)); resSig <- resSig[order(resSig$log2FoldChange),]
      resSig.top <- as.data.frame(subset(resSig, abs(resSig$log2FoldChange) > fc))
      if (dim(resSig)[1] == 0 | dim(resSig.top)[1] == 0 ){
        next
      }
      
      resSig$smvar <- factor(resSig$smvar)
      lims <- max(abs(resSig$log2FoldChange) + abs(resSig$lfcSE))*1.0
      resSig.top$deseq2<- gsub('down',basline, gsub('up',smvar, resSig.top$deseq2))
      resSig.top$deseq2 <- factor(resSig.top$deseq2, levels = c(as.character(basline), "NS", as.character(smvar)))
      
      if (type == "taxonomy") {
        resSig.top$diff_abn[is.na(resSig.top$diff_abn)] <- FALSE
      }
      


      p1 <- ggplot(resSig.top, aes(x=reorder(taxa,log2FoldChange), y=log2FoldChange, color=deseq2)) + 
        geom_hline(yintercept=0) + coord_flip() + #theme_classic() + 
        scale_color_manual(values=dircolors, labels=legend.labs, drop = FALSE)  +  # drop = FALSE keep legend even there is no value
        geom_errorbar(aes(x=taxa, ymin=log2FoldChange-lfcSE, max=log2FoldChange+lfcSE), width=0.2)  + 
        ylim(c(-lims, lims))+ xlab("Taxa") + ylab("log2FoldChange")+
        theme(text = element_text(size=font), plot.title = element_text(size=font, hjust = 1),
              axis.text.y = element_text(angle=0, vjust=0.5, hjust=1, size=font,face = "italic")) + #hjust =1
        theme(panel.grid = element_blank(),
              panel.background = element_rect(fill = "white", colour = "Black",size = 0.5, linetype = "solid"), 
              aspect.ratio = 1/0.7)
      
      
      if(type == "taxonomy" | type == "taxanomy"){
        p1 <- p1 + scale_x_discrete(breaks = as.character(resSig$taxa), labels = sprintf("%s__%s__%s",as.character(resSig$Phylum),as.character(resSig$Family), as.character(resSig$ShortName))) 
      } else if(type == "bacmet" ){
        p1 <- p1 + scale_x_discrete(breaks = as.character(resSig$taxa), labels = sprintf("%s",as.character(resSig$ShortName))) 
      }
    }
    
    
    if(plot == "volcano" |  plot == "maplot"){
      p1 <- p1 + scale_color_manual(values=dircolors, labels=legend.labs, drop = FALSE) + 
        theme(text = element_text(size=font+8), 
              plot.title = element_text(size=font+8), 
              legend.text=element_text(size=font+8),  
              legend.position="bottom",
              legend.justification = "left",
              legend.box = "vertical")+
        theme(panel.grid = element_blank(),
              panel.background = element_rect(fill = "white", colour = "Black",size = 0.5, linetype = "solid"), 
              aspect.ratio = 1/1.5) 

      # + theme_bw()
      
      if(type == "taxonomy" | type == "taxanomy" |type == "bacmet" ){
        p1 <- p1 +  geom_text_repel(aes(label=ifelse(ShortName != "NA" & df.na$padj < fdr & abs(df.na$log2FoldChange) > fc, as.character(ShortName),'')), size=font, segment.fdr = 0.25, fontface="italic",max.overlaps = overlaps )
      }else if(type == "function"){
        p1 <- p1 +  geom_text_repel(aes(label=ifelse(KOName != "NA" & df.na$padj < fdr & abs(df.na$log2FoldChange) > fc, as.character(KOName),'')), size=font,max.overlaps = overlaps)
      }else if(type == "RNAseq"){
        p1 <- p1 +  geom_text_repel(aes(label=ifelse(symbol != "NA" & df.na$padj < fdr & abs(df.na$log2FoldChange) > fc, as.character(symbol),'')), size=font,max.overlaps = overlaps)
      }
    }

    
    if(!is.null(df.na$des)){
      des <- unique(df.na$des)
      p1 <- p1 + ggtitle(sprintf("%s-%s, (padj < %s,cutoff=%s) ", mvar, des, fdr, fc))
    }else{
      p1 <- p1 + ggtitle(sprintf("%s, (padj < %s,cutoff=%s) ", mvar,  fdr, fc))
    }

    
    if (type == "taxonomy") {
      p1 <- p1 + geom_point(aes(shape=diff_abn), size=font-1.5) + scale_shape_manual(values = diff_abnshape) 
    }else{
      p1 <- p1 + geom_point()  #     geom_point(aes(shape=deseq2), size=font-1.5) 
    }
    

    
    pdf(sprintf("%s/%s%s.(%s.vs.%s).%s.%s(%s.%s).%s.pdf", out_DA, 
                ifelse(is.null(plot), "", paste(plot, ".", sep = "")), 
                mvar,
                basline, 
                smvar,
                project, 
                ifelse(is.null(name), "", paste(name, ".", sep = "")), 
                fdr, 
                fc, 
                format(Sys.Date(), "%y%m%d")), height = height, width = width)
    print(p1)
    dev.off()
    
    
    #-----------------------------#
    #          ANCOM plot         #
    #-----------------------------#
    #===== ANCOM plot 
    if (type == "taxonomy") {
      out_AC <- file.path(sprintf("%s_%s/pdf/AC_plot",project, format(Sys.Date(), "%y%m%d"))) 
      if(!file_test("-d", out_AC)) dir.create(out_AC)
      
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
      
      df.sel <- subset(df.na, diff_abn == T);dim(df.sel)
      if(dim(df.sel)[1] == 0){
        next
      }
      
      
      df.sel$ancom <- as.factor(df.sel$ancom)
      
      
      levels(df.sel$ancom)
      
      p = df.sel %>%
        ggplot(aes(x=reorder(taxa,lfc_ancombc), y = lfc_ancombc, fill = ancom)) + 
        geom_bar(stat = "identity", width = 0.55, color = "black", size = 0.25,
                 position = position_dodge(width = 0.4)) +geom_hline(yintercept=0, size=0.3) +
        geom_errorbar(aes(ymin = lfc_ancombc - se_ancombc, ymax = lfc_ancombc + se_ancombc), 
                      width = 0.2, position = position_dodge(0.05), color = "black") + 
        labs(x = NULL, y = "Log fold change", 
             title = "ANCOM Log fold changes") + 
        scale_fill_manual(values=dircolors, labels=legend.labs, drop = FALSE) +
        scale_color_discrete(name = NULL)+
        theme(text = element_text(size=font*1.5), 
              plot.title = element_text(size=font*1.5, hjust = 0.5),  #hjust = 1
              legend.key = element_rect(fill = "transparent"),
              legend.key.size = unit(0.3, 'cm'), # legend size
              panel.grid.minor.y = element_blank(),
              legend.title= element_blank(),
              axis.text.x = element_text(angle = 60, hjust = 1),
              axis.text.y = element_text(angle=0, vjust=0.5, hjust=1, size=font*1.5,face = "italic")) 
      
      
      
      # plot size ratio
      if (length(df.sel$ancom) < 4){
        num.subgroup <- 1
      }else{
        num.subgroup <- 1 + length((df.sel$ancom))*0.1
      }
      
      
      if (num.subgroup == 1 ){
        height.ac = 2
        width.ac = 2
      }else{
        height.ac = 2+length((df.sel$ancom))*0.005
        width.ac = 2+length((df.sel$ancom))*0.01
      }
      
      
      #print(height)
      #print(width)
      
      p1 <- p + scale_x_discrete(breaks = as.character(df.sel$taxa), labels = as.character(df.sel$ShortName)) + 
        theme(legend.position = c(0.85, 0.15))+
        theme(panel.grid = element_blank(),
              panel.background = element_rect(fill = "white", colour = "Black",size = 0.5, linetype = "solid"), 
              aspect.ratio = 1/num.subgroup)
      
      
      
      pdf(sprintf("%s/%s%s.(%s.vs.%s).%s.%s.%s.pdf", out_AC, 
                  ifelse(is.null(plot), "", paste("DA_ancom", ".", sep = "")), 
                  mvar,
                  basline, 
                  smvar,
                  project, 
                  ifelse(is.null(name), "", paste(name, ".", sep = "")), 
                  format(Sys.Date(), "%y%m%d")), height = height.ac, width = width.ac)
      print(p1)
      dev.off()
    }
  } 
  print("plot for volcano, maplot and forest")
}

