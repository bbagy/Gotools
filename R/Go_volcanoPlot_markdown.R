#' Create Volcano Plots from Data Tables
#'
#' This function generates volcano plots for each data frame in a list of data frames.
#' It supports different types of data analysis tools (e.g., DESeq2, ALDEx2, ANCOMBC) and
#' provides options for color customization and text display.
#'
#' @param DAtab A list of data frames, each representing a dataset for which a volcano plot will be generated.
#' @param fc A numeric value representing the fold-change threshold for significance in the volcano plot.
#' @param mycols A character vector specifying custom colors for plotting.
#' @param name An optional character string to specify the name or identifier for the plot.
#' @param overlaps An integer value indicating the maximum number of overlaps for text labels in the plot.
#' @param font A numeric value specifying the font size for text elements in the plot.
#'
#' @return The function generates and prints volcano plots but does not return a value.
#'
#' @examples
#' # Assuming 'data_list' is a list of data frames suitable for volcano plotting
#' Go_volcanoPlot_markdown(DAtab = data_list, fc = 2, mycols = c("#00BFC4", "#F8766D"), name = "ExamplePlot", overlaps = 10, font = 12)
#'
#' @export
#' @importFrom ggplot2 ggplot aes_string xlab ylab geom_vline scale_color_manual theme element_text geom_text_repel ggtitle geom_point scale_shape_manual

Go_volcanoPlot_markdown <- function(DAtab,
                                    fc,
                                    mycols=NULL, 
                                    name, 
                                    overlaps=10, 
                                    font){

    
  
  # add input files
  plot = "volcano"

  # out file
  # "name" definition
  if (class(name) == "function"){
    name <- NULL
  }
  
  tt <- try(mycols,T)
  if(class(tt) == "try-error"){
    mycols <- NULL
  }

  for (fn in 1:length(DAtab)) {
      df <- DAtab[[fn]]
    
    df$aldex2.FDR

    # tool recognizing
    get_tool <- function(df) {
      if (any(sapply(colnames(df), function(col) grepl("deseq2", col)))) return("deseq2")
      if (any(sapply(colnames(df), function(col) grepl("aldex2", col)))) return("aldex2")
      if (any(sapply(colnames(df), function(col) grepl("ancom2", col)))) return("ancom2")
      return(NA) # if none of the tools are matched
    }
    
    unique_tools <- get_tool(df)

    
    if (unique_tools == "deseq2") {
      x_var <- "log2FoldChange"
      y_var <- "-log10(pvalue)"
      z_var <- "log2(baseMean +1)"
      pval <- "deseq2.P"
      p <- "pvalue" 
      fdr <- "deseq2.FDR"
      padj <- "padj"
      tool <- "deseq2"
      model <- NULL
      vx <- "log2 fold change"
      vy <- "-log10 (p-value)"
      
    } else if (unique_tools == "aldex2") {
      tool <- "aldex2"
      vy <- "-log10 (p-value)"

      get_tool <- function(df) {
        if (any(sapply(df$model, function(col) grepl("t-test", col)))) return("t-test")
        if (any(sapply(df$model, function(col) grepl("GLM", col)))) return("GLM")
        if (any(sapply(df$model, function(col) grepl("corr", col)))) return("corr")
        return(NA) # if none of the tools are matched
      }
      
      unique_model <- get_tool(df)
      
      
      if (unique_model=="t-test"){
        x_var <- "diff.btw"
        y_var <- "-log10(wi.ep)"
        z_var <- "diff.bwt"
        pval <- "aldex2.P"
        p <- "wi.ep"
        fdr <- "aldex2.FDR"
        padj <- "wi.eBH"
        model <- "t-test"
        vx <- "Difference Between Groups"
        
      }else if (unique_model=="GLM"){
        est_col_name <- sprintf("%s%s.Est",unique(df$mvar), unique(df$smvar))
        est_col_name <- gsub("\\+", ".", est_col_name)
        
        pvalue_col_name <- sprintf("%s%s.pval", unique(df$mvar), unique(df$smvar))
        pvalue_col_name <- gsub("\\+", ".", pvalue_col_name)
        
        holm_col_name <- sprintf("%s%s.pval.holm",unique(df$mvar), unique(df$smvar))
        holm_col_name <- gsub("\\+", ".", holm_col_name)
      
        x_var <- est_col_name
        y_var <- sprintf("-log10(%s)", pvalue_col_name)
        z_var <-sprintf("-log10(%s)", pvalue_col_name)
        p <- pvalue_col_name
        pval <- "aldex2.P"
        fdr <- "aldex2.FDR"
        padj <- holm_col_name
        model <- "GLM"
        vx <- "Effect size"
      }else if (unique_model=="corr"){
        x_var <- "kendall.etau"
        y_var <- "-log10(kendall.ep)"
        pval <- "aldex2.kP"
        fdr <- "aldex2.kFDR"
        p <- "kendall.ep"
        padj <- "kendall.eBH"
        
        vx <- "correlation"
        vy <- "-log10 (p-value)"
        
        model <- "corr-kendall"
        
        a_var <- "spearman.erho"
        b_var <- "-log10(spearman.ep)" 
        spearman <- "aldex2.sP"
        sp <- "spearman.ep"
        spadj <- "spearman.eBH"
        
        model1 <- "corr-spearman"
      }
    } else if (unique_tools == "ancom2") {
      x_var <- "lfc_ancombc"
      y_var <- "-log10(pvalue_ancombc)"
      pval <- "ancom2.P"
      p <- "pvalue_ancombc"
      fdr <- "ancom2.FDR"
      padj <- "qvalue_ancombc"
      tool <- "ancom2"
      model <- NULL
      vx <- "log fold change"
      vy <- "-log10 (p-value)"
    }
    
    # Check if any filename has "confounder"
    # has_confounder <- any(grepl("with_confounder", filename1));has_confounder

    
    # if (has_confounder) {
    #  confounder <- "confounder"
    #  print(confounder)
    #} else {
    #  confounder <- NULL
    #}
    
    # get data tyep
    taxtab.col <- colnames(df)
    
    if (any(grepl("Species", taxtab.col))){
      type <- "taxonomy"
    }else if(any(grepl("KO", taxtab.col))){
      type <- "function"
    }else if(any(grepl("pathway", taxtab.col))){
      type <- "function"
    }else if(any(grepl("symbol", taxtab.col))){
      type <- "RNAseq"
    }
    
    # Clean the dataframe
    df[df == ""] <- "NA"
    df.na <- df[!is.na(df[, pval]), ]
    

    
    if (model == "t-test" || model == "GLM" || is.null(model)) {
      # Get unique values for some columns
      basline <- unique(df.na$basline)
      smvar <- unique(df.na$smvar)
      mvar <- unique(df.na$mvar)
      
      # Transform p-value column
      df.na[, pval] <- gsub('down', basline, gsub('up', smvar, df.na[, pval]))
      df.na[, pval] <- factor(df.na[, pval], levels = c(as.character(basline), "NS", as.character(smvar)))
      
      # Set color schemes
      if (!is.null(mycols)) {
        dircolors <- c(mycols[1], "grey", mycols[2])
      } else {
        dircolors <- c("#f8766d", "grey", "#7cae00")
      }
      names(dircolors) <- c(as.character(basline), "NS", as.character(smvar))
      
      # Set legend labels
      legend.labs <- c(
        paste(basline, " (n=", unique(df.na$bas.count), ")", sep = ""),
        "NS",
        paste(smvar, " (n=", unique(df.na$smvar.count), ")", sep = "")
      )
      names(legend.labs) <- c(as.character(basline), "NS", as.character(smvar))
    }else if(model == "corr-kendall") {
      # Get unique values for some columns
      basline <- "Negative"
      smvar <- "Positive"
      mvar <- unique(df.na$mvar)
      
      # Transform p-value column
      df.na[, pval] <- gsub('down', basline, gsub('up', smvar, df.na[, pval]))
      df.na[, pval] <- factor(df.na[, pval], levels = c(as.character(basline), "NS", as.character(smvar)))
      
      # Set color schemes
      if (!is.null(mycols)) {
        dircolors <- c(mycols[1], "grey", mycols[2])
      } else {
        dircolors <- c("#f8766d", "grey", "#7cae00")
      }
      names(dircolors) <- c(as.character(basline), "NS", as.character(smvar))
      
      # Set legend labels
      legend.labs <- c(basline, "NS", smvar)
      names(legend.labs) <- c(as.character(basline), "NS", as.character(smvar))
      
    }

    
    # Adjust FDR column
    df.na$dirPadj <- ifelse(df.na[, padj] < 0.05, TRUE, FALSE)
    padj_shape <- c(19, 1)
    names(padj_shape) <- c(TRUE, FALSE)

    #-----------------------------#
    #   Volcano plot     #
    #-----------------------------#
    p1 <- ggplot(data=df.na, aes_string(x=x_var, y=y_var, colour=pval)) + 
      xlab(vx) + 
      ylab(vy) + 
      geom_vline(xintercept = c(-log2(1), 0, log2(fc)), col = dircolors, linetype = "dotted", size = 1)
    
    p1 <- p1 + scale_color_manual(values=dircolors,  labels=legend.labs, drop = FALSE) +  # 
      theme(text = element_text(size=font+8), 
            plot.title = element_text(size=font+8), 
            legend.text=element_text(size=font+8),  
            legend.position="bottom",
            legend.justification = "left",
            legend.box = "vertical")+
      theme(panel.grid = element_blank(),
            panel.background = element_rect(fill = "white", colour = "Black",size = 0.5, linetype = "solid"), 
            aspect.ratio = 1/1.5) 
    
    if (type == "taxonomy" | type == "taxanomy" | type == "bacmet") {
      label_name <- "Species"
    } else if (type == "function") {
      label_name <- "KOName"
    } else if (type == "RNAseq") {
      label_name <- "symbol"
    }
    
    # Construct the label condition as a string using sprintf
    label_condition <- sprintf(
      "ifelse(%s != 'NA' & df.na[, '%s'] < 0.05 & abs(df.na[, '%s']) > fc, as.character(%s), '')", 
      label_name, p, x_var, label_name
    )
    
    # Use the constructed label_condition in your geom_text_repel function
    p1 <- p1 + geom_text_repel(aes_string(label=label_condition), size=font, fontface="italic", max.overlaps=overlaps)
    p1 <- p1 + ggtitle(sprintf("%s, %s%s (p < 0.05, cutoff=%s) ", mvar,  tool, ifelse(is.null(model), "", paste("-",model, sep = "")), fc))
    p1 <- p1 + geom_point(aes(shape=dirPadj), size=font-1.5)+  scale_shape_manual(values = padj_shape, drop = FALSE) 

    
    print(p1)
    
    if (model == "t-test" || model == "GLM" || is.null(model)) {
      next
    } else if(model == "corr-kendall"){
      # Get unique values for some columns
      basline <- "Negative"
      smvar <- "Positive"
      mvar <- unique(df.na$mvar)
      
      # Transform p-value column
      df.na[, spearman] <- gsub('down', basline, gsub('up', smvar, df.na[, spearman]))
      df.na[, spearman] <- factor(df.na[, spearman], levels = c(as.character(basline), "NS", as.character(smvar)))
      
      # Set color schemes
      if (!is.null(mycols)) {
        dircolors <- c(mycols[1], "grey", mycols[2])
      } else {
        dircolors <- c("#f8766d", "grey", "#7cae00")
      }
      names(dircolors) <- c(as.character(basline), "NS", as.character(smvar))
      
      # Set legend labels
      legend.labs <- c(basline, "NS", smvar)
      names(legend.labs) <- c(as.character(basline), "NS", as.character(smvar))
      # Adjust FDR column
      df.na$dirPadj <- ifelse(df.na[, spadj] < 0.05, TRUE, FALSE)
      padj_shape <- c(19, 1)
      names(padj_shape) <- c(TRUE, FALSE)
      
      p1 <- ggplot(data=df.na, aes_string(x=a_var, y=b_var, colour=spearman)) + 
        xlab(vx) + 
        ylab(vy) + 
        geom_vline(xintercept = c(-log2(fc), 0, log2(fc)), col = dircolors, linetype = "dotted", size = 1)
      
      p1 <- p1 + scale_color_manual(values=dircolors,  labels=legend.labs, drop = FALSE) +  # 
        theme(text = element_text(size=font+8), 
              plot.title = element_text(size=font+8), 
              legend.text=element_text(size=font+8),  
              legend.position="bottom",
              legend.justification = "left",
              legend.box = "vertical")+
        theme(panel.grid = element_blank(),
              panel.background = element_rect(fill = "white", colour = "Black",size = 0.5, linetype = "solid"), 
              aspect.ratio = 1/1.5) 
      
      # Construct the label condition as a string using sprintf
      label_condition <- sprintf(
        "ifelse(%s != 'NA' & df.na[, '%s'] < 0.05 & abs(df.na[, '%s']) > fc, as.character(%s), '')", 
        label_name, sp, a_var, label_name
      )
      
      # Use the constructed label_condition in your geom_text_repel function
      p1 <- p1 + geom_text_repel(aes_string(label=label_condition), size=font, fontface="italic", max.overlaps=overlaps)
      p1 <- p1 + ggtitle(sprintf("%s, %s%s (p < 0.05, cutoff=%s) ", mvar,  tool, ifelse(is.null(model), "", paste("-",model, sep = "")), fc))
      p1 <- p1 + geom_point(aes(shape=dirPadj), size=font-1.5)+  scale_shape_manual(values = padj_shape, drop = FALSE) 
      
      print(p1)

    }
  } 
}

