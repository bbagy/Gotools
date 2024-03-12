#' Generate Volcano Plots for Multiple Datasets
#'
#' This function creates volcano plots for a series of datasets specified by file names.
#' It supports various data analysis tools (e.g., DESeq2, ALDEx2, ANCOMBC) and includes options
#' for color customization, labeling, and output file generation.
#'
#' @param project A character string representing the project name or identifier.
#' @param file_path The path to the directory containing the data files.
#' @param files A pattern to match the filenames for processing.
#' @param fc A numeric value representing the fold-change threshold for significance in the volcano plot.
#' @param mycols A character vector specifying custom colors for plotting.
#' @param name An optional character string to specify the name or identifier for the plot.
#' @param overlaps An integer value indicating the maximum number of overlaps for text labels in the plot.
#' @param font A numeric value specifying the font size for text elements in the plot.
#' @param height The height of the output PDF file for the plot.
#' @param width The width of the output PDF file for the plot.
#'
#' @return The function generates and saves volcano plot PDF files but does not return a value.
#'
#' @examples
#' # Assuming appropriate data files are in the specified directory
#' Go_volcanoPlot(project = "MyProject", file_path = "data/", files = "deseq2_results.csv",
#'                fc = 2, mycols = c("#00BFC4", "#F8766D"), name = "ExampleVolcanoPlot",
#'                overlaps = 10, font = 12, height = 7, width = 7)
#'
#' @export
#' @importFrom ggplot2 ggplot aes_string xlab ylab geom_vline scale_color_manual theme element_text geom_text_repel ggtitle geom_point scale_shape_manual


Go_volcanoPlot <- function(project,
                       file_path,files,
                       fc,
                       mycols=NULL,
                       name,
                       overlaps=10,
                       font,
                       height, width){

  if(!is.null(dev.list())) dev.off()

  # out dir
  out <- file.path(sprintf("%s_%s",project, format(Sys.Date(), "%y%m%d")))
  if(!file_test("-d", out)) dir.create(out)
  out_path <- file.path(sprintf("%s_%s/pdf",project, format(Sys.Date(), "%y%m%d")))
  if(!file_test("-d", out_path)) dir.create(out_path)
  out_DA <- file.path(sprintf("%s_%s/pdf/DA_plot",project, format(Sys.Date(), "%y%m%d")))
  if(!file_test("-d", out_DA)) dir.create(out_DA)


  # add input files
  plot = "volcano"
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
    filename1 <- sprintf("%s/%s", path, filenames[fn]);filename1
    df <- read.csv(filename1, row.names=NULL ,check.names=FALSE,quote = "")

    df$aldex2.FDR

    # tool recognizing
    tools <- sapply(filename1, function(filename1) {
      if (grepl("deseq2", filename1)) return("deseq2")
      if (grepl("aldex2", filename1)) return("aldex2")
      if (grepl("ancom2", filename1)) return("ancom2")
      return(NA) # if none of the tools are matched
    })

    unique_tools <- unique(tools[!is.na(tools)])
    print(unique_tools)

    if (unique_tools == "deseq2") {
      x_var <- "log2FoldChange"
      y_var <- "-log10(pvalue)"
      z_var <- "log2(baseMean +1)"
      pval <- "deseq2.P"
      p <- "pvalue"
      fdr <- "deseq2.FDR"
      padj <- "padj"
      tool <- "deseq2"
      print(tool)
      model <- NULL
      vx <- "log2 fold change"
      vy <- "-log10 (p-value)"

    } else if (unique_tools == "aldex2") {
      tool <- "aldex2"
      vy <- "-log10 (p-value)"
      print(tool)
      model <- sapply(filename1, function(filename1) {
        if (grepl("t-test", filename1)) return("t-test")
        if (grepl("GLM", filename1)) return("GLM")
        if (grepl("corr", filename1)) return("corr")
        return(NA) # if none of the tools are matched
      })
      unique_model <- unique(model[!is.na(model)])

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
      print(tool)
      vx <- "log fold change"
      vy <- "-log10 (p-value)"
    }

    # Check if any filename has "confounder"
    has_confounder <- any(grepl("with_confounder", filename1));has_confounder

    if (has_confounder) {
      confounder <- "confounder"
      print(confounder)
    } else {
      confounder <- NULL
    }

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
      geom_vline(xintercept = c(-log2(fc), 0, log2(fc)), col = dircolors, linetype = "dotted", size = 1)

    p1 <- p1 + scale_color_manual(values=dircolors,  labels=legend.labs, drop = FALSE)

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
    p2 <- p1 + geom_point(aes(shape=dirPadj), size=font-1.5)+  scale_shape_manual(values = padj_shape, drop = FALSE) +
      labs(shape = "FDR < 0.05", color = sprintf("%s p < 0.05",tool)) +  #
      theme(text = element_text(size=font+8),
            plot.title = element_text(size=font+8),
            legend.text=element_text(size=font+8),
            legend.position="bottom",
            legend.justification = "left",
            legend.box.just = "left",
            legend.box = "vertical",
            legend.key = element_blank(),
            panel.grid = element_blank(),
            panel.background = element_rect(fill = "white", colour = "Black",size = 0.5, linetype = "solid"),
            aspect.ratio = 1/1.5)


    pdf(sprintf("%s/%s.%s%s.(%s.vs.%s).%s.%s%s%s(cutoff=%s).%s.pdf", out_DA,
                tool,
                ifelse(is.null(plot), "", paste(plot, ".", sep = "")),
                mvar,
                basline,
                smvar,
                project,
                ifelse(is.null(model), "", paste(model, ".", sep = "")),
                ifelse(is.null(confounder), "", paste(confounder, ".", sep = "")),
                ifelse(is.null(name), "", paste(name, ".", sep = "")),
                fc,
                format(Sys.Date(), "%y%m%d")), height = height, width = width)
    print(p2)
    dev.off()

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

      p1 <- p1 + scale_color_manual(values=dircolors,  labels=legend.labs, drop = FALSE)

      # Construct the label condition as a string using sprintf
      label_condition <- sprintf(
        "ifelse(%s != 'NA' & df.na[, '%s'] < 0.05 & abs(df.na[, '%s']) > fc, as.character(%s), '')",
        label_name, sp, a_var, label_name
      )

      # Use the constructed label_condition in your geom_text_repel function
      p1 <- p1 + geom_text_repel(aes_string(label=label_condition), size=font, fontface="italic", max.overlaps=overlaps)
      p1 <- p1 + ggtitle(sprintf("%s, %s%s (p < 0.05, cutoff=%s) ", mvar,  tool, ifelse(is.null(model), "", paste("-",model, sep = "")), fc))
      p2 <- p1 + geom_point(aes(shape=dirPadj), size=font-1.5)+  scale_shape_manual(values = padj_shape, drop = FALSE) +
        labs(shape = "FDR < 0.05", color = sprintf("%s p < 0.05",tool)) +  #
        theme(text = element_text(size=font+8),
              plot.title = element_text(size=font+8),
              legend.text=element_text(size=font+8),
              legend.position="bottom",
              legend.justification = "left",
              legend.box.just = "left",
              legend.box = "vertical",
              legend.key = element_blank(),
              panel.grid = element_blank(),
              panel.background = element_rect(fill = "white", colour = "Black",size = 0.5, linetype = "solid"),
              aspect.ratio = 1/1.5)



      pdf(sprintf("%s/%s.%s%s.(%s.vs.%s).%s.%s%s%s(cutoff=%s).%s.pdf", out_DA,
                  tool,
                  ifelse(is.null(plot), "", paste(plot, ".", sep = "")),
                  mvar,
                  basline,
                  smvar,
                  project,
                  ifelse(is.null(model1), "", paste(model1, ".", sep = "")),
                  ifelse(is.null(confounder), "", paste(confounder, ".", sep = "")),
                  ifelse(is.null(name), "", paste(name, ".", sep = "")),
                  fc,
                  format(Sys.Date(), "%y%m%d")), height = height, width = width)
      print(p2)
      dev.off()
    }
  }
}

