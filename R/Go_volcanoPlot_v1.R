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
#' @param patchwork Logical. If \code{TRUE}, skip saving and return the plot object(s) for use with \code{Gg_patchwork()} or the \pkg{patchwork} package. Default \code{FALSE}.
#' @export
#' @importFrom ggplot2 ggplot aes_string xlab ylab geom_vline scale_color_manual theme element_text ggtitle geom_point scale_shape_manual
#' @importFrom ggrepel geom_text_repel


Go_volcanoPlot <- function(project,
                       file_path = NULL,
                       files = NULL,
                       result = NULL,
                       fc,
                       mycols=NULL,
                       name = NULL,
                       overlaps=10,
                       font,
                       height, width,
                       patchwork = FALSE){

  build_taxonomy_plot_label <- function(df) {
    if (nrow(df) == 0) {
      return(character(0))
    }

    pick_first_valid <- function(candidates) {
      vals <- as.character(candidates)
      vals <- vals[!is.na(vals) & nzchar(trimws(vals)) & !vals %in% c("NA", "NA NA", "__", "s__", "g__", "f__", "o__", "c__", "p__", "k__")]
      if (length(vals) == 0) {
        return(NA_character_)
      }
      vals[1]
    }

    taxonomy_priority <- intersect(c("Species", "Genus", "Family", "Order", "Class", "Phylum", "TaxaName"), colnames(df))
    asv_col <- if ("ASV" %in% colnames(df)) {
      "ASV"
    } else if ("Row.names" %in% colnames(df)) {
      "Row.names"
    } else {
      NULL
    }

    labels <- vapply(seq_len(nrow(df)), function(i) {
      vals <- vapply(taxonomy_priority, function(col) df[i, col], character(1))
      lbl <- pick_first_valid(vals)
      if (is.na(lbl) && !is.null(asv_col)) {
        lbl <- as.character(df[i, asv_col])
      }
      if (is.na(lbl) || !nzchar(trimws(lbl))) {
        lbl <- paste0("Feature_", i)
      }
      lbl
    }, character(1))

    labels
  }

  build_optional_metric_note <- function(df) {
    first_metric <- function(candidates) {
      hits <- candidates[candidates %in% colnames(df)]
      if (length(hits) == 0) {
        return(NA_real_)
      }
      for (col in hits) {
        vals <- suppressWarnings(as.numeric(df[[col]]))
        vals <- vals[is.finite(vals) & !is.na(vals)]
        if (length(vals) > 0) {
          return(vals[1])
        }
      }
      NA_real_
    }

    auroc <- first_metric(c("auroc", "AUROC", "mean_auroc", "roc", "ROC"))
    auprc <- first_metric(c("auprc", "AUPRC", "mean_auprc", "prc", "PRC"))

    parts <- c(
      if (is.finite(auroc)) sprintf("ROC=%.3f", auroc) else NULL,
      if (is.finite(auprc)) sprintf("PRC=%.3f", auprc) else NULL
    )
    if (length(parts) == 0) {
      return(NULL)
    }
    paste(parts, collapse = " | ")
  }

  if(!is.null(dev.list())) dev.off()

  # out dir
  out <- file.path(sprintf("%s_%s",project, format(Sys.Date(), "%y%m%d")))
  if(!dir.exists(out)) dir.create(out)
  out_path <- file.path(sprintf("%s_%s/pdf",project, format(Sys.Date(), "%y%m%d")))
  if(!dir.exists(out_path)) dir.create(out_path)
  # build file list
  plot = "volcano"
  if (!is.null(result)) {
    # result is a directory path returned by Go_Deseq2 / Go_Aldex2 / Go_Ancom2
    if (!is.character(result) || !dir.exists(result)) {
      stop("result must be a directory path returned by Go_Deseq2 / Go_Aldex2 / Go_Ancom2.")
    }
    tool_files <- list.files(result, pattern = "\\.csv$", full.names = FALSE)
    if (length(tool_files) == 0) stop(sprintf("No CSV files found in: %s", result))
    file_list <- data.frame(path = result, file = tool_files, stringsAsFactors = FALSE)

  } else if (!is.null(file_path) && !is.null(files)) {
    filenames <- list.files(file_path, pattern = files)
    file_list <- data.frame(path = file_path, file = filenames, stringsAsFactors = FALSE)

    } else {
    tool_dirs <- c("deseq2", "aldex2", "ancom2", "corncob", "maaslin2")
    file_list <- do.call(rbind, lapply(tool_dirs, function(tool) {
      tool_path <- sprintf("%s_%s/table/%s", project, format(Sys.Date(), "%y%m%d"), tool)
      if (!dir.exists(tool_path)) return(NULL)
      tool_files <- list.files(tool_path, pattern = "\\.csv$", full.names = FALSE)
      if (length(tool_files) == 0) return(NULL)
      data.frame(path = tool_path, file = tool_files, stringsAsFactors = FALSE)
    }))
  }
  if (is.null(file_list) || nrow(file_list) == 0) {
    message("No DA result files found.")
    return(invisible(NULL))
  }
  print(file_list)

  # "name" definition
  if (is.function(name)){
    name <- NULL
  }

  tt <- try(mycols,T)
  if(inherits(tt, "try-error")){
    print("mycols is not defined.")
    mycols <- NULL
  }

  plotlist_pw <- list()
  for (fn in seq_len(nrow(file_list))) {
    filename1 <- file.path(file_list$path[fn], file_list$file[fn])
    df <- read.csv(filename1, row.names=NULL ,check.names=FALSE)
    colnames(df) <- gsub('^"|"$', "", colnames(df))

    df$aldex2.FDR

    # tool recognizing
    tools <- sapply(filename1, function(filename1) {
      if (grepl("condadist", filename1)) return("condadist")
      if (grepl("wilcoxon", filename1)) return("wilcoxon")
      if (grepl("deseq2", filename1)) return("deseq2")
      if (grepl("aldex2", filename1)) return("aldex2")
      if (grepl("maaslin2", filename1)) return("maaslin2")
      if (grepl("maaslin", filename1)) return("maaslin2")
      if (grepl("corncob", filename1)) return("corncob")
      if (grepl("ancombc2", filename1)) return("ancom2")
      if (grepl("ancom2", filename1)) return("ancom2")
      return(NA) # if none of the tools are matched
    })

    unique_tools <- unique(tools[!is.na(tools)])
    print(unique_tools)

    if (unique_tools == "deseq2") {
      out_DA <- file.path(sprintf("%s_%s/pdf/DA_plot",project, format(Sys.Date(), "%y%m%d")))
      if(!dir.exists(out_DA)) dir.create(out_DA, recursive = TRUE, showWarnings = FALSE)
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
      out_DA <- file.path(sprintf("%s_%s/pdf/DA_plot",project, format(Sys.Date(), "%y%m%d")))
      if(!dir.exists(out_DA)) dir.create(out_DA, recursive = TRUE, showWarnings = FALSE)
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
      # fallback: detect model from column names when filename has no model tag
      if (length(unique_model) == 0) {
        if ("diff.btw" %in% colnames(df) || "wi.ep" %in% colnames(df)) {
          unique_model <- "t-test"
        } else if (any(grepl("\\.Est$", colnames(df)))) {
          unique_model <- "GLM"
        } else {
          unique_model <- "t-test"  # default fallback
        }
      }

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

        # holm_col_name <- sprintf("%s%s.pval.holm",unique(df$mvar), unique(df$smvar))
        holm_col_name <- sprintf("%s%s.pval.padj",unique(df$mvar), unique(df$smvar))
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
      out_DA <- file.path(sprintf("%s_%s/pdf/DA_plot",project, format(Sys.Date(), "%y%m%d")))
      if(!dir.exists(out_DA)) dir.create(out_DA, recursive = TRUE, showWarnings = FALSE)
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
    } else if (unique_tools == "corncob") {
      out_DA <- file.path(sprintf("%s_%s/pdf/DA_plot",project, format(Sys.Date(), "%y%m%d")))
      if(!dir.exists(out_DA)) dir.create(out_DA, recursive = TRUE, showWarnings = FALSE)
      x_var <- "corncob_coef"
      y_var <- "-log10(corncob_pvalue)"
      pval <- "corncob.P"
      p <- "corncob_pvalue"
      fdr <- "corncob.FDR"
      padj <- "corncob_qvalue"
      tool <- "corncob"
      model <- NULL
      print(tool)
      vx <- "coefficient"
      vy <- "-log10 (p-value)"
    } else if (unique_tools == "maaslin2") {
      out_DA <- file.path(sprintf("%s_%s/pdf/DA_plot",project, format(Sys.Date(), "%y%m%d")))
      if(!dir.exists(out_DA)) dir.create(out_DA, recursive = TRUE, showWarnings = FALSE)
      x_var <- "maaslin2_coef"
      y_var <- "-log10(maaslin2_pvalue)"
      pval <- "maaslin2.P"
      p <- "maaslin2_pvalue"
      fdr <- "maaslin2.FDR"
      padj <- "maaslin2_qvalue"
      tool <- "maaslin2"
      model <- NULL
      print(tool)
      vx <- "coefficient"
      vy <- "-log10 (p-value)"
    } else if (unique_tools == "wilcoxon") {
      out_DA <- file.path(sprintf("%s_%s/pdf/DA_plot",project, format(Sys.Date(), "%y%m%d")))
      if(!dir.exists(out_DA)) dir.create(out_DA, recursive = TRUE, showWarnings = FALSE)
      x_var <- "log2FoldChange"
      y_var <- "-log10(pvalue)"
      z_var <- "log2FoldChange"
      pval <- "wilcoxon.P"
      p <- "pvalue"
      fdr <- "wilcoxon.FDR"
      padj <- "padj"
      tool <- "wilcoxon"
      model <- NULL
      print(tool)
      vx <- "log2 fold change"
      vy <- "-log10 (p-value)"
    } else if (unique_tools == "condadist") {
      out_DA <- file.path(sprintf("%s_%s/pdf/ConDa_plot",project, format(Sys.Date(), "%y%m%d")))
      if(!dir.exists(out_DA)) dir.create(out_DA, recursive = TRUE, showWarnings = FALSE)
      x_var <- "median_effect_size"
      y_var <- "-log10(fisher_combined_p)"
      pval <- "condadist.P"
      p <- "fisher_combined_p"
      fdr <- "condadist.FDR"
      padj <- "fisher_combined_q"
      tool <- "ConDA-dist"
      model <- NULL
      print(tool)
      vx <- "Consensus median effect size"
      vy <- "-log10 (combined p-value)"
    }

    # Check if any filename has "confounder"
    has_confounder <- any(grepl("with_confounder", filename1));has_confounder

    if (has_confounder) {
      confounder <- "confounder"
      print(confounder)
    } else {
      confounder <- NULL
    }

    small_n_warn <- {
      bc <- if ("bas.count" %in% colnames(df)) suppressWarnings(as.integer(unique(df$bas.count)[1])) else NA_integer_
      sc <- if ("smvar.count" %in% colnames(df)) suppressWarnings(as.integer(unique(df$smvar.count)[1])) else NA_integer_
      if (!is.na(bc) && !is.na(sc) && (bc <= 5L || sc <= 5L)) {
        "Warning: small sample size (n\u22645) - interpret with caution"
      } else NULL
    }
    wilcoxon_fallback_note <- if (identical(unique_tools, "wilcoxon") && "notes" %in% colnames(df)) {
      note_vals <- unique(df$notes[!is.na(df$notes) & nzchar(df$notes)])
      if (length(note_vals) > 0) {
        # extract original method from note e.g. "wilcoxon fallback (ancombc2 skipped: ...)"
        m <- regmatches(note_vals[1], regexpr("wilcoxon fallback \\(([^\\s]+) skipped", note_vals[1]))
        orig <- if (length(m) > 0) sub("wilcoxon fallback \\(", "", sub(" skipped.*", "", m)) else NULL
        if (!is.null(orig) && nzchar(orig)) {
          paste0("Wilcoxon fallback (", orig, " failed)")
        } else {
          "Wilcoxon fallback"
        }
      } else NULL
    } else NULL
    subtitle_parts <- c(
      wilcoxon_fallback_note,
      if (!is.null(confounder)) "confounder-adjusted DA" else NULL,
      build_optional_metric_note(df),
      small_n_warn
    )
    subtitle_text <- if (length(subtitle_parts) > 0) paste(subtitle_parts, collapse = " | ") else NULL

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

    if (type == "taxonomy") {
      df.na$plot_label <- build_taxonomy_plot_label(df.na)
    }
    table_name_token <- NULL
    if ("name_token" %in% colnames(df.na)) {
      table_name_vals <- unique(as.character(df.na$name_token))
      table_name_vals <- table_name_vals[!is.na(table_name_vals) & table_name_vals != "" & table_name_vals != "NA"]
      if (length(table_name_vals) > 0) {
        table_name_token <- gsub("^\\(|\\)$", "", table_name_vals[1])
      }
    }

    call_name_token <- if (is.null(name)) NULL else gsub("^\\(|\\)$", "", as.character(name))

    comparison_token <- NULL
    if ("comparison_token" %in% colnames(df.na)) {
      comparison_vals <- unique(as.character(df.na$comparison_token))
      comparison_vals <- comparison_vals[!is.na(comparison_vals) & comparison_vals != "" & comparison_vals != "NA"]
      if (length(comparison_vals) > 0) {
        comparison_token <- gsub("^\\(|\\)$", "", comparison_vals[1])
      }
    }

    merged_name_token <- if (!is.null(call_name_token)) call_name_token else table_name_token



    if (model == "t-test" || model == "GLM" || is.null(model)) {
      # Get unique values for some columns
      basline_vals <- unique(as.character(df.na$basline))
      basline_vals <- basline_vals[!is.na(basline_vals) & basline_vals != "" & basline_vals != "NA"]
      smvar_vals <- unique(as.character(df.na$smvar))
      smvar_vals <- smvar_vals[!is.na(smvar_vals) & smvar_vals != "" & smvar_vals != "NA"]
      mvar_vals <- unique(as.character(df.na$mvar))
      mvar_vals <- mvar_vals[!is.na(mvar_vals) & mvar_vals != "" & mvar_vals != "NA"]

      basline <- if (length(basline_vals) > 0) basline_vals[1] else "group1"
      smvar <- if (length(smvar_vals) > 0) smvar_vals[1] else "group2"
      mvar <- if (length(mvar_vals) > 0) mvar_vals[1] else "Group"

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
      if (is.null(comparison_token)) {
        comparison_token <- sprintf("%s.vs.%s%s",
                                    basline,
                                    smvar,
                                    if (is.null(merged_name_token)) "" else paste(".", merged_name_token, sep = ""))
      } else if (!is.null(merged_name_token) &&
                 !(comparison_token == merged_name_token || endsWith(comparison_token, paste0(".", merged_name_token)))) {
        comparison_token <- paste(comparison_token, merged_name_token, sep = ".")
      }
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
      if (is.null(comparison_token)) {
        comparison_token <- sprintf("%s.vs.%s%s",
                                    basline,
                                    smvar,
                                    if (is.null(merged_name_token)) "" else paste(".", merged_name_token, sep = ""))
      } else if (!is.null(merged_name_token) &&
                 !(comparison_token == merged_name_token || endsWith(comparison_token, paste0(".", merged_name_token)))) {
        comparison_token <- paste(comparison_token, merged_name_token, sep = ".")
      }

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
      geom_vline(xintercept = c(-fc, 0, fc), col = dircolors, linetype = "dotted", size = 1)

    p1 <- p1 + scale_color_manual(values=dircolors,  labels=legend.labs, drop = FALSE)

    if (type == "taxonomy" | type == "taxanomy" | type == "bacmet") {
      label_name <- "plot_label"
      label_condition <- sprintf(
        "ifelse(df.na[, '%s'] != 'NA' & df.na[, '%s'] < 0.05 & abs(df.na[, '%s']) > fc, as.character(df.na[, '%s']), '')",
        label_name, p, x_var, label_name
      )

    } else if (type == "function") {
      label_name <- "KOName"
      label_condition <- sprintf(
        "ifelse(%s != 'NA' & df.na[, '%s'] < 0.05 & abs(df.na[, '%s']) > fc, as.character(%s), '')",
        label_name, p, x_var, label_name
      )

    } else if (type == "RNAseq") {
      label_name <- "symbol"
      label_condition <- sprintf(
        "ifelse(%s != 'NA' & df.na[, '%s'] < 0.05 & abs(df.na[, '%s']) > fc, as.character(%s), '')",
        label_name, "padj", x_var, label_name
      )

    }



    # Use the constructed label_condition in your geom_text_repel function
    p1 <- p1 + geom_text_repel(aes_string(label=label_condition), size=font, fontface="italic", max.overlaps=overlaps)

    p1 <- p1 + labs(title = sprintf("%s, %s%s (p < 0.05, cutoff=%s) ", mvar,  tool, ifelse(is.null(model), "", paste("-",model, sep = "")), fc), subtitle = subtitle_text)
    p2 <- p1 + geom_point(aes(shape=dirPadj), size=font-1.5)+  scale_shape_manual(values = padj_shape, drop = FALSE) +
      labs(shape = "FDR < 0.05", color = sprintf("%s p < 0.05",tool)) +
      guides(color = guide_legend(nrow = 1, byrow = TRUE, title.position = "top"),
             shape = guide_legend(nrow = 1, byrow = TRUE, title.position = "left")) +
      theme(text = element_text(size=font+8),
            plot.title = element_text(size=font+8),
            plot.subtitle = element_text(size=font+6, lineheight = 0.9),
            legend.text=element_text(size=font+8),
            legend.position="bottom",
            legend.justification = c(0, 0),
            legend.box.just = "left",
            legend.box = "vertical",
            legend.margin = ggplot2::margin(0, 0, 0, 0),
            legend.box.margin = ggplot2::margin(0, 0, 0, -25),
            legend.key = element_blank(),
            panel.grid = element_blank(),
            panel.background = element_rect(fill = "white", colour = "Black",size = 0.5, linetype = "solid"),
            aspect.ratio = 1/1.5)


    pdf(sprintf("%s/%s.%s%s.(%s).%s.%s%s(cutoff=%s).%s.pdf", out_DA,
                tool,
                ifelse(is.null(plot), "", paste(plot, ".", sep = "")),
                mvar,
                comparison_token,
                project,
                ifelse(is.null(model), "", paste(model, ".", sep = "")),
                ifelse(is.null(confounder), "", paste(confounder, ".", sep = "")),
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
        geom_vline(xintercept = c(-fc, 0, fc), col = dircolors, linetype = "dotted", size = 1)

      p1 <- p1 + scale_color_manual(values=dircolors,  labels=legend.labs, drop = FALSE)

      # Construct the label condition as a string using sprintf
      label_condition <- sprintf(
        "ifelse(%s != 'NA' & df.na[, '%s'] < 0.05 & abs(df.na[, '%s']) > fc, as.character(%s), '')",
        label_name, sp, a_var, label_name
      )

      # Use the constructed label_condition in your geom_text_repel function
      p1 <- p1 + geom_text_repel(aes_string(label=label_condition), size=font, fontface="italic", max.overlaps=overlaps)
      p1 <- p1 + labs(title = sprintf("%s, %s%s (p < 0.05, cutoff=%s) ", mvar,  tool, ifelse(is.null(model), "", paste("-",model, sep = "")), fc), subtitle = subtitle_text)
      p2 <- p1 + geom_point(aes(shape=dirPadj), size=font-1.5)+  scale_shape_manual(values = padj_shape, drop = FALSE) +
        labs(shape = "FDR < 0.05", color = sprintf("%s p < 0.05",tool)) +
        guides(color = guide_legend(nrow = 1, byrow = TRUE, title.position = "top"),
               shape = guide_legend(nrow = 1, byrow = TRUE, title.position = "left")) +
        theme(text = element_text(size=font+8),
              plot.title = element_text(size=font+8),
              plot.subtitle = element_text(size=font+6, lineheight = 0.9),
              legend.text=element_text(size=font+8),
              legend.position="bottom",
              legend.justification = c(0, 0),
              legend.box.just = "left",
              legend.box = "vertical",
              legend.margin = ggplot2::margin(0, 0, 0, 0),
              legend.box.margin = ggplot2::margin(0, 0, 0, -25),
              legend.key = element_blank(),
              panel.grid = element_blank(),
              panel.background = element_rect(fill = "white", colour = "Black",size = 0.5, linetype = "solid"),
              aspect.ratio = 1/1.5)



      plotlist_pw[[length(plotlist_pw) + 1]] <- p2
      if (!isTRUE(patchwork)) {
        pdf(sprintf("%s/%s.%s%s.(%s).%s.%s%s(cutoff=%s).%s.pdf", out_DA,
                    tool,
                    ifelse(is.null(plot), "", paste(plot, ".", sep = "")),
                    mvar,
                    comparison_token,
                    project,
                    ifelse(is.null(model1), "", paste(model1, ".", sep = "")),
                    ifelse(is.null(confounder), "", paste(confounder, ".", sep = "")),
                    fc,
                    format(Sys.Date(), "%y%m%d")), height = height, width = width)
        print(p2)
        dev.off()
      }
    }
  }
  if (isTRUE(patchwork)) return(invisible(plotlist_pw))
}
