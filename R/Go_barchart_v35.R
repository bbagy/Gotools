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
#' @param mark mark "*" for the X-axis.
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
#'             mark = NULL,
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
                        mark=NULL,
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

  estimate_barchart_layout <- function(df_plot, tax_var, mvar, facet_vars = NULL,
                                       x_var = "SampleIDfactor", simple = FALSE,
                                       user_ncol = NULL, panel_height_base = 4.0) {
    if (!is.null(facet_vars) && mvar %in% facet_vars) {
      return(NULL)
    }

    facet_vars <- if (is.null(facet_vars)) character(0) else setdiff(as.character(facet_vars), "SampleType")
    panel_group_vars <- unique(c(facet_vars, if (isTRUE(simple)) character(0) else mvar))

    if (length(panel_group_vars) == 0) {
      panel_count <- 1
      panel_ids <- factor(rep("all", nrow(df_plot)))
    } else {
      panel_ids <- interaction(df_plot[, panel_group_vars, drop = FALSE], drop = TRUE, lex.order = TRUE)
      panel_count <- nlevels(panel_ids)
      if (panel_count == 0) panel_count <- 1
    }

    facet_ncol <- if (!is.null(user_ncol)) {
      max(1, user_ncol)
    } else {
      NA_integer_
    }
    facet_nrow <- if (is.na(facet_ncol)) 1 else max(1, ceiling(panel_count / facet_ncol))

    if (x_var %in% c("SampleID", "SampleIDfactor")) {
      samples_per_panel <- max(as.numeric(table(panel_ids)))
    } else {
      samples_per_panel <- max(tapply(as.character(df_plot[[x_var]]), panel_ids, function(v) length(unique(v))))
      if (!is.finite(samples_per_panel)) {
        samples_per_panel <- length(unique(df_plot[[x_var]]))
      }
    }
    if (!is.finite(samples_per_panel) || samples_per_panel < 1) {
      samples_per_panel <- 1
    }

    labels <- as.character(unique(df_plot[[tax_var]]))
    labels <- labels[!is.na(labels) & nzchar(labels)]
    colourCount <- length(labels)
    label_nchar <- if (colourCount > 0) nchar(labels) else 0
    max_label_nchar <- if (colourCount > 0) max(label_nchar, na.rm = TRUE) else 0

    legend_ncol_cap <- if (max_label_nchar <= 12) {
      6
    } else if (max_label_nchar <= 20) {
      5
    } else if (max_label_nchar <= 30) {
      4
    } else if (max_label_nchar <= 45) {
      3
    } else {
      2
    }
    legend_ncol <- max(1, min(colourCount, legend_ncol_cap))
    legend_rows <- max(1, ceiling(colourCount / legend_ncol))

    legend_text_size <- 8
    legend_key_size <- if (max_label_nchar > 40) {
      0.14
    } else if (max_label_nchar > 25) {
      0.18
    } else {
      0.20
    }

    list(
      facet_ncol = facet_ncol,
      facet_nrow = facet_nrow,
      legend_ncol = legend_ncol,
      legend_text_size = legend_text_size,
      legend_key_size = legend_key_size
    )
  }

  build_barchart_pdf_path <- function(rank_name) {
    rank_tag <- if (length(taxanames) > 1) paste(rank_name, ".", sep = "") else ""
    if (relative == TRUE) {
      sprintf("%s/barchart.relative.%s.%s%s(%s).%s%s.pdf", out_path,
              project,
              ifelse(is.null(facet), "", paste(facet, ".", sep = "")),
              ifelse(is.null(name), "", paste(name, ".", sep = "")),
              cutoff,
              rank_tag,
              format(Sys.Date(), "%y%m%d"))
    } else {
      sprintf("%s/barchart.absolute.%s.%s%s(%s).%s%s.pdf", out_path,
              project,
              ifelse(is.null(facet), "", paste(facet, ".", sep = "")),
              ifelse(is.null(name), "", paste(name, ".", sep = "")),
              cutoff,
              rank_tag,
              format(Sys.Date(), "%y%m%d"))
    }
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

  for(i in seq_along(taxanames)){

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

    agg <- aggregate(as.formula(sprintf(". ~ %s" , taxanames[i])), otu.filt, sum, na.action=na.pass)

    if (taxanames[i] == "Species"){
      agg <- agg[grepl("NA NA", agg$Species)==F,]
    }

    if (relative == TRUE){
      genera <- agg[,taxanames[i]]
      agg[,taxanames[i]] <- NULL
      agg <- normalizeByCols(agg)
      inds_to_grey <- which(rowMeans(agg) < cutoff)
      genera[inds_to_grey] <- "[1_#Other]"
      agg[,taxanames[i]] <- genera
      agg_other_out <- subset(agg, agg[,taxanames[i]] != "[1_#Other]")
      write.csv(agg_other_out, quote = FALSE,file=sprintf("%s/%s.taxa_relative_abundance.(%s).%s.%s%s.csv", out_taxa,
                                                          project,cutoff,taxanames[i],
                                                          ifelse(is.null(name), "", paste(name, ".", sep = "")),
                                                          format(Sys.Date(),"%y%m%d")))
      df <- melt(agg, variable="SampleID")
    }else{
      genera <- agg[,taxanames[i]]
      agg[,taxanames[i]] <- NULL
      agg.rel <- normalizeByCols(agg)
      inds_to_grey <- which(rowMeans(agg.rel) < cutoff)
      genera[inds_to_grey] <- "[1_#Other]"
      agg[,taxanames[i]] <- genera
      agg_other_out <- subset(agg, agg[,taxanames[i]] != "[1_#Other]")
      write.csv(agg_other_out, quote = FALSE, file=sprintf("%s/%s.taxa_absolute_abundance.(%s).%s.%s%s.csv", out_taxa,
                                                           project,cutoff,taxanames[i],
                                                           ifelse(is.null(name), "", paste(name, ".", sep = "")),
                                                           format(Sys.Date(),"%y%m%d")))
      df <- melt(agg, variable="SampleID")
    }

    df2 <- aggregate(as.formula(sprintf("value ~ %s + SampleID" , taxanames[i])), df, sum)
    df2$SampleID <- as.character(df2$SampleID)
    df2$SampleIDfactor <- factor(df2$SampleID, levels=ordering.pc1)

    for (mvar in cate.vars) {
      df2[,mvar] <- mapping.sel[df2$SampleID, mvar]
      if (length(orders) >= 1) {
        df2[,mvar] <- factor(df2[,mvar], levels = orders)
      }
      else {
        df2[,mvar] <- factor(df2[,mvar])
      }
    }

    if (!is.null(facet)) {
      for (fa in facet){
        rownames(mapping.sel) <- as.character(rownames(mapping.sel))
        df2[,fa] <- mapping.sel[df2$SampleID, fa]
        if (length(orders) >= 1) {
          df2[,fa] <- factor(df2[,fa], levels = orders)
        } else {
          df2[,fa] <- factor(df2[,fa])
        }
      }
    }

    if (x_label == "SampleID" | x_label == "SampleIDfactor") {
      df2 <- df2
    } else if (length(x_label) >= 1) {
      df2[, x_label] <- mapping.sel[df2$SampleID, x_label]

      if (!is.null(mark)) {
        df2[, mark] <- mapping.sel[df2$SampleID, mark]
        df2$x_label_with_star <- ifelse(df2[[mark]] == "Yes", paste0(df2[[x_label]], " *"), as.character(df2[[x_label]]))
        new_orders <- unlist(lapply(orders, function(order) c(order, paste0(order, " *"))))
        df2[, x_label] <- factor(df2$x_label_with_star, levels = new_orders)
      } else {
        df2[, x_label] <- factor(df2[, x_label], levels = orders)
      }
    }

    colourCount <- length(unique(df2[,taxanames[i]]));colourCount
    if(!is.null(mycols)){
      getPalette <- colorRampPalette(mycols)
    }

    layout_plan <- list()
    for (mvar in cate.vars) {
      layout_info <- estimate_barchart_layout(
        df_plot = df2,
        tax_var = taxanames[i],
        mvar = mvar,
        facet_vars = facet,
        x_var = x_label,
        simple = simple,
        user_ncol = ncol,
        panel_height_base = height
      )
      if (!is.null(layout_info)) {
        layout_plan[[mvar]] <- layout_info
      }
    }
    if (length(layout_plan) == 0) {
      next
    }

    pdf_width <- width
    pdf_height <- height
    pdf(build_barchart_pdf_path(taxanames[i]), height = pdf_height, width = pdf_width)

    for (mvar in cate.vars) {
      current_layout <- layout_plan[[mvar]]
      if (is.null(current_layout)) {
        next
      }

      p <- ggplot(df2, aes(x = !!sym(x_label), y = value, fill = !!sym(taxanames[i]), order = !!sym(taxanames[i]))) +
        geom_bar(stat = "identity", position = "stack") +
        theme_classic() +
        labs(fill = NULL) +
        theme(
          legend.position = "bottom",
          legend.text = element_text(face = "italic", size = current_layout$legend_text_size),
          legend.key.height = grid::unit(current_layout$legend_key_size, "cm"),
          legend.key.width = grid::unit(current_layout$legend_key_size, "cm"),
          legend.spacing.y = grid::unit(0.02, "cm"),
          legend.margin = ggplot2::margin(0, 0, 0, 0, "cm"),
          axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8)
        ) +
        guides(fill = ggplot2::guide_legend(ncol = current_layout$legend_ncol, byrow = TRUE))

      if(!is.null(y_axis)){
        p <- p + labs(y = y_axis)
      }else if (relative == TRUE){
        p <- p + labs(y = "Relative abundance") + ylim(c(-.1, 1.01))
      }else{
        p <- p + labs(y = "Absolute abundance")
      }

      if(!is.null(mycols)){
        p <- p + scale_fill_manual(values = getPalette(colourCount))
      }

      p <- p + ggtitle(sprintf("%s barplots overall of %s%s (cut off < %s) %s",
                               taxanames[i],
                               mvar,
                               ifelse(is.null(name), "", paste0("-", name)),
                               cutoff,
                               ifelse(is.null(mark), "", paste0("mark by - ", mark))))

      if (!is.null(facet)) {
        facet_vars <- setdiff(as.character(facet), "SampleType")
        facet_formula <- as.formula(sprintf("~ %s", paste(c(facet_vars, mvar), collapse = " + ")))
        if (is.na(current_layout$facet_ncol)) {
          p <- p + facet_grid(facet_formula, scales = "free_x", space = "free_x")
        } else {
          p <- p + facet_wrap(facet_formula, scales = "free_x", ncol = current_layout$facet_ncol)
        }
      } else if (simple == FALSE) {
        if (is.na(current_layout$facet_ncol)) {
          p <- p + facet_grid(as.formula(sprintf("~ %s", mvar)), scales = "free_x", space = "free_x")
        } else {
          p <- p + facet_wrap(as.formula(sprintf("~ %s", mvar)), scales = "free_x", ncol = current_layout$facet_ncol)
        }
      }

      print(p)
    }
    dev.off()
  }
}
