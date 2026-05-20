

#' Generate Heatmaps from Phyloseq Data
#'
#' This function creates heatmaps from a Phyloseq object, allowing visualization of
#' complex microbiome data. It supports grouping, color customization, and hierarchical clustering.
#'
#' @param psIN Phyloseq object containing the data for analysis.
#' @param project Name of the project or analysis.
#' @param title Title of the heatmap.
#' @param group1 First group for annotation and coloring in the heatmap.
#' @param group2 Second group for annotation and coloring in the heatmap.
#' @param group3 Third group for annotation and coloring in the heatmap.
#' @param group4 Fourth group for annotation and coloring in the heatmap.
#' @param Ntax Number of top taxa to include in the heatmap (if NULL, all taxa are considered).
#' @param name Optional name for the analysis.
#' @param col_orders Custom order of columns in the heatmap.
#' @param show_rownames Logical value to show or hide row names.
#' @param show_colnames Logical value to show or hide column names.
#' @param cutree_rows Number of clusters for row dendrogram cutting.
#' @param cutree_cols Number of clusters for column dendrogram cutting.
#' @param cluster_rows Logical value to cluster or not cluster rows.
#' @param cluster_cols Logical value to cluster or not cluster columns.
#' @param showPhylum Logical value to show or hide phylum-level annotation.
#' @param width Width of the output heatmap.
#'
#' @return Saves heatmap plots as PDF files in a specified directory.
#'
#' @details
#' The function preprocesses the Phyloseq object, potentially considering only the top N taxa,
#' and then generates heatmaps with options for clustering, grouping, and annotation.
#' It can display phylum-level annotations and supports various customization options.
#'
#' @examples
#' # psIN is a Phyloseq object
#' # Example usage:
#' Go_pheatmap(psIN = psIN,
#'             project = "MyProject",
#'             title = "Sample Heatmap",
#'             group1 = "Treatment",
#'             group2 = "Condition",
#'             Ntax = 50,
#'             name = "Analysis1",
#'             col_orders = c("Sample1", "Sample2"),
#'             show_rownames = TRUE,
#'             show_colnames = FALSE,
#'             cluster_rows = TRUE,
#'             cluster_cols = TRUE,
#'             showPhylum = TRUE,
#'             width = 10)
#'
#' @param patchwork Logical. If \code{TRUE}, skip saving and return the plot object(s) for use with \code{Gg_patchwork()} or the \pkg{patchwork} package. Default \code{FALSE}.
#' @export

Go_pheatmap <- function(psIN,project, title = NULL,
                        group1=NULL, group2=NULL, group3=NULL, group4=NULL,
                        Ntax=NULL,
                        name=NULL,
                        col_orders=NULL,
                        show_rownames = T,show_colnames = F,
                        cutree_rows = NA, cutree_cols = NA,
                        cluster_rows = T, cluster_cols = T,
                        row_gap = NULL, column_gap = NULL,
                        showPhylum = T,
                        width,
                        patchwork = FALSE){
  # BiocManager::install("ComplexHeatmap")
  # install.packages("Cairo")
  if(!is.null(dev.list())) dev.off()

  # out dir
  out <- file.path(sprintf("%s_%s",project, format(Sys.Date(), "%y%m%d")))
  if(!dir.exists(out)) dir.create(out)
  out_pdf <- file.path(sprintf("%s_%s/pdf",project, format(Sys.Date(), "%y%m%d")))
  if(!dir.exists(out_pdf)) dir.create(out_pdf)
  out_tab <- file.path(sprintf("%s_%s/table",project, format(Sys.Date(), "%y%m%d")))
  if(!dir.exists(out_tab)) dir.create(out_tab)
  out_pheatmapTab <- file.path(sprintf("%s_%s/table/pheatmapTab",project, format(Sys.Date(), "%y%m%d")))
  if(!dir.exists(out_pheatmapTab)) dir.create(out_pheatmapTab)


  otu_raw <- as.matrix(otu_table(psIN))
  otu_has_negative <- any(otu_raw < 0, na.rm = TRUE)

  #----- normalization relative abundant ---#
   if(otu_has_negative){
     ps.rel <- psIN
     print("Signed table data detected. Skip log normalization.")
   }else if(max(otu_raw, na.rm = TRUE) < 1){
     ps.rel <- psIN
     print("Percentage abundant table data.")
   }else{
     print("Counts or cpm table data. Data normalized by log1p.")
     psIN.prune <- prune_samples(sample_sums(psIN) > 1, psIN);psIN.prune
     ps.rel <- transform_sample_counts(psIN.prune, function(x) log1p(x) );ps.rel
   }

print("Check the psIN")

  otu_rank <- as.matrix(otu_table(ps.rel))
  if (!phyloseq::taxa_are_rows(ps.rel)) otu_rank <- t(otu_rank)
  taxa_score <- if (otu_has_negative) {
    rowSums(abs(otu_rank), na.rm = TRUE)
  } else {
    rowSums(otu_rank, na.rm = TRUE)
  }

  if(is.null(Ntax)){
    Ntax <- nrow(otu_rank)
    print(sprintf("number of taxa is %s",Ntax))
    ps.rel.sel <- prune_taxa(names(sort(taxa_score, TRUE)[1:Ntax]), ps.rel);ps.rel.sel
  }else{
    ps.rel.sel <- prune_taxa(names(sort(taxa_score, TRUE)[1:Ntax]), ps.rel);ps.rel.sel
  }


  # for height

   if ( dim(tax_table(ps.rel.sel))[1] > 30){
     h=8
   }else if(dim(tax_table(ps.rel.sel))[1] <= 30 & dim(tax_table(ps.rel.sel))[1] > 25){
     h=7.4
   }else if(dim(tax_table(ps.rel.sel))[1] <= 25 & dim(tax_table(ps.rel.sel))[1] > 20){
     h=6.3
   }else if(dim(tax_table(ps.rel.sel))[1] <= 20 & dim(tax_table(ps.rel.sel))[1] > 15){
     h=5
   }else if(dim(tax_table(ps.rel.sel))[1] <= 15 & dim(tax_table(ps.rel.sel))[1] > 10){
     h=5
   }else if(dim(tax_table(ps.rel.sel))[1] <= 10){
     h=4.5
   }


  matrix <- as.matrix(otu_table(ps.rel.sel))
  if (!phyloseq::taxa_are_rows(ps.rel.sel)) matrix <- t(matrix)
  matrix <- data.frame(matrix);head(matrix)
  # normalization for log2
  is.na(matrix)<-sapply(matrix, is.infinite)
  matrix[is.na(matrix)]<-0

  matrix <- matrix[,colSums(is.na(matrix))<nrow(matrix)]
  colnames(matrix) <- gsub("X","",colnames(matrix))




 # get data tyep
  print("Check the data type")
  taxa_tab_df <- data.frame(tax_table(ps.rel.sel), stringsAsFactors = FALSE, check.names = FALSE)
  taxtab.col <- colnames(taxa_tab_df)
  taxonomy_levels <- c("Species", "Genus", "Family", "Order", "Class", "Phylum", "Kingdom")
  taxonomy_label_col <- taxonomy_levels[taxonomy_levels %in% taxtab.col][1]

  if (!is.na(taxonomy_label_col) && !is.null(taxonomy_label_col)){
    taxaTab <- data.frame(Rank = taxa_tab_df[, taxonomy_label_col], stringsAsFactors = FALSE)
    type <- "taxonomy"
    print(sprintf("%s (%s)", type, taxonomy_label_col))
  }else if("KO" %in% taxtab.col){
    taxaTab <- data.frame(Rank = taxa_tab_df[, "KO"], stringsAsFactors = FALSE)
    type <- "kegg"
    print(type)
  }else if("pathway" %in% taxtab.col){
    taxaTab <- data.frame(Rank = taxa_tab_df[, "pathway"], stringsAsFactors = FALSE)
    type <- "pathway"
    print(type)
  }else if("BileAcid" %in% taxtab.col){
    taxaTab <- data.frame(Rank = taxa_tab_df[, "BileAcid"], stringsAsFactors = FALSE)
    type <- "BileAcid"
    print(type)
  }else if("SCFA" %in% taxtab.col){
    taxaTab <- data.frame(Rank = taxa_tab_df[, "SCFA"], stringsAsFactors = FALSE)
    type <- "SCFA"
    print(type)
  }else if("symbol" %in% taxtab.col){
    taxaTab <- data.frame(Rank = taxa_tab_df[, "symbol"], stringsAsFactors = FALSE)
    type <- "RNAseq"
    print(type)
  }else if("Gene" %in% taxtab.col){
    taxaTab <- data.frame(Rank = taxa_tab_df[, "Gene"], stringsAsFactors = FALSE)
    type <- "Gene"
    print(type)
  } else {
    use_col <- if (length(taxtab.col) >= 1) taxtab.col[1] else NA
    taxaTab <- if (!is.na(use_col)) {
      data.frame(Rank = taxa_tab_df[, use_col], stringsAsFactors = FALSE)
    } else {
      data.frame(Rank = taxa_names(ps.rel.sel), stringsAsFactors = FALSE)
    }
    type <- "other"
    message(sprintf("[Go_pheatmap] Using '%s' as row labels.",
                    if (!is.na(use_col)) use_col else "rownames"))
  }

  taxaTab$Rank[is.na(taxaTab$Rank) | taxaTab$Rank == ""] <- taxa_names(ps.rel.sel)[is.na(taxaTab$Rank) | taxaTab$Rank == ""]
  taxaTab.print <- taxa_tab_df

  # 파일 이름만 생성
  file_name_csv <- sprintf("pheatmap.%s.tap(%s).%s%s%s.csv",
                           project,
                           Ntax,
                           ifelse(is.null(title), "", paste0(title, ".")),
                           ifelse(is.null(name), "", paste0(name, ".")),
                           format(Sys.Date(), "%y%m%d"))

  # 전체 경로는 file.path()로 안전하게 연결
  write.csv(taxaTab.print,
            file = file.path(out_pheatmapTab, file_name_csv),
            quote = FALSE,
            row.names = TRUE)


  print("test1")
  # map 정리
  mapping <- data.frame(sample_data(ps.rel.sel));dim(mapping)

  tt <- try(sel <- intersect(rownames(mapping), colnames(matrix)), T)

  if (length(tt) == 0){
    tt1 <- try(sel <- intersect(rownames(mapping), rownames(matrix)),T)
    if (length(tt1) == 0){
      rownames(mapping) <- gsub("-",".",rownames(mapping))
      sel <- intersect(rownames(mapping), colnames(matrix))
    }else{
      sel <- intersect(rownames(mapping), rownames(matrix))
    }
  }else{
    sel <- intersect(rownames(mapping), colnames(matrix)); head(sel, 3)
  }

  mapping.sel <- mapping[sel,, drop=F];dim(mapping.sel)


  print("test2")
  head(rownames(mapping))
  head(rownames(matrix))


  # phylum annotation
  # Define a function to assign colors to Phylum/Path
  phylumcolor <- c("#E7695DFF", "#6B8993FF", "#F6F0D4FF", "#95CE8AFF", "#D2D2D2FF", "#94D4D4FF", "#969696FF", "#F1F3E8FF", "#88775FFF")

  assign_colors <- function(x, colors) {
    getPalette = colorRampPalette(colors)
    col <- getPalette(length(unique(x)))
    names(col) = levels(x)
    return(col)
  }

  # Based on the type, generate annotation row and assign colors
  annotation_row <- NULL
  Path_col <- NULL
  phylum_col <- NULL

  if(type %in% c("taxonomy", "taxanomy") && "Phylum" %in% taxtab.col) {
    annotation_row <- data.frame(Phylum = as.factor(taxa_tab_df[, "Phylum"]))
    rownames(annotation_row) <- taxa_names(ps.rel.sel)
    phylum_col <- assign_colors(annotation_row$Phylum, phylumcolor)
  } else if(type %in% c("kegg", "pathway")) {
    if("KO" %in% taxtab.col){
      annotation_row <- data.frame(Path = as.factor(taxa_tab_df[, "KO"]))
    } else if("pathway" %in% taxtab.col){
      annotation_row <- data.frame(Path = as.factor(taxa_tab_df[, "pathway"]))
    } else {
      stop("Neither 'KO' nor 'pathway' columns are present in the taxonomy table.")
    }
    rownames(annotation_row) <- taxa_names(ps.rel.sel)
    Path_col <- assign_colors(annotation_row$Path, phylumcolor)
  }else if(type %in% c("RNAseq")) {
    annotation_row <- data.frame(Path = as.factor(taxa_tab_df[, "symbol"]))
    rownames(annotation_row) <- taxa_names(ps.rel.sel)
    Path_col <- assign_colors(annotation_row$Path, phylumcolor)
  }else if(type %in% c("BileAcid")) {
    annotation_row <- data.frame(Path = as.factor(taxa_tab_df[, "BileAcid"]))
    rownames(annotation_row) <- taxa_names(ps.rel.sel)
    Path_col <- assign_colors(annotation_row$Path, phylumcolor)
  }else if(type %in% c("SCFA")) {
    annotation_row <- data.frame(Path = as.factor(taxa_tab_df[, "SCFA"]))
    rownames(annotation_row) <- taxa_names(ps.rel.sel)
    Path_col <- assign_colors(annotation_row$Path, phylumcolor)
  }

  print("test3")
  if (!is.null(annotation_row)) {
    annotation_row <- annotation_row[rownames(matrix), , drop = FALSE]
  }



  # phylum colors
  # phylum_col = RColorBrewer::brewer.pal(length(levels(annotation_row$Phylum)), "Dark2")


  # add group(s) and color list
  # Function to generate colors
  generate_colors <- function(group, hardcode_colors=NULL) {
    unique_vals <- unique(mapping.sel[, group])
    if(is.null(hardcode_colors)) {
      stop(paste("Group", group, "not found in hardcoded_colors"))
    }
    color <- head(hardcode_colors, length(unique_vals))
    names(color) <- levels(as.factor(unique_vals))
    return(color)
  }

  # Function to create annotation_col data frame
  generate_annotation_col <- function(groups) {
    annotation_col <- as.data.frame(lapply(groups, function(x) as.factor(mapping.sel[,x])))
    colnames(annotation_col) <- groups
    return(annotation_col)
  }

  # Larger list of hardcoded colors
  all_hardcoded_colors <- list(
    color_set_1 = c("#EB5291FF", "#6351A0FF","#FEF79EFF", "#FBBB68FF", "#F5BACFFF", "#9DDAF5FF","#ECF1F4FF",  "#1794CEFF","#972C8DFF"),

    color_set_2 = c("#B15928", "#CAB2D6", "#FF7F00", "#FDBF6F", "#E31A1C", "#FB9A99", "#33A02C", "#B2DF8A", "#1F78B4", "#A6CEE3", "#1170aa", "#fc7d0b", "#76B7B2", "#B07AA1", "#E15759", "#59A14F", "#EDC948", "#FF9DA7", "#9C755F","#BAB0AC","#C84248"),
    color_set_3 = c("#2366C0FF", "#E9D738FF", "#B91226FF", "#A3DA4BFF", "#FF6435FF"),
    color_set_4 = c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666")
  )

  print("test4")

  # Define groups based on the inputs
  groups <- c(group1, group2, group3, group4)
  groups <- groups[!sapply(list(group1, group2, group3, group4), is.null)]

  row_ann_colors <- if (type %in% c("taxonomy", "taxanomy") && !is.null(phylum_col)) {
    list(Phylum = phylum_col)
  } else if (type %in% c("kegg", "pathway", "RNAseq", "BileAcid", "SCFA") && !is.null(Path_col)) {
    list(Path = Path_col)
  } else {
    list()
  }

  if (length(groups) == 0) {
    annotation_col <- NULL
    ann_colors <- row_ann_colors
  } else {
    hardcoded_colors <- all_hardcoded_colors[1:length(groups)]
    names(hardcoded_colors) <- groups
    group_colors <- lapply(groups, function(x) generate_colors(x, hardcoded_colors[[x]]))
    annotation_col <- generate_annotation_col(groups)
    ann_colors <- c(group_colors, row_ann_colors)
    if (length(ann_colors) > 0) {
      names(ann_colors) <- c(groups, names(row_ann_colors))
    }
  }



  #===== data process
  matrix <- as.matrix(matrix)

  colSums(matrix)

  bk <- c(0,0.5,1)
  print("p0")

  cutree_rows_use <- if (!is.na(cutree_rows) && cutree_rows < 2) NA else cutree_rows
  cutree_cols_use <- if (!is.na(cutree_cols) && cutree_cols < 2) NA else cutree_cols

  common_args <- list(
    mat = matrix,
    fontsize = 8,
    main = title,
    show_rownames = show_rownames,
    show_colnames = show_colnames,
    cluster_rows = cluster_rows,
    cluster_cols = cluster_cols,
    labels_row = taxaTab$Rank,
    cutree_rows = cutree_rows_use,
    cutree_cols = cutree_cols_use
  )
  if (!is.null(row_gap))    common_args$row_gap    <- grid::unit(row_gap, "mm")
  if (!is.null(column_gap)) common_args$column_gap <- grid::unit(column_gap, "mm")
  if (!is.null(annotation_col)) common_args$annotation_col <- annotation_col
  if (length(ann_colors) > 0) common_args$annotation_colors <- ann_colors

  if (showPhylum ==TRUE && !is.null(annotation_row)){
    print("with annotation_row")
    common_args$annotation_row <- annotation_row
    p <- do.call(ComplexHeatmap::pheatmap, common_args)

  } else{
    print("without annotation_row")

    if(!is.null(col_orders)){
      order <- col_orders[col_orders %in% colnames(matrix)] # match order to matrix column names
      matrix.orderd <- matrix[, order]
    }else{
      matrix.orderd <- matrix
    }

    common_args$mat <- matrix.orderd
    p <- do.call(ComplexHeatmap::pheatmap, common_args)



      print("p3")
  }


  # logic for out file
  file_name_pdf <- sprintf("pheatmap.%s.%s%s%s%s",
                           project,
                           ifelse(is.null(col_orders), "", "ordered."),
                           ifelse(is.null(title), "", paste0(title, ".")),
                           ifelse(is.null(name), "", paste0(name, ".")),
                           format(Sys.Date(), "%y%m%d"))

  if (isTRUE(patchwork)) return(invisible(p))
  pdf(file = file.path(out_pdf, paste0(file_name_pdf, ".pdf")),
      height = h, width = width)

  print("p4")
  print(p)

  dev.off()
}
