#' Go_Groupheatmap
#'
#' This function generates a heatmap for the top `n` features of a data matrix, grouped by a specified category. It allows for normalization of the data and customization of plot appearance.
#'
#' @param df A data frame containing the data to be visualized. Rows represent features, and columns represent samples.
#' @param SampleData A data frame with sample metadata. Row names should match the column names of `df`.
#' @param project A character string specifying the project name, used for file naming.
#' @param Group A character string specifying the column name in `SampleData` that contains group information for the samples.
#' @param orders A character vector specifying the order of groups to display on the heatmap. If NULL, groups are displayed in descending order of their sizes.
#' @param title A character string specifying the title of the heatmap.
#' @param name An optional character string for additional naming in the output file.
#' @param top_n An integer specifying the number of top features to include in the heatmap based on mean depth. Default is 30. If NULL, all features are included.
#' @param normalization A character string specifying the normalization method. Options are "log" (log transformation and scaling), "Zscore" (z-score normalization), or "none" (no normalization). Default is "log".
#' @param width An integer specifying the width of the output PDF file in inches. Default is 16.
#' @param height An integer specifying the height of the output PDF file in inches. Default is 8.
#' @param x_label An optional character string specifying the column name in `SampleData` to be used as x-axis labels in the heatmap. If NULL, the group information is used as labels.
#'
#' @return The function saves a PDF file with the heatmap and returns nothing.
#'
#' @details
#' The function performs the following steps:
#' \itemize{
#'   \item Selects the top `n` features based on mean depth or uses all features if `top_n` is NULL.
#'   \item Reorders the data based on common sample IDs.
#'   \item Applies the specified normalization method.
#'   \item Generates a heatmap and saves it as a PDF file with the specified dimensions.
#' }
#'
#' @importFrom dplyr arrange as.data.frame filter mutate select slice table
#' @importFrom RColorBrewer brewer.pal
#' @importFrom ggplot2 ggplot aes geom_tile scale_fill_gradientn theme element_blank
#' @importFrom magrittr %>%
#' @importFrom utils write.csv
#'
#' @examples
#' \dontrun{
#' df <- your_data_matrix
#' SampleData <- your_sample_metadata
#' Go_Groupheatmap(df, SampleData, project = "Project1", Group = "Group", orders = c("Group1", "Group2"), title = "Heatmap Example", top_n = 20, normalization = "log", width = 12, height = 10)
#' }
#' @export
#'
Go_Groupheatmap <- function(df, SampleData, project, Group,
                            orders,
                            title=NULL,
                            mycols=NULL,
                            name = NULL,
                            top_n = 30, normalization = "log", width = 16, height = 8,
                            x_label = NULL) {

  if (!is.null(dev.list())) dev.off()

  # out dir
  out <- file.path(sprintf("%s_%s",project, format(Sys.Date(), "%y%m%d")))
  if(!file_test("-d", out)) dir.create(out)
  out_path <- file.path(sprintf("%s_%s/pdf",project, format(Sys.Date(), "%y%m%d")))
  if(!file_test("-d", out_path)) dir.create(out_path)

  tt <- try(mycols,T)
  if(class(tt) == "try-error"){
    print("mycols is not defined.")
    mycols <- NULL
  }

  # 상위 n개의 항목 선택
  if (!is.null(top_n)) {
    df2_top <- df %>%
      dplyr::mutate(mean_depth = base::rowMeans(dplyr::select(., dplyr::everything()), na.rm = TRUE)) %>%
      dplyr::arrange(desc(mean_depth)) %>%
      dplyr::slice(1:top_n) %>%
      dplyr::select(-mean_depth)
  } else {
    df2_top <- df  # top_n이 NULL이면 전체 데이터를 사용
  }

  # 공통된 SampleID 찾기
  common_samples <- intersect(rownames(SampleData), colnames(df2_top))

  # df2의 열을 공통된 SampleID 순서에 따라 재정렬
  df2_ordered <- df2_top[, match(common_samples, colnames(df2_top))]

  # SampleData도 공통된 SampleID로 필터링
  ordered_SampleData <- SampleData[common_samples, ]

  # 그룹 정보 추출
  group_info <- ordered_SampleData[, Group]

  # 그룹별 샘플 개수를 계산하고 정렬
  group_sizes <- table(group_info)
  sorted_groups <- names(sort(group_sizes, decreasing = TRUE))

  # group_info를 orders 순서로 재정렬
  if (!is.null(orders)) {
    sorted_groups <- orders
  }

  sorted_groups <- intersect(orders, unique(group_info))
  group_info <- factor(group_info, levels = sorted_groups)



  # group_info에 맞는 SampleID 순서로 df2 재정렬
  df2_ordered <- df2_ordered[, match(rownames(ordered_SampleData), colnames(df2_ordered))]

  # 정규화 적용
  if (normalization == "log") {
    df2_ordered <- log1p(df2_ordered)
    df2_ordered <- scale(df2_ordered)
  } else if (normalization == "Zscore") {
    df2_ordered <- t(scale(t(df2_ordered)))
  }

  print("There are three options for normalization: 'log', 'Zscore', or 'none'.")

  # 열 라벨 설정: x_label이 제공되면 그 열(column)을 사용
  if (!is.null(x_label)) {
    sample_labels <- as.character(ordered_SampleData[, x_label])
  } else {
    sample_labels <- NULL  # 기본값은 group_info를 사용
  }

  # 색상 팔레트 설정 (Set3)


  if(!is.null(mycols)){
    palette <- mycols
  }else{
    palette <- brewer.pal(length(sorted_groups), "Set3")
  }



  # PDF 저장 경로 설정

  pdf(sprintf("%s/heatmap.%s.%s%s%s.%s.pdf", out_path, project,
              ifelse(is.null(name), "", paste(name, ".", sep = "")),
              ifelse(is.null(title), "", paste(title, ".", sep = "")),
              normalization, format(Sys.Date(), "%y%m%d")), width = width, height = height)

  # heatmap 생성
  if(!is.null(x_label)){
    ht <- Heatmap(as.matrix(df2_ordered),
                  name = "expression",
                  column_title = sprintf("Profile of Antibiotic Resistance Genes in %s Bacteria", title),
                  cluster_rows = TRUE,
                  cluster_columns = FALSE,
                  show_column_names = TRUE,
                  column_labels = sample_labels,  # 샘플 라벨을 변경
                  show_row_names = TRUE,
                  column_split = group_info,
                  top_annotation = HeatmapAnnotation(
                    foo = anno_block(
                      gp = gpar(fill = palette),
                      labels = as.character(sorted_groups)
                    )
                  )
    )
  }else{
    ht <- Heatmap(as.matrix(df2_ordered),
                  name = "expression",
                  column_title = sprintf("Profile of Antibiotic Resistance Genes in %s Bacteria", title),
                  cluster_rows = TRUE,
                  cluster_columns = FALSE,
                  show_column_names = TRUE,
                  show_row_names = TRUE,
                  column_split = group_info,
                  top_annotation = HeatmapAnnotation(
                    foo = anno_block(
                      gp = gpar(fill = palette),
                      labels = as.character(sorted_groups)
                    )
                  )
    )

  }

  draw(ht)
  dev.off()
}
