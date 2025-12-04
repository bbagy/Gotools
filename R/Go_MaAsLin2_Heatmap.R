#' Go_Maaslin2_heatmap
#'
#' Generate a longitudinal heatmap for MaAsLin2 regression coefficients
#' across Timepoints, highlighting statistically meaningful patterns
#' (p-value or q-value thresholds) with optional row clustering.
#'
#' This function takes a MaAsLin2 `all_results.tsv` table, extracts
#' only the Timepoint-associated coefficients, filters features based on
#' user-defined significance thresholds (`p_cut`, `q_cut`), expands them
#' across all timepoints, and visualizes the effect-size trajectories
#' using a clean ggplot2 heatmap.
#'
#' Significant tiles are annotated using a significance symbol system:
#' \itemize{
#'   \item `"**"` for p/q < 0.01
#'   \item `"*"`  for p/q < 0.05
#'   \item `"."`  for p/q < 0.08
#'   \item `""`   otherwise
#' }
#'
#' Timepoint coefficients are plotted by `coef` (log-effect size), while
#' significance is reflected only by tile annotation (text), not by color.
#'
#' When `cluster_rows = TRUE`, features are hierarchically clustered
#' based on Euclidean distance over their coefficient trajectories,
#' grouping together bacteria with similar longitudinal patterns.
#'
#' Output is automatically saved under:
#' \preformatted{
#'   <project>_<YYMMDD>/pdf/MaAsLin2_Heatmap_<name>_<project>_<pX>_<qY>.pdf
#' }
#'
#' @param df A MaAsLin2 result table (typically `all_results.tsv`) containing
#'        at least: `feature`, `metadata`, `value` (timepoint),
#'        `coef`, `pval`, `qval`.
#'
#' @param project Character. Project name used for directory creation and
#'        output file prefix (e.g., `"IBD"`, `"CM"`).
#'
#' @param p_cut Numeric or NULL. Optional p-value cutoff.  
#'        If provided, features with *any* timepoint p-value ≤ `p_cut`
#'        are retained.  
#'        If `NULL`, p-value filtering is disabled.
#'
#' @param q_cut Numeric or NULL. Optional q-value cutoff.  
#'        Logic is identical to `p_cut`.  
#'        If both `p_cut` and `q_cut` are supplied, filtering is applied
#'        to **either** (logical OR).
#'
#' @param orders Optional vector specifying the desired order of timepoints.
#'        If NULL, timepoints are sorted numerically.
#'
#' @param name Optional character used as a user-defined prefix for output
#'        PDF naming (e.g., `"AtoM"`, `"MtoA"`, `"Common"`).  
#'        If `NULL`, no name prefix is added.
#'
#' @param height Numeric. Height of the saved PDF plot (in inches).
#' @param width  Numeric. Width of the saved PDF plot.
#'
#' @param cluster_rows Logical. If TRUE (default), features are hierarchically
#'        clustered using their longitudinal coefficient profiles.
#'
#' @param title Optional title for the heatmap.  
#'        If NULL, title is auto-generated based on p/q thresholds.
#'
#' @return Invisibly returns the ggplot object used for the heatmap.
#'         A PDF file is saved as a side effect.
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Subsets rows where `metadata == "Timepoint"`.
#'   \item Applies p- or q-value significance filtering (feature-level).
#'   \item Expands filtered features across *all* timepoints
#'         (missing combinations are filled with coef=0).
#'   \item Optionally clusters features by effect-size similarity.
#'   \item Annotates tiles with significance markers.
#'   \item Generates a ggplot2 heatmap (`geom_tile + geom_text`).
#'   \item Saves the figure to a dated project directory.
#' }
#'
#' @examples
#' \dontrun{
#' df <- read.table("all_results.tsv", sep="\t", header=TRUE)
#'
#' Go_Maaslin2_heatmap(df,
#'     project = "CM",
#'     p_cut   = 0.1,
#'     q_cut   = NULL,
#'     orders  = c(14,30,44,58,74,88),
#'     name    = "AtoM",
#'     height  = 4,
#'     width   = 5
#' )
#' }
#'
#' @import dplyr tidyr ggplot2
#' @export

Go_Maaslin2_heatmap <- function(df,
                                project      = "Proj",
                                p_cut        = NULL,       # p-value cutoff, NULL이면 OFF
                                q_cut        = NULL,       # q-value cutoff, NULL이면 OFF
                                orders       = NULL,       # timepoint 순서 (예: c(14,30,44,58,74,88))
                                name         = NULL,       # 파일명 prefix
                                height       = 3,
                                width        = 4,
                                cluster_rows = TRUE,       # row를 coef 패턴으로 클러스터링
                                title        = NULL)       # 제목 수동 지정 가능
{
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  
  ## --------------------------
  ## 1) Output directory setup
  ## --------------------------
  date_tag <- format(Sys.Date(), "%y%m%d")
  
  out_root <- sprintf("%s_%s", project, date_tag)
  if (!utils::file_test("-d", out_root)) dir.create(out_root)
  
  out_pdf <- file.path(out_root, "pdf")
  if (!utils::file_test("-d", out_pdf)) dir.create(out_pdf)
  
  ## 파일 이름 태그
  name_tag <- if (is.null(name)) "" else paste0(name,"_")
  p_tag    <- if (is.null(p_cut)) "pOFF" else paste0("p", p_cut)
  q_tag    <- if (is.null(q_cut)) "qOFF" else paste0("q", q_cut)
  
  pdf_file <- file.path(
    out_pdf,
    sprintf("MaAsLin2_Heatmap_%s%s_%s_%s.pdf", name_tag, project, p_tag, q_tag)
  )
  
  ## --------------------------
  ## 2) Timepoint만 사용
  ## --------------------------
  df_time <- df %>%
    dplyr::filter(metadata == "Timepoint") %>%
    dplyr::mutate(value = as.numeric(as.character(value)))
  
  ## Timepoint 순서 설정
  if (!is.null(orders)) {
    tp_levels <- orders
  } else {
    tp_levels <- sort(unique(df_time$value))
  }
  
  df_time <- df_time %>%
    dplyr::mutate(value = factor(value, levels = tp_levels))
  
  ## --------------------------
  ## 3) p_cut / q_cut 필터링
  ## --------------------------
  df_sig <- df_time
  
  if (!is.null(p_cut)) {
    df_sig <- df_sig %>% dplyr::filter(pval <= p_cut)
  }
  if (!is.null(q_cut)) {
    df_sig <- df_sig %>% dplyr::filter(qval <= q_cut)
  }
  
  if (nrow(df_sig) == 0) {
    stop(sprintf(
      "No features remain after thresholds: p_cut=%s, q_cut=%s",
      ifelse(is.null(p_cut), "NULL", p_cut),
      ifelse(is.null(q_cut), "NULL", q_cut)
    ))
  }
  
  sig_features <- unique(df_sig$feature)
  
  ## --------------------------
  ## 4) 선택된 feature × 모든 timepoint 확장
  ## --------------------------
  df_expanded <- df_time %>%
    dplyr::filter(feature %in% sig_features) %>%
    dplyr::select(feature, value, coef, pval, qval) %>%
    tidyr::complete(feature, value, fill = list(
      coef = 0,
      pval = NA,
      qval = NA
    ))
  
  ## --------------------------
  ## 5) row clustering (coef 패턴 기준 feature 순서 정렬)
  ## --------------------------
  if (isTRUE(cluster_rows)) {
    mat <- df_expanded %>%
      dplyr::select(feature, value, coef) %>%
      tidyr::pivot_wider(names_from = value, values_from = coef) %>%
      as.data.frame()
    
    rownames(mat) <- mat$feature
    mat$feature <- NULL
    
    d  <- dist(mat)                 # distance: euclidean
    hc <- hclust(d, method = "complete")
    row_order <- hc$labels[hc$order]
    
    df_expanded$feature <- factor(df_expanded$feature,
                                  levels = row_order)
  } else {
    df_expanded$feature <- factor(df_expanded$feature)
  }
  
  ## --------------------------
  ## 6) Significance 코드 생성
  ##    (p_cut이 있으면 pval, 없고 q_cut만 있으면 qval)
  ## --------------------------
  df_expanded <- df_expanded %>%
    dplyr::mutate(
      p_active = dplyr::case_when(
        !is.null(p_cut) ~ pval,
        !is.null(q_cut) ~ qval,
        TRUE ~ NA_real_
      ),
      Significance = dplyr::case_when(
        is.na(p_active) ~ "",
        TRUE ~ as.character(cut(
          p_active,
          breaks = c(-Inf, 0.01, 0.05, 0.08, Inf),
          labels = c("**", "*", ".", "")
        ))
      )
    )
  
  ## --------------------------
  ## 7) Title 자동 생성
  ## --------------------------
  if (is.null(title)) {
    ptxt <- if (is.null(p_cut)) "p:OFF" else sprintf("p<%.3f", p_cut)
    qtxt <- if (is.null(q_cut)) "q:OFF" else sprintf("q<%.3f", q_cut)
    title <- sprintf("MaAsLin2 Heatmap (%s, %s)", ptxt, qtxt)
  }
  
  ## --------------------------
  ## 8) Heatmap 그리기
  ## --------------------------
  p <- ggplot(df_expanded,
              aes(x = value, y = feature, fill = coef)) +
    geom_tile(color = "grey80") +
    geom_text(aes(label = Significance),
              color = "black",
              size = 5,
              fontface = "bold") +
    scale_fill_gradient2(
      low = "blue", mid = "white", high = "red",
      midpoint = 0, name = "coef"
    ) +
    labs(
      x = "Timepoint",
      y = "Feature",
      title = title
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid = element_blank()
    )
  
  ## --------------------------
  ## 9) PDF 저장 (ggplot2::ggsave로 강제)
  ## --------------------------
  ggplot2::ggsave(filename = pdf_file, plot = p,
                  height = height, width = width)
  
  message("[Saved] ", pdf_file)
  invisible(p)
}
