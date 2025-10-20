#' Go_patternPlot
#'
#' Draw a “patient pattern” plot across ordered timepoints (or visits), showing
#' each unique pattern of observations and how many patients exhibit it.
#'
#' @description
#' Given a long data frame with an ID, an x-axis variable (e.g., timepoint), and
#' a fill/legend variable (e.g., status), the function:
#' \itemize{
#'   \item widens to one row per \code{StudyID} and one column per \code{yinfor},
#'   \item collapses to unique patterns across \code{yinfor.order},
#'   \item counts how many IDs share each pattern,
#'   \item plots patterns (left) and counts (right) in a combined figure,
#'   \item saves a PDF under \verb{<project>_<YYMMDD>/pdf/}.
#' }
#'
#' @param df A data.frame in long format containing at least the columns named
#'   by \code{StudyID}, \code{yinfor}, and \code{fillinfor}.
#' @param project Character scalar used for output folder and file naming. Default \code{"project"}.
#' @param yinfor Character scalar. Column name (string) for the x-axis grouping (e.g., timepoint).
#' @param fillinfor Character scalar. Column name (string) used for point fill/legend.
#' @param StudyID Character scalar. Column name (string) for subject/sample ID.
#' @param yinfor.order Character vector giving the display/order of \code{yinfor} levels (left-to-right on the plot).
#' @param mycols Optional character vector of colors for levels of \code{fillinfor}.
#'   If \code{NULL}, a default palette is used. Length must be \eqn{\ge} number of
#'   unique non-NA levels in \code{fillinfor}.
#' @param width Numeric. PDF width in inches. Default \code{6}.
#' @param height Numeric. PDF height in inches. Default \code{15}.
#'
#' @return (Invisibly) a list with:
#' \itemize{
#'   \item \code{plot}: the combined \code{patchwork} plot,
#'   \item \code{pattern_table}: a data.frame with unique patterns and their counts,
#'   \item \code{long_data}: the long-format pattern table used for plotting,
#'   \item \code{file}: the output PDF filepath.
#' }
#'
#' @details
#' The function creates dated output directories:
#' \verb{<project>_<YYMMDD>/pdf} and \verb{<project>_<YYMMDD>/table/dist}.
#' Missing values in \code{fillinfor} are labeled as \code{"NA"} and shown with a distinct
#' legend entry “No Sample” (shape 4). The right-hand bar chart shows how many
#' IDs share each pattern.
#'
#' @section Expected data format:
#' \itemize{
#'   \item \code{df[[StudyID]]}: identifier (character or factor).
#'   \item \code{df[[yinfor]]}: categorical (must match \code{yinfor.order}).
#'   \item \code{df[[fillinfor]]}: categorical (status/label per visit).
#' }
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' df <- data.frame(
#'   ID = rep(paste0("S", 1:8), each = 3),
#'   Visit = rep(c("V1","V2","V3"), times = 8),
#'   Status = sample(c("A","B", NA), size = 24, replace = TRUE),
#'   stringsAsFactors = FALSE
#' )
#'
#' Go_patternPlot(
#'   df,
#'   project = "Demo",
#'   yinfor = "Visit",
#'   fillinfor = "Status",
#'   StudyID = "ID",
#'   yinfor.order = c("V1","V2","V3"),
#'   mycols = c(A = "#1f77b4", B = "#ff7f0e"),
#'   width = 7, height = 10
#' )
#' }
#'
#' @importFrom dplyr select group_by summarise across mutate %>%
#' @importFrom tidyr pivot_wider pivot_longer replace_na
#' @importFrom ggplot2 ggplot aes geom_line geom_point geom_col geom_text
#' @importFrom ggplot2 scale_shape_manual scale_fill_manual scale_x_discrete
#' @importFrom ggplot2 theme_classic theme element_text labs coord_flip
#' @importFrom ggplot2 ggsave
#' @importFrom patchwork plot_layout plot_annotation
#' @importFrom rlang sym
#' @importFrom scales hue_pal
#' @importFrom grDevices cairo_pdf
#' @export


Go_patternPlot <- function(
    df,
    project = "project",
    yinfor,                 # x축(가로) 구분 변수명 (문자열)
    fillinfor,              # 점의 채움/범례 변수명 (문자열)
    StudyID,                # 환자/샘플 ID 변수명 (문자열)
    yinfor.order,           # x축 순서 (문자 벡터)
    mycols = NULL,          # 색상 팔레트 (선택)
    width = 6,              # 저장 폭(inch)
    height = 15             # 저장 높이(inch)
){
  require(dplyr)
  require(tidyr)
  require(ggplot2)
  require(patchwork)
  require(rlang)

  # 날짜 폴더 자동 생성
  stamp <- format(Sys.Date(), "%y%m%d")
  base <- sprintf("%s_%s", project, stamp)
  pdf_dir  <- file.path(base, "pdf")
  pattern_dir <- file.path(base, "table", "pattern")
  dir.create(pdf_dir,  recursive = TRUE, showWarnings = FALSE)
  dir.create(pattern_dir, recursive = TRUE, showWarnings = FALSE)

  # 채움 값
  fill_vals <- sort(unique(na.omit(df[[fillinfor]])))
  if (is.null(mycols)) {
    mycols <- scales::hue_pal()(length(fill_vals))
  } else if (length(mycols) < length(fill_vals)) {
    stop("mycols 길이가 fill 범주의 개수보다 적습니다.")
  }
  fill_map   <- setNames(mycols[seq_along(fill_vals)], fill_vals)
  shape_map  <- setNames(rep(21, length(fill_vals)), fill_vals)
  legend_lab <- setNames(fill_vals, fill_vals)

  # wide 변환
  data_wide <- df %>%
    select(all_of(c(StudyID, yinfor, fillinfor))) %>%
    pivot_wider(id_cols = all_of(StudyID),
                names_from = all_of(yinfor),
                values_from = all_of(fillinfor))

  # 패턴 카운트
  unique_patterns <- data_wide %>%
    group_by(across(all_of(yinfor.order))) %>%
    summarise(pattern_count = n(), .groups = "drop") %>%
    mutate(across(all_of(yinfor.order), ~ replace_na(as.character(.), "NA"))) %>%
    mutate(ID = row_number())

  # long 변환
  data_long <- unique_patterns %>%
    pivot_longer(cols = all_of(yinfor.order),
                 names_to = yinfor,
                 values_to = fillinfor)

  # factor 처리
  data_long[[yinfor]] <- factor(data_long[[yinfor]], levels = yinfor.order)
  desired_order <- unique(data_long$ID)
  data_long$ID <- factor(data_long$ID, levels = rev(desired_order))
  unique_patterns$ID <- factor(unique_patterns$ID, levels = rev(desired_order))

  # NA 항목 범례 추가
  fill_all_levels  <- c(fill_vals, "NA")
  shape_all        <- c(rep(21, length(fill_vals)), 4)
  color_all        <- c(unname(fill_map), "black")
  legend_all       <- c(fill_vals, "No Sample")
  names(shape_all) <- names(color_all) <- names(legend_all) <- fill_all_levels

  # p1
  p1 <- ggplot(data_long, aes(x = !!sym(yinfor), y = ID, group = ID)) +
    geom_line(color = "grey70", linewidth = 0.8) +
    geom_point(aes(shape = !!sym(fillinfor), fill = !!sym(fillinfor)),
               color = "black", size = 3, stroke = 0.6) +
    scale_shape_manual(values = shape_all, labels = legend_all) +
    scale_fill_manual(values = color_all, labels = legend_all, drop = FALSE) +
    scale_x_discrete(position = "top") +
    labs(y = "Patterns", x = NULL) +
    theme_classic(base_size = 11) +
    theme(axis.text.x  = element_text(hjust = 0.5),
          legend.title = element_blank())

  # p3
  unique_patterns$pattern_count <- as.numeric(unique_patterns$pattern_count)
  p3 <- ggplot(unique_patterns, aes(x = ID, y = pattern_count)) +
    geom_col(fill = "#149BEDFF", width = 0.85) +
    geom_text(aes(label = pattern_count), hjust = 1.1, vjust = 0.5, size = 3) +
    coord_flip() +
    labs(y = "Count of Patients", x = NULL) +
    theme_classic(base_size = 11) +
    theme(axis.text.y  = element_blank(),
          axis.title.y = element_blank())

  # 합성
  title_txt <- sprintf("%s (n = %s)", fillinfor, sum(unique_patterns$pattern_count, na.rm = TRUE))
  p_final <- (p1 + p3) +
    plot_layout(widths = c(3, 1), guides = "collect") +
    plot_annotation(title = title_txt,
                    theme = theme(plot.title = element_text(hjust = 0.5, face = "bold")))

  # 저장 (PDF)
  outfile <- file.path(pdf_dir, sprintf("pattern_plot.%s.%s.%s.pdf", fillinfor, project, stamp))
  ggsave(outfile, plot = p_final, device = cairo_pdf, width = width, height = height, units = "in", limitsize = FALSE)

  write.csv(pattern, sprintf("%s/pattern_plot.%s.%s.%s.csv", pattern_dir, fillinfor, project, stamp))

  invisible(list(plot = p_final, pattern_table = unique_patterns, long_data = data_long, file = outfile))
}
