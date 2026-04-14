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
#' @param orders Optional character vector for categorical ordering
#'   (used for \code{fillinfor} and \code{multi.sites} via \code{intersect()}).
#' @param multi.sites Optional character scalar. Column name (string) for site subgroup
#'   (e.g., Nose/Mouth/Stool). If \code{NULL}, default single-point mode is used.
#'   If provided, each timepoint is shown as side-by-side squares per site:
#'   filled = sample present, empty = no sample.
#' @param pattern Logical. If \code{TRUE} (default), aggregate to unique patterns and count patients.
#'   If \code{FALSE}, plot each \code{StudyID} row directly (no pattern collapsing).
#' @param site.compact.step Numeric. Horizontal spacing step for site squares in
#'   \code{multi.sites + pattern=TRUE} mode. Larger values spread squares more.
#' @param site.dodge.width Numeric. Dodge width for site squares in
#'   \code{multi.sites + pattern=FALSE} mode. Larger values spread squares more.
#' @param mycols Optional character vector of colors for levels of \code{fillinfor}.
#'   If \code{NULL}, a default palette is used. Length must be \eqn{\ge} number of
#'   unique non-NA levels in \code{fillinfor}.
#' @param subtitle Optional character scalar shown under the main title.
#'   If \code{multi.sites} is used, this text is appended after the automatic
#'   site-order subtitle.
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
#' @importFrom ggplot2 theme_classic theme element_text labs coord_flip guide_legend
#' @importFrom ggplot2 ggsave
#' @importFrom patchwork plot_layout plot_annotation
#' @importFrom rlang sym
#' @importFrom scales hue_pal
#' @importFrom grDevices cairo_pdf
#' @param patchwork Logical. If \code{TRUE}, skip saving and return the plot object(s) for use with \code{Gg_patchwork()} or the \pkg{patchwork} package. Default \code{FALSE}.
#' @export


Go_patternPlot <- function(
    df,
    project = "project",
    yinfor,                 # x축(가로) 구분 변수명 (문자열)
    fillinfor = NULL,       # 점의 채움/범례 변수명 (문자열)
    StudyID,                # 환자/샘플 ID 변수명 (문자열)
    yinfor.order,           # x축 순서 (문자 벡터)
    orders = NULL,          # 정렬 통합 설정(list: yinfor/fillinfor/site)
    multi.sites = NULL,     # 사이트 서브그룹 변수명 (문자열, 선택)
    pattern = TRUE,         # TRUE: 패턴 집계 / FALSE: StudyID별 직접 표시
    site.compact.step = 0.14, # multi.sites + pattern=TRUE 네모 간격
    site.dodge.width = 0.6,   # multi.sites + pattern=FALSE dodge 간격
    mycols = NULL,          # 색상 팔레트 (선택)
    subtitle = NULL,        # 부제목 (선택)
    name = NULL,
    width = 6,              # 저장 폭(inch)
    height = 15,            # 저장 높이(inch)
    patchwork = FALSE
){
  tt <- try(orders, TRUE)
  if (inherits(tt, "try-error")) {
    print("orders is not defined.")
    orders <- NULL
  }

  # --- NEW 1: fillinfor = NULL이면 yinfor를 자동 fill로 사용 ---
  fillinfor_user <- fillinfor
  auto_fill_mode <- FALSE
  if (is.null(fillinfor)) {
    df$.__auto_fill__ <- df[[yinfor]]
    fillinfor <- ".__auto_fill__"
    auto_fill_mode <- TRUE
  }
  site_levels_for_subtitle <- NULL

  # 날짜 폴더 자동 생성
  name_user <- if (!is.null(name) && !is.na(name) && name != "") name else NULL
  name <- ifelse(is.null(name), "", paste(name, ".", sep = ""))
  stamp <- format(Sys.Date(), "%y%m%d")
  base <- sprintf("%s_%s", project, stamp)
  pdf_dir  <- file.path(base, "pdf")
  pattern_dir <- file.path(base, "table", "pattern")
  dir.create(pdf_dir,  recursive = TRUE, showWarnings = FALSE)
  dir.create(pattern_dir, recursive = TRUE, showWarnings = FALSE)

  if (is.null(multi.sites)) {
    # --- 기본 모드: 기존 점 패턴 플롯 ---
    df <- as.data.frame(df[, c(StudyID, yinfor, fillinfor), drop = FALSE], stringsAsFactors = FALSE)
    df <- df %>%
      distinct(across(all_of(c(StudyID, yinfor))), .keep_all = TRUE)

    # 채움 값 정렬
    fill_vals_raw <- unique(na.omit(as.character(df[[fillinfor]])))
    if (!is.null(orders)) {
      fill_vals_in_user <- intersect(orders, fill_vals_raw)
      fill_vals <- c(fill_vals_in_user, setdiff(fill_vals_raw, fill_vals_in_user))
    } else {
      fill_vals_in_yorder <- yinfor.order[yinfor.order %in% fill_vals_raw]
      fill_vals <- c(fill_vals_in_yorder, setdiff(fill_vals_raw, fill_vals_in_yorder))
    }
        if (is.null(mycols)) {
            mycols <- scales::hue_pal()(length(fill_vals))
        } else if (length(mycols) < length(fill_vals)) {
            stop("Length of 'mycols' is smaller than the number of fill categories.")
        }
        if (!is.null(names(mycols)) && all(fill_vals %in% names(mycols))) {
            fill_map <- mycols[fill_vals]
        } else {
            fill_map <- setNames(mycols[seq_along(fill_vals)], fill_vals)
        }

    # wide 변환
    data_wide <- df %>%
      select(all_of(c(StudyID, yinfor, fillinfor))) %>%
      pivot_wider(id_cols = all_of(StudyID),
                  names_from = all_of(yinfor),
                  values_from = all_of(fillinfor))

    # 패턴 카운트 또는 StudyID별 표시
    if (isTRUE(pattern)) {
      unique_patterns <- data_wide %>%
        group_by(across(all_of(yinfor.order))) %>%
        summarise(pattern_count = n(), .groups = "drop") %>%
        mutate(across(all_of(yinfor.order), ~ replace_na(as.character(.), "NA"))) %>%
        mutate(ID = row_number())
    } else {
      unique_patterns <- data_wide %>%
        mutate(across(all_of(yinfor.order), ~ replace_na(as.character(.), "NA")))
      sample_mat <- as.data.frame(unique_patterns[, yinfor.order, drop = FALSE], stringsAsFactors = FALSE)
      unique_patterns$pattern_count <- rowSums(sample_mat != "NA")
      unique_patterns$ID <- as.character(unique_patterns[[StudyID]])
    }

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
    if (!isTRUE(pattern)) {
      id_levels <- levels(data_long$ID)
      n_ids <- length(id_levels)
      ncol_split <- if (n_ids <= 1) 1L else 2L
      ids_per_col <- ceiling(n_ids / ncol_split)
      id_pos <- match(as.character(data_long$ID), id_levels)
      data_long$panel_col <- factor(ifelse(id_pos <= ids_per_col, "Col 1", "Col 2"),
                                    levels = if (ncol_split == 1L) "Col 1" else c("Col 1", "Col 2"))
    }

    # NA 항목 범례 추가
    fill_all_levels  <- c(fill_vals, "NA")
    shape_all        <- c(rep(21, length(fill_vals)), 4)
    color_all        <- c(unname(fill_map), "black")
    legend_all       <- c(fill_vals, "No Sample")
    names(shape_all) <- names(color_all) <- names(legend_all) <- fill_all_levels
    data_long[[fillinfor]] <- factor(data_long[[fillinfor]], levels = fill_all_levels)

    # p1
    p1 <- ggplot(data_long, aes(x = !!sym(yinfor), y = ID, group = ID)) +
      geom_line(color = "grey70", linewidth = 0.8) +
      geom_point(aes(shape = !!sym(fillinfor), fill = !!sym(fillinfor)),
                 color = "black", size = 3, stroke = 0.6) +
      scale_shape_manual(values = shape_all, labels = legend_all, breaks = fill_all_levels) +
      scale_fill_manual(values = color_all, labels = legend_all, breaks = fill_all_levels, drop = FALSE) +
      guides(
        fill = guide_legend(order = 1, override.aes = list(shape = unname(shape_all), size = 3)),
        shape = "none"
      ) +
      scale_x_discrete(position = "top") +
      labs(y = "Patterns", x = NULL) +
      theme_classic(base_size = 11) +
      theme(axis.text.x  = element_text(hjust = 0.5),
            legend.title = element_blank(),
            legend.position = "right")
    if (!isTRUE(pattern)) {
      p1 <- p1 +
        ggplot2::facet_wrap(ggplot2::vars(panel_col), nrow = 1, scales = "free_y") +
        theme(strip.background = element_blank(),
              strip.text.x = element_blank())
    }

  } else {
    # --- multi.sites 모드: 시점마다 사이트 네모를 나란히 표시 ---
    if (!multi.sites %in% names(df)) {
      stop("The 'multi.sites' column is not found in df: ", multi.sites)
    }

    site_vec <- df[[multi.sites]]
    if (is.factor(site_vec)) {
      site_levels <- levels(site_vec)
      site_levels <- site_levels[site_levels %in% as.character(na.omit(site_vec))]
    } else {
      site_levels <- unique(as.character(na.omit(site_vec)))
    }
    if (!is.null(orders)) {
      site_in_user <- intersect(orders, site_levels)
      site_levels <- c(site_in_user, setdiff(site_levels, site_in_user))
    }
    if (length(site_levels) == 0) {
      stop("No valid (non-NA) values were found in 'multi.sites'.")
    }
    site_levels_for_subtitle <- site_levels

    color_by_fill <- !isTRUE(auto_fill_mode)
    color_source_col <- if (color_by_fill) fillinfor else yinfor

    # 존재 여부(1/0) 테이블 생성
    df_ms <- as.data.frame(df[, c(StudyID, yinfor, multi.sites, color_source_col), drop = FALSE], stringsAsFactors = FALSE) %>%
      mutate(.ID = as.character(.data[[StudyID]]),
             .TIME = as.character(.data[[yinfor]]),
             .SITE = as.character(.data[[multi.sites]]),
             .COLOR = as.character(.data[[color_source_col]]),
             present = 1L) %>%
      select(.ID, .TIME, .SITE, .COLOR, present) %>%
      group_by(.ID, .TIME, .SITE) %>%
      summarise(
        .COLOR = {
          x <- .COLOR[!is.na(.COLOR) & .COLOR != ""]
          if (length(x) == 0) NA_character_ else x[1]
        },
        present = as.integer(any(present == 1L)),
        .groups = "drop"
      )

    if (color_by_fill) {
      fill_levels_raw <- unique(na.omit(df_ms$.COLOR))
      if (!is.null(orders)) {
        fill_in_user <- intersect(orders, fill_levels_raw)
        fill_levels <- c(fill_in_user, setdiff(fill_levels_raw, fill_in_user))
      } else {
        fill_levels <- fill_levels_raw
      }
      if (length(fill_levels) == 0) {
        color_by_fill <- FALSE
        fill_levels <- yinfor.order
      }
    } else {
      fill_levels <- yinfor.order
    }
    if (is.null(mycols)) {
      mycols <- scales::hue_pal()(length(fill_levels))
    } else if (length(mycols) < length(fill_levels)) {
      stop("Length of 'mycols' is smaller than the number of color categories.")
    }
    if (!is.null(names(mycols)) && all(fill_levels %in% names(mycols))) {
      fill_cols <- mycols[fill_levels]
    } else {
      fill_cols <- setNames(mycols[seq_along(fill_levels)], fill_levels)
    }

    data_long <- df_ms %>%
      tidyr::complete(.ID,
                      .TIME = yinfor.order,
                      .SITE = site_levels,
                      fill = list(present = 0L, .COLOR = NA_character_))
    data_long$present <- as.integer(data_long$present)
    data_long$.COLOR <- as.character(data_long$.COLOR)

    if (color_by_fill) {
      data_wide <- data_long %>%
        mutate(.COLOR = ifelse(present == 1L, .COLOR, NA_character_)) %>%
        pivot_wider(id_cols = .ID,
                    names_from = c(.TIME, .SITE),
                    values_from = .COLOR,
                    names_sep = "___",
                    values_fn = function(x) {
                      y <- x[!is.na(x)]
                      if (length(y) == 0) NA_character_ else as.character(y[1])
                    })

      pat_cols <- setdiff(names(data_wide), ".ID")
      if (isTRUE(pattern)) {
        unique_patterns <- data_wide %>%
          group_by(across(all_of(pat_cols))) %>%
          summarise(pattern_count = n(), .groups = "drop") %>%
          mutate(ID = row_number())
      } else {
        unique_patterns <- data_wide
        sample_mat <- as.data.frame(unique_patterns[, pat_cols, drop = FALSE], stringsAsFactors = FALSE)
        unique_patterns$pattern_count <- rowSums(!is.na(sample_mat))
        unique_patterns$ID <- as.character(unique_patterns$.ID)
      }

      data_long <- unique_patterns %>%
        pivot_longer(cols = all_of(pat_cols),
                     names_to = "panel_key",
                     values_to = "fill_key") %>%
        tidyr::separate(panel_key, into = c(".TIME", ".SITE"), sep = "___", remove = TRUE)
      data_long$present <- ifelse(is.na(data_long$fill_key), 0L, 1L)
    } else {
      data_wide <- data_long %>%
        pivot_wider(id_cols = .ID,
                    names_from = c(.TIME, .SITE),
                    values_from = present,
                    names_sep = "___",
                    values_fn = max)

      pat_cols <- setdiff(names(data_wide), ".ID")
      if (isTRUE(pattern)) {
        unique_patterns <- data_wide %>%
          group_by(across(all_of(pat_cols))) %>%
          summarise(pattern_count = n(), .groups = "drop") %>%
          mutate(ID = row_number())
      } else {
        unique_patterns <- data_wide
        sample_mat <- as.data.frame(unique_patterns[, pat_cols, drop = FALSE], stringsAsFactors = FALSE)
        unique_patterns$pattern_count <- rowSums(sample_mat > 0)
        unique_patterns$ID <- as.character(unique_patterns$.ID)
      }

      data_long <- unique_patterns %>%
        pivot_longer(cols = all_of(pat_cols),
                     names_to = "panel_key",
                     values_to = "present") %>%
        tidyr::separate(panel_key, into = c(".TIME", ".SITE"), sep = "___", remove = TRUE)
      data_long$fill_key <- ifelse(data_long$present == 1L, as.character(data_long$.TIME), NA_character_)
    }

    data_long$.TIME <- factor(data_long$.TIME, levels = yinfor.order)
    data_long$.SITE <- factor(data_long$.SITE, levels = site_levels)
    desired_order <- unique(data_long$ID)
    data_long$ID <- factor(data_long$ID, levels = rev(desired_order))
    unique_patterns$ID <- factor(unique_patterns$ID, levels = rev(desired_order))
    if (!isTRUE(pattern)) {
      id_levels <- levels(data_long$ID)
      n_ids <- length(id_levels)
      ncol_split <- if (n_ids <= 1) 1L else 2L
      ids_per_col <- ceiling(n_ids / ncol_split)
      id_pos <- match(as.character(data_long$ID), id_levels)
      data_long$panel_col <- factor(ifelse(id_pos <= ids_per_col, "Col 1", "Col 2"),
                                    levels = if (ncol_split == 1L) "Col 1" else c("Col 1", "Col 2"))
    }

    # 시점 내부에서 site 순서대로 나란히 배치
    n_site <- length(site_levels)
    # 시점 내 네모 간격
    compact_step <- site.compact.step
    compact_span <- compact_step * (n_site - 1)
    site_offset <- setNames(seq(-compact_span / 2, compact_span / 2, length.out = n_site), site_levels)
    data_long$xnum <- as.numeric(data_long$.TIME) + site_offset[as.character(data_long$.SITE)]
    data_long$fill_key <- factor(data_long$fill_key, levels = fill_levels)

    if (isTRUE(pattern)) {
      p1 <- ggplot(data_long, aes(x = xnum, y = ID, group = interaction(ID, .SITE))) +
        geom_line(color = "grey80", linewidth = 0.5) +
        geom_point(aes(fill = fill_key),
                   shape = 22, color = "black", size = 3.1, stroke = 0.5) +
        scale_fill_manual(values = fill_cols, breaks = fill_levels, drop = FALSE, na.translate = FALSE) +
        guides(fill = guide_legend(order = 1, override.aes = list(shape = 22, size = 3.2))) +
        scale_x_continuous(breaks = seq_along(yinfor.order),
                           labels = yinfor.order,
                           position = "top",
                           expand = ggplot2::expansion(mult = c(0.04, 0.04))) +
        labs(y = "Patterns", x = NULL) +
        theme_classic(base_size = 11) +
        theme(axis.text.x  = element_text(hjust = 0.5),
              legend.title = element_blank(),
              legend.position = "right")
    } else {
      p1 <- ggplot(data_long, aes(x = .TIME, y = ID)) +
        geom_line(aes(group = ID), color = "grey80", linewidth = 0.5) +
        geom_point(aes(fill = fill_key, group = .SITE),
                   shape = 22,
                   position = ggplot2::position_dodge(width = site.dodge.width),
                   color = "black", size = 3.8, stroke = 0.5) +
        scale_fill_manual(values = fill_cols, breaks = fill_levels, drop = FALSE, na.translate = FALSE) +
        guides(fill = guide_legend(order = 1, override.aes = list(shape = 22, size = 3.8))) +
        scale_x_discrete(position = "top") +
        labs(y = "Patterns", x = NULL) +
        theme_classic(base_size = 11) +
        theme(axis.text.x  = element_text(hjust = 0.5),
              axis.text.y  = element_text(size = 5),
              legend.title = element_blank(),
              legend.position = "right")
    }
    if (!isTRUE(pattern)) {
      p1 <- p1 +
        ggplot2::facet_wrap(ggplot2::vars(panel_col), nrow = 1, scales = "free_y") +
        theme(strip.background = element_blank(),
              strip.text.x = element_blank())
    }
  }

  # p3
  unique_patterns$pattern_count <- as.numeric(unique_patterns$pattern_count)
  count_label <- if (isTRUE(pattern)) "Count of Patients" else "Count of Samples"
  label_pos <- 0.4
  p3 <- ggplot(unique_patterns, aes(x = ID, y = pattern_count)) +
    geom_col(fill = "#149BEDFF", width = 0.85) +
    geom_text(aes(y = label_pos, label = pattern_count), hjust = 0, vjust = 0.5, size = 3) +
    coord_flip() +
    labs(y = count_label, x = NULL) +
    theme_classic(base_size = 11) +
    theme(axis.text.y  = element_blank(),
          axis.title.y = element_blank())

  fill_label_clean <- if (isTRUE(auto_fill_mode)) yinfor else fillinfor
  plot_label <- if (is.null(multi.sites)) fill_label_clean else multi.sites
  file_label <- if (is.null(multi.sites)) {
    fill_label_clean
  } else if (!isTRUE(auto_fill_mode) && !is.null(fillinfor) && fillinfor != "") {
    paste0(multi.sites, ".", fillinfor)
  } else {
    multi.sites
  }
  file_label <- paste0(file_label, if (isTRUE(pattern)) ".patternT" else ".patternF")
  file_label <- gsub("[^A-Za-z0-9._-]", "_", file_label)

  # 합성
  n_for_title <- if (isTRUE(pattern)) sum(unique_patterns$pattern_count, na.rm = TRUE) else nrow(unique_patterns)
  fill_title_tag <- if (!is.null(multi.sites) && !isTRUE(auto_fill_mode) && !is.null(fillinfor_user) && fillinfor_user != "") {
    paste0(" | fill: ", fillinfor_user)
  } else {
    ""
  }
  subtitle_auto <- if (!is.null(multi.sites) && !is.null(site_levels_for_subtitle) && length(site_levels_for_subtitle) > 0) {
    paste0("Square order (left-to-right within each visit): ", paste(site_levels_for_subtitle, collapse = ", "))
  } else {
    NULL
  }
  subtitle_parts <- c(name_user, subtitle, subtitle_auto)
  subtitle_parts <- subtitle_parts[!is.na(subtitle_parts) & nzchar(subtitle_parts)]
  subtitle_txt <- if (length(subtitle_parts) == 0) NULL else paste(subtitle_parts, collapse = "\n")
  title_txt <- sprintf(
    "%s%s (n = %s)",
    plot_label,
    fill_title_tag,
    n_for_title
  )


  if (isTRUE(pattern)) {
    p_final <- (p1 + p3) +
      plot_layout(widths = c(3, 1), guides = "collect") +
      plot_annotation(title = title_txt,
                      subtitle = subtitle_txt,
                      theme = theme(plot.title = element_text(hjust = 0.5, face = "bold"),
                                    plot.subtitle = element_text(hjust = 0.5, size = 10)))
  } else {
    p_final <- p1 +
      plot_annotation(title = title_txt,
                      subtitle = subtitle_txt,
                      theme = theme(plot.title = element_text(hjust = 0.5, face = "bold"),
                                    plot.subtitle = element_text(hjust = 0.5, size = 10)))
  }

  # 저장 (PDF)
  outfile <- file.path(pdf_dir, sprintf("pattern_plot.%s.%s.%s%s.pdf", file_label, project, name, stamp))
  if (isTRUE(patchwork)) return(invisible(p_final))
  ggsave(outfile, plot = p_final, device = cairo_pdf, width = width, height = height, units = "in", limitsize = FALSE)

  write.csv(unique_patterns, sprintf("%s/pattern_table.%s.%s.%s%s.csv", pattern_dir, file_label, project, name, stamp))
  write.csv(data_wide, sprintf("%s/pattern_wide.%s.%s.%s%s.csv", pattern_dir, file_label, project, name, stamp))

  invisible(list(plot = p_final, pattern_table = unique_patterns, long_data = data_long, file = outfile))
}
