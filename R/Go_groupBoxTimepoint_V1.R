#' Go_groupBoxTimepoint
#'
#' Create grouped box plots across timepoints with optional baseline subtraction.
#'
#' This function generates grouped box plots with overlaid lines and points
#' representing group means across multiple timepoints. It supports baseline
#' subtraction for paired data, allowing visualization of changes relative
#' to a specified baseline timepoint.
#'
#' @param df A data frame containing the data to be plotted.
#' @param project A character string specifying the project name. Used for
#'   file naming and directory structure.
#' @param name An optional character string for additional naming in the output file.
#' @param maingroup A character string specifying the column in `df` defining
#'   the main grouping variable (e.g., treatment or condition).
#' @param timepoint A character string specifying the column in `df` defining
#'   the timepoints.
#' @param outcomes A character vector specifying the outcome variables to be plotted.
#'   These columns must exist in `df`. If `baseline` and `paired` are specified,
#'   subtraction is automatically applied to all existing outcomes.
#' @param mycols An optional named character vector of colors to be used for
#'   different groups in `maingroup`. If NULL, default ggplot colors are used.
#' @param orders An optional character vector specifying the order of timepoints
#'   and/or groups in the plots. If NULL, the default order in the data is used.
#' @param fill_labels An optional character vector providing legend labels
#'   corresponding to `mycols`.
#' @param title An optional character string specifying the plot title.
#'   If NULL, the main group name and test method are used as the title.
#' @param x.axis An optional character string specifying the x-axis label.
#'   If NULL, the `timepoint` names are used.
#' @param height Numeric value specifying the height (in inches) of the output PDF.
#'   Default is 3.5.
#' @param width Numeric value specifying the width (in inches) of the output PDF.
#'   Default is 5.
#' @param baseline An optional character string specifying the baseline
#'   timepoint (e.g., "Baseline"). If provided together with `paired`,
#'   outcome values at this baseline are subtracted from all corresponding
#'   values in the same paired subject.
#' @param paired An optional character string specifying the subject ID column
#'   used to match baseline values (e.g., "SubjectID"). Required for
#'   baseline subtraction.
#'
#' @return
#' Saves box plots for each outcome variable as PDF files in the directory
#' returned by `Go_path(project, pdf="yes")`. Each plot is also printed
#' to the R graphical device.
#'
#' @details
#' For each variable in `outcomes`, the function:
#' \itemize{
#'   \item Creates box plots grouped by `maingroup` and `timepoint`.
#'   \item Overlays mean values with connecting lines and points.
#'   \item Performs Wilcoxon tests between groups and annotates p-values.
#'   \item (If `baseline` and `paired` are specified) Calculates difference-from-baseline
#'         values for each outcome. New columns are created in the format:
#'         \code{<outcome>_D0} for baseline values and
#'         \code{Sub0_<outcome>} for subtracted values.
#' }
#'
#' Only existing outcomes in the dataframe are processed; missing outcomes
#' are automatically skipped.
#'
#' @importFrom dplyr group_by summarise ungroup
#' @importFrom ggplot2 ggplot aes_string geom_boxplot geom_line geom_point theme_bw labs scale_fill_manual scale_colour_manual ggtitle
#' @importFrom ggpubr stat_compare_means
#'
#' @examples
#' \dontrun{
#' df <- data.frame(
#'   SubjectID = rep(1:10, each = 3),
#'   Time = rep(c("Baseline", "Day5", "Day10"), times = 10),
#'   Group = rep(c("Control", "Treatment"), each = 15),
#'   Chao1 = rnorm(30, 50, 10),
#'   Shannon = rnorm(30, 3, 0.5)
#' )
#'
#' # Without baseline subtraction
#' Go_groupBoxTimepoint(
#'   df = df,
#'   project = "IBS_Study",
#'   name = "AlphaDiversity",
#'   maingroup = "Group",
#'   timepoint = "Time",
#'   outcomes = c("Chao1", "Shannon"),
#'   mycols = c("Control" = "#1f77b4", "Treatment" = "#d62728"),
#'   orders = c("Baseline", "Day5", "Day10"),
#'   title = "Alpha Diversity"
#' )
#'
#' # With baseline subtraction (e.g., delta from Baseline)
#' Go_groupBoxTimepoint(
#'   df = df,
#'   project = "IBS_Study",
#'   name = "AlphaDiversity_Delta",
#'   maingroup = "Group",
#'   timepoint = "Time",
#'   outcomes = c("Chao1", "Shannon"),
#'   mycols = c("Control" = "#1f77b4", "Treatment" = "#d62728"),
#'   orders = c("Baseline", "Day5", "Day10"),
#'   baseline = "Baseline",
#'   paired = "SubjectID",
#'   title = "Alpha Diversity (Δ from Baseline)"
#' )
#' }
#'
#' @export

Go_groupBoxTimepoint <- function(df = NULL,
                                 project = NULL,
                                 name = NULL,
                                 maingroup = NULL,
                                 timepoint = NULL,
                                 outcomes = NULL,
                                 mycols = NULL,
                                 orders = NULL,
                                 fill_labels = NULL,
                                 title = NULL,
                                 x.axis = NULL,
                                 height = NULL,
                                 width = NULL,
                                 baseline = NULL,   # 추가
                                 paired = NULL) {   # 추가

  library(dplyr)
  library(ggplot2)
  library(ggpubr)

  if (is.null(df)) stop("Input dataframe (df) is required.")
  if (is.null(project)) stop("Project name is required.")

  df.sel <- df

  if (is.null(mycols)) {
    print("mycols is not defined. Default colors will be used.")
  }
  if (is.null(orders)) {
    print("orders is not defined.")
  }

  if (is.null(height)) height <- 3.5
  if (is.null(width)) width <- 5

  #--------------------------------------------
  # baseline 보정(subtraction) 수행
  #--------------------------------------------
  if (!is.null(baseline) && !is.null(paired)) {
    if (baseline %in% unique(df.sel[[timepoint]])) {
      message(sprintf("Performing baseline subtraction using baseline = '%s' and paired = '%s'", baseline, paired))

      df.base <- subset(df.sel, df.sel[[timepoint]] == baseline)
      matched_indices <- match(df.sel[[paired]], df.base[[paired]])

      # outcomes 중 실제 존재하는 것만 사용
      outcomes.exist <- outcomes[outcomes %in% colnames(df.sel)]

      for (var in outcomes.exist) {
        baseline_col <- paste0(var, "_D0")
        sub_col <- paste0("Sub0_", var)

        df.sel[[baseline_col]] <- df.base[[var]][matched_indices]
        df.sel[[sub_col]] <- df.sel[[var]] - df.sel[[baseline_col]]
      }

      # outcomes 갱신 (Sub0_ prefix 적용)
      outcomes <- paste0("Sub0_", outcomes.exist)
    } else {
      warning("Baseline value not found in timepoint column. Skipping baseline subtraction.")
    }
  }

  #--------------------------------------------
  # 변수별 boxplot 작성
  #--------------------------------------------
  for (variable in outcomes) {
    if (!is.null(dev.list())) dev.off()

    if (!variable %in% colnames(df.sel)) {
      warning(sprintf("Variable '%s' not found in dataframe. Skipping.", variable))
      next
    }

    average_data <- df.sel %>%
      dplyr::group_by(!!sym(timepoint), !!sym(maingroup)) %>%
      dplyr::summarise(
        !!paste0("mean_", variable) := mean(!!sym(variable), na.rm = TRUE),
        .groups = "drop"
      ) %>%
      ungroup()

    if (!is.null(orders)) {
      df.sel[, timepoint] <- factor(df.sel[, timepoint], levels = intersect(orders, df.sel[, timepoint]))
      df.sel[, maingroup] <- factor(df.sel[, maingroup], levels = intersect(orders, df.sel[, maingroup]))
    }

    max_label_y <- if (all(is.na(df.sel[[variable]]))) 1 else max(df.sel[[variable]], na.rm = TRUE) + 0.2

    p <- ggplot(df.sel, aes_string(x = timepoint, y = variable, fill = maingroup)) +
      geom_boxplot(outlier.shape = NA, lwd = 0.3) +
      geom_line(data = average_data, aes_string(x = timepoint, y = sprintf("mean_%s", variable), group = maingroup, colour = maingroup), size = 1, alpha = 0.4) +
      geom_point(aes_string(colour = maingroup), position = position_dodge(width = 0.75), size = 1.5, alpha = 0.8) +
      theme_bw() +
      theme(legend.position = "bottom", legend.title = element_blank()) +
      labs(y = variable, x = NULL) +
      stat_compare_means(method = "wilcox.test", label = "p.format",
                         aes_string(group = maingroup),
                         label.y = max_label_y)

    if (!is.null(mycols)) {
      p <- p + scale_fill_manual(values = mycols, labels = fill_labels) +
        scale_colour_manual(values = mycols, labels = fill_labels)
    }

    if (!is.null(title)) {
      p <- p + ggtitle(sprintf("%s (%s)", title, "wilcox.test"))
    } else {
      p <- p + ggtitle(sprintf("%s (%s)", maingroup, "wilcox.test"))
    }

    dir <- Go_path(project, pdf = "yes", table = "no", path = NULL)

    pdf(sprintf("%s/groupBox3.%s.%s.%s%s%s.pdf", dir$pdf,
                project,
                maingroup,
                ifelse(is.null(variable), "", paste(variable, ".", sep = "")),
                ifelse(is.null(name), "", paste(name, ".", sep = "")),
                format(Sys.Date(), "%y%m%d")), height = height, width = width)

    print(p)
    dev.off()
  }
}
