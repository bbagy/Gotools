#' Go_groupBoxTimepoint
#'
#' This function creates grouped box plots with overlaid points and lines representing averages for different timepoints and main groups. It allows customization of plot appearance, including colors and dimensions.
#'
#' @param df A dataframe containing the data to be plotted.
#' @param project A character string specifying the project name, used for file naming and directory structure.
#' @param name An optional character string for additional naming in the output file.
#' @param maingroup A character string specifying the column name in `df` that defines the main groups.
#' @param timepoint A character string specifying the column name in `df` that defines the timepoints.
#' @param outcomes A character vector specifying the names of the outcomes to be plotted.
#' @param mycols An optional character vector of colors to be used for the main groups; if NULL, default colors are used.
#' @param orders An optional character vector specifying the order of groups and/or timepoints; if NULL, the default order is used.
#' @param fill_labels An optional character vector providing labels for the legend corresponding to `mycols`.
#' @param title An optional character string specifying the title of the plot; if NULL, the main group name is used as the title.
#' @param x.axis An optional character string specifying the label for the x-axis; if NULL, timepoint names are used.
#' @param height An optional numeric value specifying the height of the output PDF in inches; default is 3.5 inches.
#' @param width An optional numeric value specifying the width of the output PDF in inches; default is 5 inches.
#'
#' @return Saves the box plots as PDF files in the specified directory and prints each plot to the R graphical device.
#'
#' @details
#' The function iteratively processes each outcome specified in `outcomes`, creating a box plot for each. The box plots display data distribution across different groups and timepoints. Lines and points on the plots represent the average values for each group and timepoint. Statistical comparisons are made using the Wilcoxon test, and p-values are displayed on the plots.
#'
#' @importFrom dplyr group_by summarise ungroup
#' @importFrom ggplot2 ggplot aes_string geom_boxplot geom_line geom_point theme_bw labs scale_fill_manual scale_colour_manual ggtitle
#' @importFrom ggpubr stat_compare_means
#'
#' @examples
#' \dontrun{
#' df <- data.frame(
#'   Time = rep(c("T1", "T2", "T3"), each = 20),
#'   Group = rep(c("Control", "Treatment"), each = 30),
#'   Measurement = rnorm(60)
#' )
#' Go_groupBoxTimepoint(
#'   df = df,
#'   project = "ExampleProject",
#'   name = "ExamplePlot",
#'   maingroup = "Group",
#'   timepoint = "Time",
#'   outcomes = c("Measurement"),
#'   mycols = c("Control" = "blue", "Treatment" = "red"),
#'   orders = c("T1", "T2", "T3"),
#'   title = "Example Group Box Plot",
#'   height = 4,
#'   width = 6
#' )
#' }
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
