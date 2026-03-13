#' Extended Barplot from Phyloseq Data
#'
#' This function generates extended bar plots from a Phyloseq object for comparative analysis
#' between two groups at specified taxonomic or functional levels. It calculates statistical
#' significance using the Wilcoxon test and includes options for multiple transformations and
#' confidence interval calculations.
#'
#' @param psIN Phyloseq object containing microbiome data for analysis.
#' @param project Name of the project or analysis to label outputs.
#' @param group1 First group for comparison, specified as a level within the metadata variable.
#' @param group2 Second group for comparison, opposed to group1.
#' @param name Optional name for saving the plot, influences output file name.
#' @param mvar Metadata variable in Phyloseq object to divide data for comparison.
#' @param func Variable for dynamic grouping based on taxonomy or functional data ('KO', 'KO.des', 'pathway', 'path.des').
#' @param wilcox.p Cutoff for the p-value from the Wilcoxon test to determine statistical significance.
#' @param transform Normalization/transformation method. One of "clr", "relative", "tss".
#' @param p_adjust_method Multiple-testing correction method passed to \code{stats::p.adjust}.
#' @param use_adjusted_p Logical; if TRUE, filter/significance uses adjusted p-values (q-values).
#' @param height Height of the output plot in inches.
#' @param width Width of the output plot in inches.
#'
#' @return Invisibly returns a list containing plot object, Wilcoxon table, and plotted data.
#'
#' @details
#' The function preprocesses the data to perform normalization transformations such as
#' Total Sum Scaling (TSS), relative abundance calculation, and Centered Log Ratio (CLR).
#' It then computes Wilcoxon tests for the specified groups across the data grouped by
#' 'func' criteria. Results are visualized in an extended bar plot showing mean differences,
#' confidence intervals, and p-values for significant differences.
#'
#' @examples
#' # psIN is a Phyloseq object
#' # Example usage:
#' Go_extendedBarplot(psIN = psIN,
#'                    project = "MicrobiomeAnalysis",
#'                    group1 = "Control",
#'                    group2 = "Treatment",
#'                    name = "ExtendedBarplotExample",
#'                    mvar = "Group",
#'                    func = "pathway",
#'                    wilcox.p = 0.05,
#'                    height = 8,
#'                    width = 11)
#'
#' @export

Go_extendedBarplot <- function(psIN,
                               project,
                               group1,
                               group2,
                               name=NULL,
                               mvar,
                               func, #  "KO", "KO.des", "pathway", "path.des"
                               wilcox.p = 0.05,
                               transform = c("clr", "relative", "tss"),
                               p_adjust_method = "BH",
                               use_adjusted_p = TRUE,
                               height,
                               width){

  if(!is.null(dev.list())) dev.off()
  # out dir
  out <- file.path(sprintf("%s_%s",project, format(Sys.Date(), "%y%m%d")))
  if(!file_test("-d", out)) dir.create(out, recursive = TRUE)
  out_path <- file.path(sprintf("%s_%s/pdf",project, format(Sys.Date(), "%y%m%d")))
  if(!file_test("-d", out_path)) dir.create(out_path, recursive = TRUE)
  out_table_path <- file.path(sprintf("%s_%s/table/extendedBarplot",project, format(Sys.Date(), "%y%m%d")))
  if(!file_test("-d", out_table_path)) dir.create(out_table_path, recursive = TRUE)

  # input

  map <- data.frame(phyloseq::sample_data(psIN))
  if(!mvar %in% names(map)) {
    stop("mvar is not a valid column in the sample data.")
  }
  if (!group1 %in% map[[mvar]] || !group2 %in% map[[mvar]]) {
    stop("group1/group2 are not both present in sample_data(psIN)[[mvar]].")
  }
  if (!is.numeric(wilcox.p) || length(wilcox.p) != 1 || wilcox.p <= 0 || wilcox.p >= 1) {
    stop("wilcox.p must be a single numeric value between 0 and 1.")
  }
  if (!is.numeric(height) || !is.numeric(width) || height <= 0 || width <= 0) {
    stop("height and width must be positive numeric values.")
  }

  p_adjust_method <- match.arg(
    p_adjust_method,
    c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")
  )

  # Keep only the two groups to compare.
  samples_to_keep <- phyloseq::sample_data(psIN)[[mvar]] %in% c(group1, group2)
  if (sum(samples_to_keep, na.rm = TRUE) < 4) {
    stop("Not enough samples after filtering to group1/group2.")
  }
  ps1.sel <- phyloseq::prune_samples(samples_to_keep, psIN)


  df <- phyloseq::psmelt(ps1.sel)
  if (!func %in% names(df)) {
    stop(sprintf("func '%s' is not a valid column in psmelt(psIN).", func))
  }
  if (!mvar %in% names(df)) {
    stop(sprintf("mvar '%s' is not found in psmelt(psIN).", mvar))
  }

  transform <- match.arg(transform)

  if (transform == "tss") {
    df_top_norm <- df %>%
      dplyr::group_by(Sample) %>%
      dplyr::mutate(Abundance = Abundance / sum(Abundance, na.rm = TRUE) * 1e6) %>%
      dplyr::ungroup()
  } else if (transform == "relative") {
    df_top_norm <- df %>%
      dplyr::group_by(Sample) %>%
      dplyr::mutate(Abundance = Abundance / sum(Abundance, na.rm = TRUE)) %>%
      dplyr::ungroup()
  } else if (transform == "clr") {
    df_top_norm <- df %>%
      dplyr::group_by(Sample) %>%
      dplyr::mutate(
        Rel_Abundance = (Abundance + 1) / sum(Abundance + 1, na.rm = TRUE),
        Log_Abundance = log(Rel_Abundance),
        Abundance = Log_Abundance - mean(Log_Abundance, na.rm = TRUE)
      ) %>%
      dplyr::ungroup()
  }



  #===== wilcox test
  func_map <- list(
    pathway = c("pathway", "path.des"),
    path.des = c("pathway", "path.des"),
    KO = c("KO", "KO.des"),
    KO.des = c("KO", "KO.des"),
    EC = c("KO", "EC"),
    ec_number = c("KO", "ec_number"),
    Function = c("Function", NA_character_)
  )

  if (!func %in% names(func_map)) {
    stop("func must be one of: pathway, path.des, KO, KO.des, EC, ec_number, Function")
  }
  func1 <- func_map[[func]][1]
  func2 <- func_map[[func]][2]
  if (!is.na(func2) && !func2 %in% names(df_top_norm)) {
    func2 <- NULL
  }

  # Function to perform Wilcoxon tests with dynamic group names
  perform_wilcox_test <- function(data, group1, group2, func1, func2 = NULL, mvar) {
    grouped <- data %>%
      dplyr::group_by(!!rlang::sym(func1)) %>%
      dplyr::summarise(
        list_group1 = list(Abundance[.data[[mvar]] == group1]),
        list_group2 = list(Abundance[.data[[mvar]] == group2]),
        count_group1 = length(Abundance[.data[[mvar]] == group1]),
        count_group2 = length(Abundance[.data[[mvar]] == group2]),
        .groups = 'drop'
      ) %>%
      rowwise() %>%
      dplyr::mutate(
        p_value = ifelse(
          count_group1 > 1 & count_group2 > 1,
          stats::wilcox.test(unlist(list_group1), unlist(list_group2), exact = FALSE)$p.value,
          NA_real_
        )
      ) %>%
      ungroup()

    grouped <- grouped %>%
      dplyr::mutate(p_adj = stats::p.adjust(p_value, method = p_adjust_method))

    if (!is.null(func2) && !is.na(func2) && func2 != func1) {
      desc_df <- data %>% dplyr::select(!!rlang::sym(func1), !!rlang::sym(func2)) %>% dplyr::distinct()
      grouped <- dplyr::left_join(grouped, desc_df, by = func1)
    }

    if (!is.null(func2) && !is.na(func2) && func2 != func1) {
      grouped %>% dplyr::select(!!rlang::sym(func1), !!rlang::sym(func2), p_value, p_adj, count_group1, count_group2)
    } else {
      grouped %>% dplyr::select(!!rlang::sym(func1), p_value, p_adj, count_group1, count_group2)
    }
  }

  wilcox_results <- perform_wilcox_test(df_top_norm, group1, group2,func1, func2, mvar)

  #===== calculate quantile interval
  df_summary <- df_top_norm %>%
    dplyr::group_by(!!rlang::sym(func), !!rlang::sym(mvar)) %>%
    dplyr::summarise(
      Mean_Abundance = mean(Abundance, na.rm = TRUE),
      Lower_CI = quantile(Abundance, probs = 0.025, na.rm = TRUE),
      Upper_CI = quantile(Abundance, probs = 0.975, na.rm = TRUE),
      .groups = 'drop'
    )


  df_wide <- df_summary %>%
    tidyr::pivot_wider(
      names_from = !!rlang::sym(mvar),
      values_from = c(Mean_Abundance, Lower_CI, Upper_CI),
      names_sep = "_"
    )

  #===== Calculate the difference in means and the confidence intervals
  mean_col_group1 <- paste("Mean_Abundance", group1, sep = "_")
  mean_col_group2 <- paste("Mean_Abundance", group2, sep = "_")
  lower_col_group1 <- paste("Lower_CI", group1, sep = "_")
  lower_col_group2 <- paste("Lower_CI", group2, sep = "_")
  upper_col_group1 <- paste("Upper_CI", group1, sep = "_")
  upper_col_group2 <- paste("Upper_CI", group2, sep = "_")


  df_wide <- df_wide %>%
    dplyr::mutate(
      Diff_Mean = !!rlang::sym(mean_col_group1) - !!rlang::sym(mean_col_group2),
      Diff_Lower = !!rlang::sym(lower_col_group1) - !!rlang::sym(upper_col_group2),
      Diff_Upper = !!rlang::sym(upper_col_group1) - !!rlang::sym(lower_col_group2)
    )

  df_wide <- df_wide %>%
    dplyr::filter(
      !is.na(!!rlang::sym(mean_col_group1)),
      !is.na(!!rlang::sym(mean_col_group2))
    )

  if (nrow(df_wide) == 0) {
    stop("No comparable features found for both groups after aggregation.")
  }

  data_long <- df_wide %>%
    tidyr::pivot_longer(
      cols = c(
        paste("Mean_Abundance_", group1, sep = ""),
        paste("Mean_Abundance_", group2, sep = ""),
        paste("Lower_CI_", group1, sep = ""),
        paste("Upper_CI_", group1, sep = ""),
        paste("Lower_CI_", group2, sep = ""),
        paste("Upper_CI_", group2, sep = "")
      ),
      names_to = c(".value", mvar),
      names_pattern = sprintf("(.+)_(%s|%s)$", group1, group2)
    )

  data_long <- data_long %>%
    dplyr::mutate(
      Higher_Group = ifelse(Diff_Mean >= 0, group1, group2)
    ) %>%
    ungroup()

  ordered_feature <- df_wide %>%
    dplyr::arrange(Diff_Mean) %>%
    dplyr::pull(!!rlang::sym(func))
  data_long$Reordered_Category <- factor(data_long[[func]], levels = unique(ordered_feature))


  data_long$Panel <- "Difference Mean with CI"
  point_data <- data_long
  point_data$Panel <- "Mean Proportion"

  pvalue_data <- data_long
  pvalue_panel_label <- if (isTRUE(use_adjusted_p) && p_adjust_method == "BH") {
    "BH-adjusted q_value"
  } else if (isTRUE(use_adjusted_p)) {
    sprintf("Adjusted q_value (%s)", p_adjust_method)
  } else {
    "p_value"
  }
  pvalue_data$Panel <- pvalue_panel_label

  plot_data <- rbind(data_long, point_data,pvalue_data)

  plot_data$Panel <- factor(plot_data$Panel, levels = c("Mean Proportion", "Difference Mean with CI", pvalue_panel_label))

  join_col <- if (func %in% names(wilcox_results)) func else func1
  merged_data <- dplyr::left_join(plot_data, wilcox_results, by = setNames(join_col, join_col))

  filter_col <- if (isTRUE(use_adjusted_p)) "p_adj" else "p_value"
  merged_data.sig <- merged_data %>%
    dplyr::filter(!is.na(.data[[filter_col]]) & .data[[filter_col]] < wilcox.p)
  merged_data.sig$stat_label <- if (isTRUE(use_adjusted_p)) {
    sprintf("q = %.3f", merged_data.sig$p_adj)
  } else {
    sprintf("p = %.3f", merged_data.sig$p_value)
  }

  if (nrow(merged_data.sig) == 0) {
    message(sprintf("No feature passed %s < %.3f.", filter_col, wilcox.p))
  }

  merged_data.sig[[mvar]] <- factor(merged_data.sig[[mvar]] , levels = c(group1, group2))


  # Save result tables for reproducibility.
  utils::write.csv(
    wilcox_results,
    file.path(out_table_path,
              sprintf("wilcox.%s.vs.%s.%s.%s.csv", group1, group2, project, format(Sys.Date(), "%y%m%d"))),
    row.names = FALSE
  )
  utils::write.csv(
    merged_data.sig,
    file.path(out_table_path,
              sprintf("plot_data_sig.%s.vs.%s.%s.%s.csv", group1, group2, project, format(Sys.Date(), "%y%m%d"))),
    row.names = FALSE
  )

  if (nrow(merged_data.sig) > 0) {
    p <- ggplot2::ggplot(merged_data.sig, ggplot2::aes_string(x = "Reordered_Category", y = "Mean_Abundance", fill = mvar))

    p1 <- p + ggplot2::geom_bar(data = subset(merged_data.sig, Panel == "Mean Proportion"),
                                stat = "identity", position = ggplot2::position_dodge(width = 0.7), width = 0.6)

    p2 <- p1 + ggplot2::geom_point(data = subset(merged_data.sig, Panel == "Difference Mean with CI"),
                                   ggplot2::aes_string(x = "Reordered_Category", y = "Diff_Mean", group = mvar, fill = "Higher_Group"),
                                   size = 3, shape = 23) +
      ggplot2::geom_errorbar(data = subset(merged_data.sig, Panel == "Difference Mean with CI"),
                             ggplot2::aes_string(ymin = "Diff_Lower", ymax = "Diff_Upper", group = mvar), width = 0.25)

    p3 <- p2 + ggplot2::geom_text(data = subset(merged_data.sig, Panel == pvalue_panel_label),
                                  ggplot2::aes(x = Reordered_Category, y = 0, label = stat_label),
                                  hjust = 0, vjust = 0, inherit.aes = FALSE, size = 3, color= "darkgrey")

    p4 <- p3 + ggplot2::facet_grid(. ~ Panel, scales = "free_x") +
      ggplot2::coord_flip() +
      ggplot2::scale_fill_manual(values = c("blue", "orange")) +
      ggplot2::labs(x = func, y = NULL, fill = "Group") +
      ggplot2::theme_minimal() +
      ggplot2::theme(legend.position = "bottom",
                     legend.title = ggplot2::element_blank())
  } else {
    p4 <- ggplot2::ggplot() +
      ggplot2::annotate("text", x = 0, y = 0,
                        label = sprintf("No feature passed %s < %.3f", filter_col, wilcox.p),
                        size = 4) +
      ggplot2::theme_void()
  }



  pdf(sprintf("%s/extended_error_barplot.(%s.vs.%s%s).%s.%s.pdf", out_path,
              group1,
              group2,
              ifelse(is.null(name), "", paste0(".", gsub("^\\(|\\)$", "", name))),
              project,
              format(Sys.Date(), "%y%m%d")), height = height, width = width)

  print(p4)

  dev.off()
  message("Recommendation: Use dedicated DA methods (e.g., ALDEx2, ANCOM-BC2, MaAsLin2) for primary statistical inference.")
  invisible(list(
    plot = p4,
    wilcox_results = wilcox_results,
    merged_data_significant = merged_data.sig,
    transform = transform,
    p_adjust_method = p_adjust_method,
    use_adjusted_p = use_adjusted_p
  ))
}
