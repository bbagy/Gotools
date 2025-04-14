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
#' @param height Height of the output plot in inches.
#' @param width Width of the output plot in inches.
#'
#' @return Generates a PDF file containing the extended bar plot visualizing differences
#' between groups along with overlaid statistical significance. The function outputs
#' the modified data with significance levels and confidence intervals as well.
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
                               height,
                               width){

  library(phyloseq)
  library(dplyr)
  library(ggplot2)
  library(tidyr)

  if(!is.null(dev.list())) dev.off()
  # out dir
  out <- file.path(sprintf("%s_%s",project, format(Sys.Date(), "%y%m%d")))
  if(!file_test("-d", out)) dir.create(out)
  out_path <- file.path(sprintf("%s_%s/pdf",project, format(Sys.Date(), "%y%m%d")))
  if(!file_test("-d", out_path)) dir.create(out_path)

  # input

  map <- data.frame(sample_data(psIN))
  if(!mvar %in% names(map)) {
    stop("mvar is not a valid column in the sample data.")
  }

  # ps1.sel <- subset_samples(psIN, map[[mvar]] %in% c(group1, group2)) # it doesn't work.
  # ps1.sel <- subset_samples(psIN, get_variable(psIN, mvar) %in% c(group1, group2))
  # ps1.sel <- subset_samples(psIN, map[,mvar] %in% c(group1,group2));ps1.sel

  # Manually define a logical condition function
  manual_in <- function(x, y) {
    sapply(x, function(xi) any(xi == y))
  }

  # Apply this in the subset_samples context
  samples_to_keep <- manual_in(sample_data(psIN)[[mvar]], c(group1, group2))
  ps1.sel <- prune_samples(samples_to_keep, psIN)


  df <- psmelt(ps1.sel)

  #===TSS
  df_top_norm <- df %>%
    dplyr::group_by(Sample) %>%
    dplyr::mutate(Abundance = Abundance / sum(Abundance) * 1e6)

  head(df_top_norm$Abundance)

  #===relative
  df_top_norm <- df %>%
    dplyr::group_by(Sample) %>%
    dplyr::mutate(Abundance = Abundance / sum(Abundance))



  #==== Centered Log Ratio (CLR) transformation
  df_top_norm <- df %>%
    dplyr::group_by(Sample) %>%
    dplyr::mutate(
      Log_Abundance = log(Abundance + 1),
      Geometric_Mean = exp(mean(Log_Abundance)),
      Abundance = Log_Abundance - log(Geometric_Mean)
    )
  head(df_top_norm)



  #===== wilcox test
  if(func == "pathway" | func == "path.des"){
    func1 <- "pathway"
    func2 <- "path.des"
  } else if(func == "KO" | func == "KO.des"){
    func1 <- "KO"
    func2 <- "KO.des"
  } else if(func == "EC"){
    func1 <- "KO"
    func2 <- "EC"
  }else if(func == "ec_number"){
    func1 <- "KO"
    func2 <- "ec_number"
  } else if(func == "Function"){
    func1 <- "Function"
  }

  # Function to perform Wilcoxon tests with dynamic group names
  perform_wilcox_test <- function(data, group1, group2, func1, func2 = NULL, mvar) {
    library(dplyr)
    library(tidyr)
    library(rlang)

    # 그룹화 대상: func1만 사용
    grouped <- data %>%
      dplyr::group_by(!!sym(func1)) %>%
      dplyr::summarise(
        list_group1 = list(Abundance[get(mvar) == group1]),
        list_group2 = list(Abundance[get(mvar) == group2]),
        count_group1 = length(Abundance[get(mvar) == group1]),
        count_group2 = length(Abundance[get(mvar) == group2]),
        .groups = 'drop'
      ) %>%
      rowwise() %>%
      mutate(
        p_value = if (count_group1 > 1 && count_group2 > 1) {
          stats::wilcox.test(unlist(list_group1), unlist(list_group2), exact = FALSE)$p.value
        } else {
          NA_real_
        }
      ) %>%
      ungroup()

    # 2. func2가 주어졌다면 description 붙이기
    if (!is.null(func2) && func2 != func1) {
      desc_df <- data %>% select(!!sym(func1), !!sym(func2)) %>% distinct()
      grouped <- left_join(grouped, desc_df, by = func1)
    }

    # 결과 정리
    if (!is.null(func2) && func2 != func1) {
      grouped %>% select(!!sym(func1), !!sym(func2), p_value, count_group1, count_group2)
    } else {
      grouped %>% select(!!sym(func1), p_value, count_group1, count_group2)
    }
  }


  # Example usage
  wilcox_results <- perform_wilcox_test(df_top_norm, group1, group2,func1, func2, mvar)

  # Print results to debug
  print(wilcox_results)





  #===== calculate confidance interval

  # Function to perform Wilcoxon tests with dynamic group names
  df_summary <- df_top_norm %>%
    dplyr::group_by(!!rlang::sym(func), !!rlang::sym(mvar)) %>%
    dplyr::summarise(
      Mean_Abundance = mean(Abundance, na.rm = TRUE),
      Lower_CI = quantile(Abundance, probs = 0.025),
      Upper_CI = quantile(Abundance, probs = 0.975),
      .groups = 'drop'
    )


  #===== table transform
  df_wide <- df_summary %>%
    pivot_wider(
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

  data_long <- df_wide %>%
    pivot_longer(
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

  colnames(data_long)


  data_long <- data_long %>%
    dplyr::group_by(!!rlang::sym(func)) %>%
    dplyr::mutate(Higher_Group = ifelse(Mean_Abundance[!!rlang::sym(mvar) == group1] > Mean_Abundance[!!rlang::sym(mvar) == group2], group1, group2)) %>%
    ungroup()


  data_long$Reordered_Category <- reorder(data_long[[func]], data_long$Diff_Mean)


  # Update 'data_long' to include a 'Panel' column for bar plots
  data_long$Panel <- "Difference Mean with CI"  # for points and their error bars
  # Create a separate data frame for the points and error bars
  point_data <- data_long
  point_data$Panel <- "Mean Proportion"  # for bar plots and their error bars

  pvalue_data <- data_long
  pvalue_data$Panel <- "p_value"

  # Combine the data for plotting
  plot_data <- rbind(data_long, point_data,pvalue_data)

  # Ensure that the Panel factor has the intended order
  plot_data$Panel <- factor(plot_data$Panel, levels = c("Mean Proportion", "Difference Mean with CI","p_value"))


  colnames(plot_data)
  head(wilcox_results)



  # Assuming wilcox_results and plot_data are already created and they share common identifiers 'pathway' and 'path.des'
  merged_data <- left_join(plot_data, wilcox_results, by = c(func))

  merged_data.sig <- subset(merged_data, p_value < wilcox.p);dim(merged_data.sig)



  merged_data.sig[mvar] <- factor(merged_data.sig[[mvar]] , levels = c(group1, group2))



  # Plotting
  p <- ggplot(merged_data.sig, aes_string(x = "Reordered_Category", y = "Mean_Abundance", fill = mvar))


  p1 <- p + geom_bar(data = subset(merged_data.sig, Panel == "Mean Proportion"),
                     stat = "identity", position = position_dodge(width = 0.7), width = 0.6)
   # geom_errorbar(data = subset(merged_data.sig, Panel == "Proportion with CI"),
   #                aes_string(ymin = "Lower_CI", ymax = "Upper_CI", group = mvar),
   #               position = position_dodge(width = 0.7), width = 0.25)

  p2 <- p1 + geom_point(data = subset(merged_data.sig, Panel == "Difference Mean with CI"),
                        aes_string(x = "Reordered_Category", y = "Diff_Mean", group = mvar, fill = "Higher_Group"),
                        size = 3, shape = 23) +
    geom_errorbar(data = subset(merged_data.sig, Panel == "Difference Mean with CI"),
                  aes_string(ymin = "Diff_Lower", ymax = "Diff_Upper", group = mvar), width = 0.25)



  # Adjusting the geom_text for p-values
  p3 <- p2 + geom_text(data = subset(merged_data.sig, Panel == "p_value"),
                       aes(x = Reordered_Category, y = 0, label = sprintf("p = %.3f", p_value)),
                       hjust = 0, vjust = 0, inherit.aes = FALSE, size = 3, color= "darkgrey")

  # Print the final plot with p-values

  p4 <- p3 + facet_grid(. ~ Panel, scales = "free_x") + #space "free_x"
    coord_flip() +
    scale_fill_manual(values = c("blue", "orange")) +
    labs(x = func, y = NULL, fill = "Group") +
    theme_minimal() +
    theme(legend.position = "bottom",
          legend.title = element_blank())



  pdf(sprintf("%s/extended_error_barplot.(%s.vs.%s).%s.%s%s.pdf", out_path,
              group1,
              group2,
              project,
              ifelse(is.null(name), "", paste(name, ".", sep = "")),
              format(Sys.Date(), "%y%m%d")), height = height, width = width)

  print(p4)
  dev.off()
  #return(merged_data.sig)
}


