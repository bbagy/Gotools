#' Generate a Pie Plot PDF for a Given Project
#'
#' This function generates a pie plot based on the provided data frame and plotting parameters for a specified project.
#' It allows for the creation of pie plots with up to three hierarchical levels and saves the plot in a structured project directory.
#'
#' @param df A data frame containing the variables to be plotted.
#' @param project A string specifying the project name, used for directory structuring and file naming.
#' @param pie1 The name of the column in `df` to be used for the primary pie slices.
#' @param pie2 The name of the column in `df` to be used for the secondary pie slices (optional).
#' @param pie3 The name of the column in `df` containing the percentage values for the slices (optional).
#' @param mycols A vector of colors for the slices; if NULL, default colors are used.
#' @param name A string for the specific plot name; if NULL, a default name is used.
#' @param height The height of the plot in inches, used for the saved PDF file.
#' @param width The width of the plot in inches, used for the saved PDF file.
#'
#' @details The function creates a directory for the project if it does not exist, and also a subdirectory for the PDF files.
#' It checks for user-provided colors (`mycols`) and handles cases where they are not properly defined.
#' The pie chart labels are determined based on the availability of `pie3` percentages.
#' The function prints the plot to the current graphics device and saves a PDF file to the project directory.
#'
#' @return No return value; the function is called for its side effects of creating a plot and saving a PDF.
#' @export
#' @examples
#' Go_piePlot(df = df, project = "MyProject", pie1 = "Category", pie3 = "Percentage")

Go_piePlot <- function(df,
                       project,
                       pie1,
                       pie2 = NULL,
                       pie3 = NULL,
                       orders,
                       mycols = NULL,
                       name = NULL,
                       height,
                       width) {

  # Close any open graphic devices to start fresh
  if(!is.null(dev.list())) dev.off()

  # Setup output directory
  out <- file.path(sprintf("%s_%s",project, format(Sys.Date(), "%y%m%d")))
  if(!file_test("-d", out)) dir.create(out)
  out_path <- file.path(sprintf("%s_%s/pdf",project, format(Sys.Date(), "%y%m%d")))
  if(!file_test("-d", out_path)) dir.create(out_path)

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

  print(0)
  # Initialize an empty list for storing data frames
  dfs <- list()

  # Define a vector of pie parameters
  pies <- Filter(Negate(is.null), list(pie1, pie2, pie3))


  # Loop through the pie parameters to create and store data frames
  for(pie in pies) {
    df_temp <- as.data.frame(table(df[[pie]]))
    names(df_temp) <- c("Category", "Count") # Standardize column names
    df_temp$pie.group <- pie
    dfs[[length(dfs) + 1]] <- df_temp # Append to the list
  }


  # Combine all data frames in the list
  combined <- do.call(rbind, dfs)


  if(!"pie.group" %in% colnames(combined)) {
    stop("Column 'pie.group' is not found in the combined dataframe.")
  }


  # Check orders
  if (!is.null(orders)){
    combined$Category <- factor(combined$Category, levels = orders)
  }else{
    combined$Category <- factor(combined$Category)
  }

  combined$pie.group <- factor(combined$pie.group, levels = pies[!sapply(pies, is.null)])


  print(1)
  combined_aggregated <- combined %>%
    group_by(pie.group, Category) %>% # Grouping by both type and Category might be needed based on your hierarchy.
    summarise(val = sum(Count), .groups = 'drop') # Explicitly drop grouping


  print(2)
  # Calculate the total values for each type and the percentage
  if(!"pie.group" %in% colnames(combined_aggregated)) {
    print(colnames(combined_aggregated))
    stop("Column 'pie.group' is not found in 'combined_aggregated'.")
  }

  combined_aggregated <- combined_aggregated %>%
    group_by(pie.group) %>%
    mutate(Total = sum(val), # Calculate total for each type
           Percentage = (val / Total) * 100)


  combined_aggregated$Label <- paste0(combined_aggregated$Category, " (", round(combined_aggregated$Percentage, 1), "%)")

  print(3)

  p <- ggplot(combined_aggregated, aes(x = pie.group, y = val, fill = Category)) +
    geom_bar(stat = "identity", position = "fill") +   theme_minimal() +
    geom_text(aes(label = Label), position = position_fill(vjust = 0.5), size = 3, color = "black")



  if(!is.null(mycols)){
    p <- p + scale_fill_manual(values = mycols)
  }else{
    p <- p
  }

  # Add title and subtitles if they are not NULL
  p <- p + labs(title = sprintf("Distribution of %s%s%s", pie1,
                                ifelse(is.null(pie2), "", paste("-",pie2, sep = "")),
                                ifelse(is.null(pie3), "", paste("-", pie3, sep = ""))))





  p1 <- p + coord_polar(theta = "y")
  print(p1)
  # Save the plot as a PDF
  pdf(sprintf("%s/pie.%s.%s%s.pdf", out_path,
              project,
              ifelse(is.null(name), "", paste(name, ".", sep = "")),
              format(Sys.Date(), "%y%m%d")), height = height, width = width)

  print(p)
  print(p1)

  dev.off()
}

