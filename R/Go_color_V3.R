
#' Generate Color Palette for Taxonomic Data
#'
#' @param cdf Data frame containing taxonomic information (phylum level) and corresponding colors.
#' @param taxaName Vector of taxonomic names for which colors are assigned.
#'
#' @details
#' This function creates a color palette for taxonomic data, primarily at the phylum level. It assigns distinct colors to
#' each phylum and ensures visual differentiation between various taxonomic groups.
#'
#' @return
#' A list containing two elements:
#'   - `$color_table`: A data frame mapping each taxa to its corresponding HSV color code.
#'   - `$coloring`: A named vector where names are taxa and values are color codes.
#'
#' @examples
#' color_data <- Go_color(cdf = taxonomic_dataframe, taxaName = vector_of_taxa_names)
#' # Access the color table
#' color_table <- color_data$color_table
#' # Access the coloring vector
#' coloring_vector <- color_data$coloring
#'
#' @export

# https://alloyui.com/examples/color-picker/hsv.html


Go_color <- function(cdf, taxaName){
  cPalette <-cdf
  num_of_final_phyla = length(unique(cdf$PhylumCol)); num_of_final_phyla
  cPalette$h = 0
  cPalette$s = 0
  cPalette$v = 0

  col1 <- "Actinobacteriota" # Actinobacteriota/ Actinobacteria
  col2 <- "Bacteroidota" # Bacteroidota/ Bacteroidetes
  col3 <- "Firmicutes" # Firmicutes/ Firmicutes
  col4 <- "Fusobacteriota"# Fusobacteriota/ Fusobacteria
  col5 <- "Proteobacteria" # Proteobacteria/Proteobacteria
  col6 <- "Verrucomicrobia" # Verrucomicrobia
  col7 <- "TM7" # Verrucomicrobia
  col8 <- ""   # Patescibacteria/

  num_col1 = length(grep(col1, cPalette$PhylumCol));num_col1
  num_col2 = length(grep(col2, cPalette$PhylumCol));num_col2
  num_col3 = length(grep(col3, cPalette$PhylumCol));num_col3
  num_col4 = length(grep(col4, cPalette$PhylumCol));num_col4
  num_col5 = length(grep(col5, cPalette$PhylumCol));num_col5
  num_col6 = length(grep(col6, cPalette$PhylumCol));num_col6
  num_col7 = length(grep(col7, cPalette$PhylumCol));num_col7

  
  # Synergistetes
  number_of_other_phyla = num_of_final_phyla - ((num_col1 > 0) + (num_col2 > 0) + (num_col3 > 0) +(num_col4 > 0)+  (num_col5 > 0) + (num_col6 > 0)+ (num_col7 > 0))
  
  #print(number_of_other_phyla)
  # col1 = green_pallete
  cPalette[grep(col1, cPalette$PhylumCol), -1] = expand.grid(h=0.4, s=seq(0.3,1,length.out=num_col1), v=0.9)
  
  # col2 = purple_pallete
  cPalette[grep(col2, cPalette$PhylumCol), -1] = expand.grid(h=0.8, s=seq(0.3,1,length.out=num_col2), v=0.9) 

  # col3 = blue_pallete
  cPalette[grep(col3, cPalette$PhylumCol), -1] = expand.grid(h=0.6, s=seq(0.3,1,length.out=num_col3), v=0.9)

  # col4 = orange_pallete
  cPalette[grep(col4, cPalette$PhylumCol), -1] = expand.grid(h=0.2, s=seq(0.3,1,length.out=num_col4), v=0.9)
  
  # col5 = red_pallete
  cPalette[grep(col5, cPalette$PhylumCol), -1] = expand.grid(h=0, s=seq(0.3,1,length.out=num_col5), v=0.9)
  
  # col6 = brown_pallete
  cPalette[grep(col6, cPalette$PhylumCol), -1] = expand.grid(h=0.1, s=seq(0.3,1,length.out=num_col6), v=1)
  
  # col7 = yellow_pallete
  cPalette[grep(col7, cPalette$PhylumCol), -1] = expand.grid(h=0.1, s=seq(0.3,1,length.out=num_col7), v=1)
  
  #print(cPalette)
  #print(number_of_other_phyla)
  
  
  print(cPalette)
  
  # Add other and species name
  cPalette$PhylumCol <- taxaName
  other<-data.frame("[1_#Other]",0,0,0.75)
  names(other)<-c("PhylumCol", "h","s","v")
  color_table <- rbind(other, cPalette)
  
  
  ## hsv to color code ##
  ## taxa위치 변경을 용이 하게 하려면 hsv to color code를 해야 한다. 
  taxa_vs_color = cbind(color_table, apply(color_table[,-1], 1, function(x) hsv(x[1],x[2],x[3])))[,c(1,5)];taxa_vs_color
  colnames(taxa_vs_color) <- c("Taxa", "Color");taxa_vs_color
  class(taxa_vs_color)
  coloring <- as.character(taxa_vs_color$Color) 
  names(coloring) <- taxa_vs_color$Taxa;coloring
  
  
  functionReturningTwoValues <- function() { 
    results <- list()
    results$color_table <- color_table
    results$coloring <-coloring
    return(results) 
  }
  cat("\n")
  print("$color_table and $coloring are returned.")
  
  functionReturningTwoValues()
}
