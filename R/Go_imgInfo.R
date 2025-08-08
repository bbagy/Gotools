
#' Retrieve and Organize Image Information
#'
#' This function organizes paths to various graphical outputs such as rarefaction curves,
#' bar charts, heatmaps, etc., and returns a list of these paths. The function checks if any
#' paths are provided and prints details about what each parameter expects if none are given.
#'
#' @param Rarefaction Path or vector of paths for rarefaction curve images.
#' @param Barchart Path or vector of paths for bar chart images.
#' @param Bac.heatmap Path or vector of paths for bacterial heatmap images.
#' @param Adivplot Path or vector of paths for alpha diversity plots.
#' @param Foreplot Path or vector of paths for forest plots.
#' @param Bdivplot Path or vector of paths for beta diversity plots.
#' @param DAplot Path or vector of paths for differential abundance plots.
#'
#' @return A list containing paths or vectors of paths for each type of graphical output.
#'
#' @examples
#' \dontrun{
#'   imgInfo <- Go_imgInfo(
#'     Rarefaction = "path/to/rarefaction.png",
#'     Barchart = c("path/to/bar1.png", "path/to/bar2.png"),
#'     Bac.heatmap = "path/to/heatmap.png"
#'   )
#'   print(imgInfo)
#' }
#'
#' @export
Go_imgInfo <- function(Overview=NA,
                       Rarefaction=NA,
                       Barchart=NA,
                       Bac.heatmap=NA,
                       RNAseq.heatmap=NA,
                       HumannHeatmap=NA,
                       Adivplot=NA,
                       Foreplot=NA,
                       Bdivplot=NA,
                       DAplot=NA,
                       EBplot=NA,
                       Network=NA,
                       SHeatmap=NA,
                       MAGs_tree=NA,
                       MAGs_upset=NA,
                       MAGs_heatmap=NA) {
  # Check if all arguments are missing and print options if they are
  if (all(is.na(c(Overview, Rarefaction, Barchart, Bac.heatmap, HumannHeatmap,
                  RNAseq.heatmap, Adivplot, Foreplot, Bdivplot, DAplot, EBplot, Network, SHeatmap, MAGs_tree, MAGs_upset, MAGs_heatmap)))) {
    cat(
      "Provide paths for the respective graphical outputs. Each parameter expects a path or a vector of paths for its respective images. \n",
      "Available parameters are: Overview, Rarefaction, Barchart, Bac.heatmap, RNAseq.heatmap, HumannHeatmap, Adivplot, Foreplot, Bdivplot, DAplot, EBplot, and Network. \n",
      "Specify paths to utilize the function fully, or review the @examples in the documentation for more details.\n"
    )
    return(invisible())
  }

  # Return a list with all the inputs, treating them directly as paths or vectors of paths
  return(list(
    overview = Overview,
    rarefaction = Rarefaction,
    barchart = Barchart,
    bac.heatmap = Bac.heatmap,
    rnaseq.heatmap = RNAseq.heatmap,
    humannheatmap = HumannHeatmap,
    adivplot = Adivplot,
    foreplot = Foreplot,
    bdivplot = Bdivplot,
    daplot = DAplot,
    ebplot = EBplot,
    ntplot = Network,
    shplot = SHeatmap,
    magstree = MAGs_tree,
    magsupset = MAGs_upset,
    magsheatmap = MAGs_heatmap
  ))
}
