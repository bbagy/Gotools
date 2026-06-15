
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
                       MAGs_heatmap=NA,
                       Extra=NULL) {
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

  # Parse helper: supports both "token" strings and list(token=..., out_width=...)
  .parse <- function(x) {
    if (is.null(x) || (length(x) == 1 && all(is.na(x)))) return(list(tokens = x, widths = NULL))
    if (is.character(x)) return(list(tokens = x, widths = NULL))
    tokens <- vapply(x, function(item) if (is.character(item)) item else as.character(item$token), character(1))
    widths <- vapply(x, function(item) if (is.character(item) || is.null(item$out_width)) NA_character_ else as.character(item$out_width), character(1))
    list(tokens = tokens, widths = if (all(is.na(widths))) NULL else widths)
  }
  adiv_p <- .parse(Adivplot); fore_p <- .parse(Foreplot)
  bdiv_p <- .parse(Bdivplot); da_p   <- .parse(DAplot)
  bc_p   <- .parse(Barchart)

  # Normalize Extra items: derive title from token if missing
  if (!is.null(Extra)) {
    Extra <- lapply(Extra, function(item) {
      if (is.null(item$title) && !is.null(item$token)) {
        item$title <- tools::toTitleCase(
          gsub("([a-z])([A-Z])", "\\1 \\2", gsub("_", " ", item$token))
        )
      }
      item
    })
  }

  # Return a list with all the inputs, treating them directly as paths or vectors of paths
  return(list(
    overview       = Overview,
    rarefaction    = Rarefaction,
    barchart       = bc_p$tokens,   barchart_width  = bc_p$widths,
    bac.heatmap    = Bac.heatmap,
    rnaseq.heatmap = RNAseq.heatmap,
    humannheatmap  = HumannHeatmap,
    adivplot       = adiv_p$tokens, adivplot_width  = adiv_p$widths,
    foreplot       = fore_p$tokens, foreplot_width  = fore_p$widths,
    bdivplot       = bdiv_p$tokens, bdivplot_width  = bdiv_p$widths,
    daplot         = da_p$tokens,   daplot_width    = da_p$widths,
    ebplot         = EBplot,
    ntplot         = Network,
    shplot         = SHeatmap,
    magstree       = MAGs_tree,
    magsupset      = MAGs_upset,
    magsheatmap    = MAGs_heatmap,
    extra          = Extra
  ))
}
