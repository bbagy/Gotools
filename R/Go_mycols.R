#' Generate custom color palettes
#'
#' Returns a color vector from a built-in `Gotools` palette, an
#' `RColorBrewer` qualitative palette, or a `yarrr::piratepal()` palette.
#' Previewing the palette is optional.
#'
#' @param palette Optional palette name. This is the recommended modern entry
#'   point and accepts names from the built-in sets, `RColorBrewer`, or
#'   `yarrr::piratepal()`.
#' @param preview Logical; if `TRUE`, draw a simple palette preview. Default
#'   `TRUE`.
#' @param ... Legacy palette selectors kept only for backward compatibility.
#'   New code should use `palette` only.
#'
#' @return A character vector of hex colors.
#'
#' @examples
#' Go_myCols(palette = "cols1")
#' Go_myCols(palette = "Set1")
#' Go_myCols(palette = "basel")
#' Go_myCols(custom_cols = "cols2", preview = FALSE)
#'
#' @export
#' @importFrom grDevices dev.list dev.off
#' @importFrom graphics barplot
#' @importFrom RColorBrewer brewer.pal brewer.pal.info display.brewer.all
#' @importFrom yarrr piratepal
Go_myCols <- function(
    palette = NULL,
    preview = TRUE,
    ...
) {
  builtin_palettes <- list(
    cols1 = c("#1170aa", "#fc7d0b", "#76B7B2", "#B07AA1", "#E15759",
              "#59A14F", "#EDC948", "#FF9DA7", "#9C755F", "#BAB0AC"),
    cols2 = c("#CBD588", "#5F7FC7", "orange", "#DA5724", "#508578",
              "#CD9BCD", "#AD6F3B", "#673770", "#D14285", "#652926",
              "#C84248", "#8569D5", "#5E738F", "#D1A33D", "#8A7C64",
              "#599861"),
    cols3 = c("#323337FF", "#FF0000", "#009E73", "#0000FF", "#FFFF00",
              "#00FFFF", "#FF00FF", "#C0C0C0", "#666666", "#D86C4FFF",
              "#66A61E", "#000080", "#808000", "#800080", "#F57C00",
              "#FFC0CB", "#DC143C", "#ADFF2F", "#00FF7F", "#00CED1",
              "#795548", "#673AB7", "#1E90FF", "#F1C40F", "#FF69B4",
              "#FF6347", "#4169E1", "#FFA07A", "#C0392B", "#95CE8AFF",
              "#20B2AA"),
    cols4 = c("#2166AC", "#74ADD1", "#990000", "#D73027", "#F46D43",
              "#FDAE61", "#66BD63", "#A6DBA0", "#01665E", "#1B9E77"),
    taxaStack = c(
      "gray", "#1170AA", "#FC7D0B", "#E15759", "#59A14F", "#EDC948", "#B07AA1",
      "#9C755F", "#5F7FC7", "#DA5724", "#673770", "#D14285", "#C84248",
      "#8569D5", "#D1A33D", "#323337", "#FF0000", "#009E73", "#0000FF",
      "#D86C4F", "#66A61E", "#000080", "#800080", "#F57C00", "#DC143C",
      "#00CED1", "#1E90FF", "#F1C40F", "#C0392B", "#20B2AA", "#2166AC",
      "#74ADD1", "#990000", "#D73027", "#F46D43", "#FDAE61", "#66BD63",
      "#A6DBA0", "#01665E", "#1B9E77", "#FF69B4", "#795548", "#4169E1",
      "#AD6F3B", "#508578", "#CD9BCD", "#8A7C64", "#95CE8A", "#76B7B2"
    )
  )

  resolve_palette_name <- function() {
    legacy <- list(...)
    allowed_legacy <- c("custom_cols", "custumCols", "RColorBrewer", "piratepal")
    unknown_legacy <- setdiff(names(legacy), allowed_legacy)
    if (length(unknown_legacy) > 0) {
      stop(
        "[Go_myCols] Unknown argument(s): ",
        paste(unknown_legacy, collapse = ", "),
        ".\nUse `palette = ...` in new code."
      )
    }

    vals <- c(
      palette = palette,
      custom_cols = legacy$custom_cols,
      custumCols = legacy$custumCols,
      RColorBrewer = legacy$RColorBrewer,
      piratepal = legacy$piratepal
    )
    vals <- vals[!vapply(vals, is.null, logical(1))]

    if (length(vals) == 0) {
      return(NULL)
    }
    if (length(vals) > 1 && length(unique(unname(vals))) > 1) {
      stop(
        "[Go_myCols] Please provide only one palette selector.\n",
        "Use one of: palette, custom_cols, custumCols, RColorBrewer, piratepal."
      )
    }
    unname(vals[[1]])
  }

  preview_palette <- function(cols, title) {
    if (!isTRUE(preview)) {
      return(invisible(NULL))
    }
    if (!is.null(grDevices::dev.list())) {
      grDevices::dev.off()
    }
    graphics::barplot(rep(1, length(cols)), col = cols, main = title, yaxt = "n", border = NA)
    invisible(NULL)
  }

  selected <- resolve_palette_name()

  if (is.null(selected)) {
    available_brewer <- rownames(RColorBrewer::brewer.pal.info[
      RColorBrewer::brewer.pal.info$category == "qual", , drop = FALSE
    ])
    available_pirate <- c(
      "basel", "pony", "xmen", "decision", "southpark", "google", "eternal",
      "evildead", "usualsuspects", "ohbrother", "appletv", "brave", "bugs",
      "cars", "nemo", "rat", "up", "espresso", "ipod", "info", "info2"
    )
    msg <- paste0(
      "[Go_myCols] Please choose one palette.\n",
      "  Built-in: ", paste(names(builtin_palettes), collapse = ", "), "\n",
      "  RColorBrewer: ", paste(available_brewer, collapse = ", "), "\n",
      "  piratepal: ", paste(available_pirate, collapse = ", ")
    )
    if (interactive()) {
      message(msg)
      if (isTRUE(preview)) {
        RColorBrewer::display.brewer.all(type = "qual")
        yarrr::piratepal("all")
      }
    } else {
      stop(msg)
    }
    return(invisible(NULL))
  }

  if (selected %in% names(builtin_palettes)) {
    cols <- builtin_palettes[[selected]]
    preview_palette(cols, selected)
    return(cols)
  }

  brewer_info <- RColorBrewer::brewer.pal.info
  if (selected %in% rownames(brewer_info)) {
    max_n <- brewer_info[selected, "maxcolors"]
    cols <- RColorBrewer::brewer.pal(max_n, selected)
    preview_palette(cols, selected)
    return(cols)
  }

  pirate_cols <- tryCatch(
    yarrr::piratepal(palette = selected),
    error = function(e) NULL
  )
  if (!is.null(pirate_cols)) {
    cols <- as.character(unname(pirate_cols))
    preview_palette(cols, selected)
    return(cols)
  }

  stop(
    "[Go_myCols] Unknown palette: ", selected, "\n",
    "Use a built-in palette (", paste(names(builtin_palettes), collapse = ", "),
    "), a qualitative RColorBrewer palette, or a yarrr pirate palette."
  )
}
