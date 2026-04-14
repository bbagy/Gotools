
#' Create a Dual Y-Axis Plot
#'
#' Combines a boxplot (left Y axis) and one or two line plots of group means
#' (right Y axis) in a single panel. The secondary axis is automatically scaled
#' to match the range of the boxplot data, so both variables are visible
#' regardless of their original units.
#'
#' @param df Data frame containing sample metadata and variables to plot.
#' @param TaxTab Either a file path (CSV) to a taxa abundance table, or a data
#'   frame with samples as rows and taxa as columns. If NULL, taxa columns are
#'   assumed to already be present in \code{df}.
#' @param cate.vars Vector of categorical variables (grouping variables).
#' @param project Project name for output file naming.
#' @param Box Column name for the boxplot (left Y axis).
#' @param Line1 Column name for the first line plot (right Y axis, group means).
#' @param Line2 Optional second column name for a second line (right Y axis).
#' @param title Optional plot title. Defaults to the grouping variable name.
#' @param name Optional suffix for output file name.
#' @param mycols Optional color palette for groups.
#' @param orders Optional vector to order factor levels.
#' @param xangle Rotation angle for x-axis labels. Default 90.
#' @param addnumber Logical. Add sample counts to group labels. Default TRUE.
#' @param height Height of the output PDF.
#' @param width Width of the output PDF.
#'
#' @return Saves a PDF and returns the last plot invisibly.
#'
#' @details
#' The right Y axis is automatically scaled so that the line plot range aligns
#' with the boxplot range. This avoids the misleading \code{*10} hardcoding and
#' makes the function work correctly regardless of data magnitude.
#'
#' @examples
#' Go_dualYplot(df = meta_df, TaxTab = "taxa_table.csv",
#'              cate.vars = "Timepoint", project = "MyProject",
#'              Box = "Bacteroides", Line1 = "IL6",
#'              height = 6, width = 5)
#'
#' @param patchwork Logical. If \code{TRUE}, skip saving and return the plot object(s) for use with \code{Gg_patchwork()} or the \pkg{patchwork} package. Default \code{FALSE}.
#' @export

Go_dualYplot <- function(df, TaxTab = NULL, cate.vars, project,
                         Box, Line1, Line2 = NULL,
                         title    = NULL,
                         name     = NULL,
                         mycols   = NULL,
                         orders   = NULL,
                         xangle   = 90,
                         addnumber = TRUE,
                         height, width,
                         patchwork = FALSE) {

  if (!is.null(dev.list())) dev.off()

  # output dirs
  out      <- file.path(sprintf("%s_%s", project, format(Sys.Date(), "%y%m%d")))
  out_path <- file.path(sprintf("%s_%s/pdf", project, format(Sys.Date(), "%y%m%d")))
  for (d in c(out, out_path)) if (!dir.exists(d)) dir.create(d, recursive = TRUE)

  set.seed(151)

  tt <- try(mycols, TRUE); if (inherits(tt, "try-error")) mycols <- NULL
  tt <- try(orders, TRUE); if (inherits(tt, "try-error")) orders <- NULL
  if (is.function(name)) name <- NULL

  # ── load and merge TaxTab ────────────────────────────────────────────────
  if (!is.null(TaxTab)) {
    if (is.character(TaxTab)) {
      df2 <- read.csv(TaxTab, header = TRUE, row.names = 1, check.names = FALSE)
    } else {
      df2 <- as.data.frame(TaxTab)
    }
    # use Species column as rownames if present
    if ("Species" %in% colnames(df2)) {
      rownames(df2) <- gsub(" ", "_", df2$Species)
      df2$Species   <- NULL
    }
    df2      <- as.data.frame(t(df2))
    base_df  <- merge(df, df2, by = "row.names")
    rownames(base_df) <- base_df$Row.names
    base_df$Row.names <- NULL
  } else {
    base_df <- df
  }
  base_df$etc <- NULL

  # dot size by plot area
  dot.size     <- if (height * width <= 6) 0.7 else if (height * width < 10) 1.0 else 1.5
  box.tickness <- if (height * width <= 6) 0.3 else if (height * width < 10) 0.4 else 0.5

  plotlist_pw <- list()
  if (!isTRUE(patchwork)) {
    pdf(sprintf("%s/dualYplot.%s.%s%s%s.pdf",
                out_path, project,
                ifelse(is.null(Line1), "", paste0(Line1, ".")),
                ifelse(is.null(name),  "", paste0(name,  ".")),
                format(Sys.Date(), "%y%m%d")),
        height = height, width = width)
  }

  for (mvar in cate.vars) {

    if (length(unique(base_df[[mvar]])) < 2) next

    # NA handling
    df_w         <- base_df
    df_w[[mvar]] <- as.character(df_w[[mvar]])
    df_w[[mvar]][df_w[[mvar]] == ""] <- NA
    df_w         <- df_w[!is.na(df_w[[mvar]]), ]
    df_w[[mvar]] <- factor(df_w[[mvar]],
                           levels = if (length(orders) >= 1) orders else sort(unique(df_w[[mvar]])))

    if (length(unique(df_w[[mvar]])) < 2) next

    # add n to labels
    if (isTRUE(addnumber)) {
      lvls <- levels(df_w[[mvar]])
      new_lvls <- sapply(lvls, function(lv) {
        paste0(lv, " (n=", sum(df_w[[mvar]] == lv), ")")
      })
      levels(df_w[[mvar]]) <- new_lvls
    }

    message(sprintf("## %s (n=%d) ##", mvar, nrow(df_w)))

    # ── auto scale factor for secondary axis ──────────────────────────────
    box_vals  <- as.numeric(df_w[[Box]])
    line1_grp <- aggregate(df_w[[Line1]], list(df_w[[mvar]]), FUN = mean, na.rm = TRUE)
    colnames(line1_grp) <- c(mvar, Line1)
    line1_vals <- line1_grp[[Line1]]

    box_range  <- range(box_vals,  na.rm = TRUE)
    line1_range <- range(line1_vals, na.rm = TRUE)

    # scale: map line1 range onto box range
    if (diff(line1_range) == 0) {
      scale_fac <- 1; offset <- 0
    } else {
      scale_fac <- diff(box_range) / diff(line1_range)
      offset    <- box_range[1] - line1_range[1] * scale_fac
    }

    line1_grp[[Line1]] <- line1_grp[[Line1]] * scale_fac + offset

    # ── base boxplot ──────────────────────────────────────────────────────
    p <- ggplot() +
      theme_bw() +
      theme(strip.background = element_blank(),
            text = element_text(size = 9),
            axis.text.x = element_text(angle = xangle, hjust = 1, vjust = 0.5),
            plot.title = element_text(size = 9, face = "bold")) +
      geom_boxplot(data = df_w,
                   mapping = aes(x = !!sym(mvar), y = !!sym(Box), colour = !!sym(mvar)),
                   outlier.shape = NA, show.legend = FALSE, linewidth = box.tickness)

    if (!is.null(mycols)) p <- p + scale_color_manual(values = mycols)

    # jitter by sample size
    n_max <- max(table(df_w[[mvar]]))
    dot_alpha <- if (n_max > 500) 0.4 else 0.8
    dot_sz    <- if (n_max > 500) dot.size / 3 else if (n_max > 250) dot.size / 2 else dot.size

    p <- p + geom_jitter(data = df_w,
                         mapping = aes(x = !!sym(mvar), y = !!sym(Box), colour = !!sym(mvar)),
                         shape = 16, alpha = dot_alpha, size = dot_sz,
                         position = position_jitter(0.2), show.legend = FALSE)

    # ── line(s) ───────────────────────────────────────────────────────────
    if (!is.null(Line2)) {
      line2_grp <- aggregate(df_w[[Line2]], list(df_w[[mvar]]), FUN = mean, na.rm = TRUE)
      colnames(line2_grp) <- c(mvar, Line2)
      line2_range <- range(line2_grp[[Line2]], na.rm = TRUE)
      if (diff(line2_range) == 0) {
        scale_fac2 <- 1; offset2 <- 0
      } else {
        scale_fac2 <- diff(box_range) / diff(line2_range)
        offset2    <- box_range[1] - line2_range[1] * scale_fac2
      }
      line2_grp[[Line2]] <- line2_grp[[Line2]] * scale_fac2 + offset2

      line_combined <- merge(line1_grp, line2_grp, by = mvar)
      line_long <- tidyr::pivot_longer(line_combined, cols = c(Line1, Line2),
                                       names_to = "variable", values_to = "value")

      last_pts <- do.call(rbind, lapply(c(Line1, Line2), function(v) {
        sub <- line_long[line_long$variable == v, ]
        sub[nrow(sub), ]
      }))

      p <- p +
        geom_line(data = line_long,
                  mapping = aes(x = !!sym(mvar), y = value,
                                group = variable, color = variable),
                  linewidth = 1, inherit.aes = FALSE) +
        geom_text_repel(data = last_pts,
                        aes(x = !!sym(mvar), y = value,
                            label = gsub("_", " ", variable)),
                        size = 3, fontface = "italic") +
        guides(color = "none") +
        scale_y_continuous(
          sec.axis = sec_axis(~ (. - offset) / scale_fac, name = Line1)
        )

    } else {
      line1_long <- tidyr::pivot_longer(line1_grp, cols = Line1,
                                        names_to = "variable", values_to = "value")
      last_pt <- line1_long[nrow(line1_long), ]

      p <- p +
        geom_line(data = line1_long,
                  mapping = aes(x = !!sym(mvar), y = value,
                                group = variable, color = variable),
                  linewidth = 1, inherit.aes = FALSE) +
        geom_text_repel(data = last_pt,
                        aes(x = !!sym(mvar), y = value,
                            label = gsub("_", " ", variable)),
                        size = 3, fontface = "italic") +
        guides(color = "none") +
        scale_y_continuous(
          sec.axis = sec_axis(~ (. - offset) / scale_fac, name = Line1)
        )
    }

    p <- p + ggtitle(if (!is.null(title)) title else mvar)
    plotlist_pw[[length(plotlist_pw) + 1]] <- p
    if (!isTRUE(patchwork)) print(p)
  }

  if (isTRUE(patchwork)) return(invisible(plotlist_pw))
  dev.off()
  invisible(p)
}
