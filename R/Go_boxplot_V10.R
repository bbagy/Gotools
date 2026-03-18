
#' Generate Boxplots for Multiple Variables
#'
#' @param df Data frame containing the data to be plotted. Provide either \code{df} or \code{psIN}.
#' @param psIN Phyloseq object. If \code{outcomes} are alpha diversity metrics (e.g. "Shannon"),
#'   \code{Go_adiv()} is called internally. If taxa names are given, rank is auto-detected from
#'   \code{tax_table}.
#' @param cate.vars Categorical variables to be used for the boxplot x-axis.
#' @param project Project name used for output file naming.
#' @param outcomes Numeric variables (or alpha metrics / taxa names when using psIN) to plot on y-axis.
#' @param orders Vector of ordered factor levels for the categorical variables.
#' @param mycols Custom color palette for the plots.
#' @param combination Number of group combinations to display in the boxplot.
#' @param ylim Y-axis limits for the plot.
#' @param title Title of the plot.
#' @param facet Optional variable for creating facetted plots.
#' @param paired Indicates if the data points are paired.
#' @param name Optional name for saving the plot.
#' @param ncol Number of columns for facet wrapping in the plot.
#' @param addnumber Boolean to add the number of samples in each group. Default TRUE.
#' @param statistics Whether to perform statistical tests. Default TRUE.
#' @param parametric Whether to use parametric tests. Default FALSE.
#' @param model Statistical engine: "nonparametric", "parametric", or "lmm". If NULL, inferred from \code{parametric}.
#' @param covariates Optional covariate column names for adjusted models (ANCOVA/LMM).
#' @param p_adjust P-value adjustment method for pairwise tests (e.g., "BH", "bonferroni"). Default "BH".
#' @param xangle Angle of x-axis labels. Default 90.
#' @param cutoff Significance level for statistical tests. Default 0.1.
#' @param min_height Minimum plot panel height in inches (excluding x-axis labels). Default 3.
#' @param min_width Minimum plot width per column in inches. Default 3.
#' @param plotCols Number of columns in the multiplot layout. Default 1.
#' @param plotRows Number of rows in the multiplot layout. Default 1.
#'
#' @details
#' PDF dimensions are calculated automatically from the data:
#' \itemize{
#'   \item \strong{Width}: \code{max(min_width, n_groups * 0.5) * plotCols}
#'   \item \strong{Panel height}: \code{max(min_height, 2 + n_comparisons * 0.3) * plotRows}
#'   \item \strong{Label height}: \code{max(0.5, max_label_chars * 0.09)} added for rotated x-axis labels
#' }
#' Panel height and label height are allocated separately via \code{cowplot::plot_grid()},
#' so long group names do not eat into the boxplot panel space.
#'
#' @return A PDF file containing the boxplots.
#'
#' @examples
#' Go_boxplot(df = my_data, cate.vars = "Group", project = "MyProject",
#'            outcomes = c("Variable1", "Variable2"))
#'
#' Go_boxplot(psIN = ps, cate.vars = "Group", project = "MyProject",
#'            outcomes = c("Shannon", "Chao1"))
#'
#' @export

Go_boxplot <- function(df          = NULL,
                       psIN        = NULL,
                       cate.vars,
                       project,
                       outcomes,
                       orders      = NULL,
                       mycols      = NULL,
                       combination = NULL,
                       ylim        = NULL,
                       title       = NULL,
                       facet       = NULL,
                       paired      = NULL,
                       name        = NULL,
                       ncol        = NULL,
                       addnumber   = TRUE,
                       statistics  = TRUE,
                       parametric  = FALSE,
                       model       = NULL,
                       covariates  = NULL,
                       p_adjust    = "BH",
                       xangle      = 90,
                       cutoff      = 0.1,
                       min_height  = 3,
                       min_width   = 3,
                       plotCols    = 1,
                       plotRows    = 1) {

  # ── inner helpers ───────────────────────────────────────────────────────────
  has_facet <- function(x) !is.null(x) && length(x) >= 1

  coerce_to_ordered_factor <- function(x, orders, var_name = "variable") {
    x_chr <- as.character(x)
    if (is.null(orders) || length(orders) == 0) return(factor(x_chr))
    ord_chr <- as.character(orders)
    bad <- setdiff(unique(x_chr), ord_chr)
    if (length(bad) > 0)
      stop(sprintf("`%s` has levels not in `orders`: %s", var_name, paste(bad, collapse = ", ")))
    factor(x_chr, levels = intersect(ord_chr, x_chr))
  }

  build_comparisons <- function(cbn_mat) {
    lapply(seq_len(ncol(cbn_mat)), function(i) cbn_mat[, i])
  }

  build_plot_title <- function(base_title, test_name, pval) {
    sprintf("%s%s%s%s",
            base_title,
            ifelse(is.null(test_name), "", paste0("\n", test_name, " ")),
            ifelse(is.null(pval),      "", "p= "),
            ifelse(is.null(pval),      "", paste0(pval, " ")))
  }

  build_plot_subtitle <- function(stat_res, p_adjust, use_covariates, covariates_label) {
    parts <- character(0)
    method_label <- if (!is.null(stat_res$test.name) &&
                        stat_res$test.name %in% c("ANCOVA", "LMM")) stat_res$test.name else NULL
    if (!is.null(method_label))
      parts <- c(parts, paste0("method=", method_label))
    if (!is.null(stat_res$annotation) && !is.null(method_label))
      parts <- c(parts, sprintf("pairwise (adjust=%s)", p_adjust))
    if (isTRUE(use_covariates))
      parts <- c(parts, sprintf("covariates=%s", covariates_label))
    if (length(parts) == 0) return(NULL)
    paste(parts, collapse = "\n")
  }

  sanitize_tag <- function(x) gsub("[^A-Za-z0-9._-]+", "-", x)

  make_base_plot <- function(dat, mvar, oc) {
    ggplot(dat, aes(x = !!sym(mvar), y = !!sym(oc))) +
      labs(y = oc, x = NULL) +
      theme(strip.background = element_blank(),
            text              = element_text(size = 8),
            axis.text.x       = element_text(angle = xangle, hjust = 1, vjust = 0.5, size = 8),
            plot.title        = element_text(size = 8),
            plot.subtitle     = element_text(size = 6, lineheight = 0.9))
  }

  add_geom_layers <- function(p1, dat, mvar, paired, mycols, dot.size, box.tickness) {
    if (!is.null(paired)) {
      if (!is.null(mycols)) p1 <- p1 + scale_color_manual(values = mycols)
      p1 +
        geom_boxplot(aes(colour = !!sym(mvar)), outlier.shape = NA,
                     linewidth = box.tickness, show.legend = FALSE) +
        geom_point(aes(colour = !!sym(mvar), group = !!sym(paired)),
                   alpha = 0.8, size = dot.size, position = position_dodge(0.3),
                   show.legend = FALSE) +
        geom_line(aes(group = !!sym(paired)), color = "grey50",
                  linewidth = 0.3, position = position_dodge(0.3)) +
        theme(legend.title = element_blank(), legend.position = "bottom",
              legend.justification = "left",
              legend.box.margin = ggplot2::margin(0, 0, 0, -1, "cm"))
    } else if (max(table(dat[[mvar]])) > 150) {
      if (!is.null(mycols)) p1 <- p1 + scale_fill_manual(values = mycols)
      p1 + geom_boxplot(aes(fill = !!sym(mvar)), outlier.shape = NA,
                        linewidth = box.tickness, show.legend = FALSE)
    } else {
      if (!is.null(mycols)) p1 <- p1 + scale_color_manual(values = mycols)
      p1 +
        geom_boxplot(aes(colour = !!sym(mvar)), outlier.shape = NA,
                     linewidth = box.tickness, show.legend = FALSE) +
        geom_jitter(aes(colour = !!sym(mvar)), shape = 16, alpha = 0.8,
                    size = dot.size, position = position_jitter(0.2),
                    show.legend = FALSE)
    }
  }

  add_facet_and_guides <- function(p1, facet, ncol, df) {
    if (has_facet(facet)) {
      nc <- if (is.null(ncol)) length(unique(df[[facet[1]]])) else ncol
      p1 + facet_wrap(as.formula(sprintf("~ %s", paste(facet, collapse = "+"))),
                      scales = "free_x", ncol = nc) +
        guides(color = "none", size = "none", shape = "none")
    } else {
      p1 + guides(color = "none", size = "none", shape = "none")
    }
  }

  # Finalize panel style and store sizing metadata as attributes
  finalize_plot <- function(p1, n_grp, n_comp, max_lbl_chars) {
    p1 <- p1 + theme(panel.grid       = element_blank(),
                     panel.background = element_rect(fill = "white", colour = "black",
                                                     linewidth = 0.5, linetype = "solid"))
    attr(p1, "n_grp")         <- n_grp
    attr(p1, "n_comp")        <- n_comp
    attr(p1, "max_lbl_chars") <- max_lbl_chars
    p1
  }

  # Split plot into panel (fixed height) + x-axis labels (variable height)
  apply_label_split <- function(p, panel_h, label_h) {
    p_panel <- p + theme(axis.text.x  = element_blank(),
                         axis.ticks.x = element_blank(),
                         plot.margin  = ggplot2::margin(5, 5, 2, 5, "pt"))

    # Extract x-axis grob from full plot
    gt       <- ggplot_gtable(ggplot_build(p))
    axis_idx <- which(gt$layout$name == "axis-b")
    if (length(axis_idx) == 0) return(p_panel)

    axis_grob <- gt$grobs[[axis_idx]]
    p_xlab    <- cowplot::ggdraw() +
                   cowplot::draw_grob(axis_grob, x = 0, y = 0, width = 1, height = 1)

    cowplot::plot_grid(p_panel, p_xlab,
                       ncol        = 1,
                       rel_heights = c(panel_h, label_h))
  }

  # ── setup ───────────────────────────────────────────────────────────────────
  if (!is.null(dev.list())) dev.off()

  if (has_facet(facet) && isTRUE(addnumber)) {
    addnumber <- FALSE
    message("facet detected: addnumber forced to FALSE")
  }

  out      <- sprintf("%s_%s",     project, format(Sys.Date(), "%y%m%d"))
  out_path <- sprintf("%s_%s/pdf", project, format(Sys.Date(), "%y%m%d"))
  for (d in c(out, out_path)) if (!dir.exists(d)) dir.create(d, recursive = TRUE)
  set.seed(151)

  if (is.function(name)) name <- NULL
  tt <- try(mycols, TRUE); if (inherits(tt, "try-error")) mycols <- NULL
  tt <- try(orders, TRUE); if (inherits(tt, "try-error")) orders <- NULL

  # ── resolve df / psIN ────────────────────────────────────────────────────────
  if (!is.null(df) && !is.null(psIN))
    stop("Provide either `df` or `psIN`, not both.")
  if (is.null(df) && is.null(psIN))
    stop("Either `df` or `psIN` must be provided.")

  if (!is.null(psIN)) {
    ALPHA_METRICS   <- c("Observed", "Chao1", "ACE", "Shannon",
                         "Simpson", "InvSimpson", "Fisher", "PD")
    alpha_requested <- outcomes[outcomes %in% ALPHA_METRICS]

    if (length(alpha_requested) > 0) {
      message(sprintf("[Go_boxplot] alpha diversity mode: %s",
                      paste(alpha_requested, collapse = ", ")))
      df <- Go_adiv(psIN = psIN, project = project, alpha_metrics = alpha_requested)

    } else {
      tax_mat       <- as.data.frame(tax_table(psIN))
      detected_rank <- NULL
      for (rnk in colnames(tax_mat)) {
        if (any(outcomes %in% tax_mat[[rnk]])) { detected_rank <- rnk; break }
      }
      if (is.null(detected_rank))
        stop("outcomes not found in tax_table. Check taxa names or use `df=` instead.")
      message(sprintf("[Go_boxplot] taxa mode: rank = %s", detected_rank))

      ps_glom  <- tax_glom(psIN, taxrank = detected_rank)
      otu_mat  <- as.data.frame(t(otu_table(ps_glom)))
      colnames(otu_mat) <- as.data.frame(tax_table(ps_glom))[[detected_rank]]
      taxa_sel <- otu_mat[, outcomes[outcomes %in% colnames(otu_mat)], drop = FALSE]
      meta     <- as.data.frame(sample_data(psIN))
      df       <- merge(meta, taxa_sel, by = "row.names")
      rownames(df)  <- df$Row.names
      df$Row.names  <- NULL
    }
  }

  # ── resolve model ────────────────────────────────────────────────────────────
  has_covariates <- !is.null(covariates) && length(covariates) > 0
  model_is_auto  <- is.null(model)
  resolved_model <- if (model_is_auto) {
    if (!is.null(paired))        "lmm"
    else if (has_covariates)     "parametric"
    else if (isTRUE(parametric)) "parametric"
    else                         "nonparametric"
  } else {
    tolower(as.character(model))
  }
  if (!resolved_model %in% c("nonparametric", "parametric", "lmm"))
    stop("`model` must be one of: 'nonparametric', 'parametric', 'lmm'.")

  covariates_label <- if (has_covariates) paste(covariates, collapse = "+") else NULL
  use_covariates   <- !is.null(covariates_label) && resolved_model %in% c("parametric", "lmm")

  if (!is.null(covariates_label) && resolved_model == "nonparametric" && !model_is_auto)
    warning("`covariates` are ignored when model='nonparametric'.")
  if (resolved_model == "lmm" && is.null(paired))
    warning("model='lmm' requested but `paired` is NULL. LMM cannot be fitted.")

  method_file_tag     <- if (resolved_model == "lmm" && !is.null(paired)) "(method=LMM)."
                         else if (isTRUE(use_covariates) && resolved_model == "parametric") "(method=ANCOVA)."
                         else ""
  covariates_file_tag <- if (isTRUE(use_covariates))
                           paste0("(cov=", sanitize_tag(covariates_label), ").")
                         else ""

  dot.size     <- if (min_height * min_width <= 6) 0.7 else if (min_height * min_width < 10) 1.0 else 1.5
  box.tickness <- if (min_height * min_width <= 6) 0.3 else if (min_height * min_width < 10) 0.4 else 0.5

  file_name <- paste0(
    "box.", project, ".",
    ifelse(is.null(facet),       "", paste0(facet,      ".")),
    ifelse(is.null(paired),      "", paste0("(paired=", paired,      ").")),
    ifelse(is.null(combination), "", paste0("(cbn=",    combination, ").")),
    method_file_tag, covariates_file_tag,
    ifelse(is.null(name), "", paste0(name, ".")),
    format(Sys.Date(), "%y%m%d"), ".pdf"
  )

  # ── pass 1: build plotlist ───────────────────────────────────────────────────
  plotlist <- list()

  for (mvar in cate.vars) {
    if (length(unique(df[[mvar]])) < 2)          next
    if (has_facet(facet) && mvar %in% facet)     next

    df[[mvar]] <- as.character(df[[mvar]])
    df[[mvar]][df[[mvar]] == ""] <- NA
    df.na      <- df[!is.na(df[[mvar]]), ]
    df.na[[mvar]] <- coerce_to_ordered_factor(df.na[[mvar]], orders, var_name = mvar)

    if (isTRUE(addnumber)) {
      for (grp_name in levels(df.na[[mvar]])) {
        new_n <- paste0(grp_name, " (n=", sum(df.na[[mvar]] == grp_name), ")")
        levels(df.na[[mvar]])[levels(df.na[[mvar]]) == grp_name] <- new_n
      }
    }

    message(sprintf("## %s (without NA: %d/%d) ##", mvar, nrow(df.na), nrow(df)))
    if (length(unique(df.na[[mvar]])) == 1) next

    run_stats <- function(dat, comparisons) {
      if (!statistics)
        return(list(test.name = NULL, pval = NULL, testmethod = NULL, annotation = NULL))
      Go_boxplot_stats_engine(
        df = dat, mvar = mvar, oc = oc,
        comparisons = comparisons, model = resolved_model,
        parametric = parametric, covariates = covariates,
        paired = paired, facet = facet, p_adjust = p_adjust
      )
    }

    assemble_plot <- function(dat, my_comparisons) {
      stat_res      <- run_stats(dat, my_comparisons)
      title_text    <- build_plot_title(ifelse(is.null(title), mvar, title),
                                        stat_res$test.name, stat_res$pval)
      subtitle_text <- if (statistics)
        build_plot_subtitle(stat_res, p_adjust, use_covariates, covariates_label) else NULL

      p1 <- make_base_plot(dat, mvar, oc)
      p1 <- p1 + labs(title = title_text, subtitle = subtitle_text)
      p1 <- add_geom_layers(p1, dat, mvar, paired, mycols, dot.size, box.tickness)

      if (statistics)
        p1 <- Go_boxplot_add_stats_layer(p1 = p1, stat_res = stat_res,
                                         my_comparisons = my_comparisons,
                                         paired = paired, cutoff = cutoff)

      if (!is.null(ylim) && oc != "Chao1") p1 <- p1 + ylim(ylim[1], ylim[2])

      p1 <- add_facet_and_guides(p1, facet, ncol, df)

      # metadata for auto-sizing
      n_grp         <- length(unique(dat[[mvar]]))
      n_comp        <- length(my_comparisons)
      max_lbl_chars <- max(nchar(as.character(levels(dat[[mvar]]))), na.rm = TRUE)

      finalize_plot(p1, n_grp, n_comp, max_lbl_chars)
    }

    # ── combination branch ───────────────────────────────────────────────────
    if (!is.null(combination)) {
      message(sprintf("Combination n=%d", combination))
      if (combination < 2 || combination > 10) {
        message("combination should be 2~10 only."); next
      }
      group_comparisons <- build_comparisons(combn(x = levels(df.na[[mvar]]), m = combination))

      for (i in seq_along(group_comparisons)) {
        group.combination <- unlist(group_comparisons[[i]])
        df.cbn  <- df.na[df.na[[mvar]] %in% group.combination[seq_len(combination)], ]
        df.cbn[[mvar]] <- factor(df.cbn[[mvar]])
        cbn            <- combn(x = levels(df.cbn[[mvar]]), m = 2)
        my_comparisons <- build_comparisons(cbn)
        if (combination != 2) my_comparisons <- my_comparisons[seq_len(combination - 1)]

        if (has_facet(facet)) {
          for (fc in facet) {
            df.cbn[[fc]] <- as.character(df.cbn[[fc]])
            df.cbn[[fc]][df.cbn[[fc]] == ""] <- NA
            df.cbn <- df.cbn[!is.na(df.cbn[[fc]]), ]
            df.cbn[[fc]] <- coerce_to_ordered_factor(df.cbn[[fc]], orders, var_name = fc)
          }
        }
        for (oc in outcomes) {
          message(oc)
          plotlist[[length(plotlist) + 1]] <- assemble_plot(df.cbn, my_comparisons)
        }
      }

    } else {
      # ── no combination branch ──────────────────────────────────────────────
      my_comparisons <- build_comparisons(combn(x = levels(df.na[[mvar]]), m = 2))

      if (has_facet(facet)) {
        for (fc in facet) {
          df.na[[fc]] <- as.character(df.na[[fc]])
          df.na[[fc]][df.na[[fc]] == ""] <- NA
          df.na <- df.na[!is.na(df.na[[fc]]), ]
          df.na[[fc]] <- coerce_to_ordered_factor(df.na[[fc]], orders, var_name = fc)
        }
      }
      for (oc in outcomes) {
        message(oc)
        plotlist[[length(plotlist) + 1]] <- assemble_plot(df.na, my_comparisons)
      }
    }
  }

  if (length(plotlist) == 0) {
    message("[Go_boxplot] No plots to render.")
    return(invisible(NULL))
  }

  # ── pass 2: calculate PDF dimensions ────────────────────────────────────────
  all_n_grp   <- sapply(plotlist, function(p) attr(p, "n_grp")         %||% 3)
  all_n_comp  <- sapply(plotlist, function(p) attr(p, "n_comp")        %||% 0)
  all_max_lbl <- sapply(plotlist, function(p) attr(p, "max_lbl_chars") %||% 10)

  max_n_grp  <- max(all_n_grp,   na.rm = TRUE)
  max_n_comp <- max(all_n_comp,  na.rm = TRUE)
  max_lbl    <- max(all_max_lbl, na.rm = TRUE)

  panel_h <- max(min_height, 2 + max_n_comp * 0.3)   # panel + stat brackets
  label_h <- max(0.5, max_lbl * 0.09)                 # rotated x-axis labels
  pdf_h   <- (panel_h + label_h) * plotRows
  pdf_w   <- max(min_width, max_n_grp * 0.5) * plotCols

  message(sprintf("[Go_boxplot] PDF: %.1f x %.1f in  (panel=%.1f, labels=%.1f, w/grp=%.2f)",
                  pdf_w, pdf_h, panel_h, label_h, max_n_grp * 0.5))

  # Apply panel / label split to each plot
  plotlist <- lapply(plotlist, function(p) {
    ph <- max(min_height, 2 + (attr(p, "n_comp") %||% 0) * 0.3)
    lh <- max(0.5, (attr(p, "max_lbl_chars") %||% 10) * 0.09)
    apply_label_split(p, ph, lh)
  })

  # ── render ───────────────────────────────────────────────────────────────────
  pdf(file.path(out_path, file_name), height = pdf_h, width = pdf_w)

  multiplot <- function(..., plotlist = NULL, cols = 1, rows = 1) {
    plots    <- c(list(...), plotlist)
    numPlots <- length(plots)
    i <- 1
    while (i <= numPlots) {
      numToPlot <- min(numPlots - i + 1, cols * rows)
      layout    <- matrix(seq(i, i + cols * rows - 1), ncol = cols, nrow = rows, byrow = TRUE)
      if (numToPlot == 1) {
        print(plots[[i]])
      } else {
        grid::grid.newpage()
        grid::pushViewport(grid::viewport(layout = grid::grid.layout(nrow(layout), ncol(layout))))
        for (j in seq(i, i + numToPlot - 1)) {
          idx <- as.data.frame(which(layout == j, arr.ind = TRUE))
          print(plots[[j]], vp = grid::viewport(layout.pos.row = idx$row,
                                                layout.pos.col = idx$col))
        }
      }
      i <- i + numToPlot
    }
  }

  multiplot(plotlist = plotlist, cols = plotCols, rows = plotRows)
  dev.off()
}

# null-coalescing operator (local to this file's scope via environment)
`%||%` <- function(a, b) if (!is.null(a)) a else b
