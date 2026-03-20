#' Generate Ranked Lollipop Plots for DA and ConDA Results
#'
#' This function follows the same directory-first contract as `Go_volcanoPlot()`.
#' It reads CSV result files or bridge CSV files, recognizes the originating
#' DA tool, and generates ranked lollipop plots into `pdf/DA_plot` or
#' `pdf/ConDa_plot`.
#'
#' @param project Project name or identifier.
#' @param file_path Optional directory path to scan for CSV files.
#' @param files Optional filename pattern used together with `file_path`.
#' @param result Optional result directory returned by a Go DA-family tool or a
#'   ConDA bridge directory.
#' @param top_n Number of top-ranked features to display.
#' @param mycols Optional two-color vector for negative / positive direction.
#' @param name Optional plot label used in output filenames.
#' @param font Base font size.
#' @param height PDF height.
#' @param width PDF width.
#'
#' @return The function saves PDF files and returns `NULL` invisibly.
#'
#' @export
Go_lollipopPlot <- function(project,
                            file_path = NULL,
                            files = NULL,
                            result = NULL,
                            top_n = 20,
                            mycols = NULL,
                            name = NULL,
                            font = 10,
                            height = 7,
                            width = 8) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 is required.")
  }

  first_or_na <- function(x) {
    if (length(x) == 0 || all(is.na(x))) {
      return(NA_character_)
    }
    x[[1]]
  }

  build_taxonomy_plot_label <- function(df) {
    if (nrow(df) == 0) {
      return(character(0))
    }

    pick_first_valid <- function(candidates) {
      vals <- as.character(candidates)
      vals <- vals[!is.na(vals) & nzchar(trimws(vals)) & !vals %in% c("NA", "NA NA", "__", "s__", "g__", "f__", "o__", "c__", "p__", "k__")]
      if (length(vals) == 0) {
        return(NA_character_)
      }
      vals[1]
    }

    taxonomy_priority <- intersect(c("Species", "Genus", "Family", "Order", "Class", "Phylum", "TaxaName"), colnames(df))
    asv_col <- if ("ASV" %in% colnames(df)) {
      "ASV"
    } else if ("Row.names" %in% colnames(df)) {
      "Row.names"
    } else {
      NULL
    }

    labels <- vapply(seq_len(nrow(df)), function(i) {
      vals <- vapply(taxonomy_priority, function(col) df[i, col], character(1))
      lbl <- pick_first_valid(vals)
      if (is.na(lbl) && !is.null(asv_col)) {
        lbl <- as.character(df[i, asv_col])
      }
      if (is.na(lbl) || !nzchar(trimws(lbl))) {
        lbl <- paste0("Feature_", i)
      }
      lbl
    }, character(1))

    labels
  }

  build_feature_label <- function(df) {
    if (nrow(df) == 0) {
      return(character(0))
    }
    if (all(c("KOName") %in% colnames(df))) {
      vals <- as.character(df$KOName)
      vals[is.na(vals) | !nzchar(trimws(vals))] <- if ("KO" %in% colnames(df)) as.character(df$KO[is.na(vals) | !nzchar(trimws(vals))]) else NA_character_
      vals[is.na(vals) | !nzchar(trimws(vals))] <- if ("ASV" %in% colnames(df)) as.character(df$ASV[is.na(vals) | !nzchar(trimws(vals))]) else paste0("Feature_", seq_len(nrow(df)))[is.na(vals) | !nzchar(trimws(vals))]
      return(vals)
    }
    if (all(c("symbol") %in% colnames(df))) {
      vals <- as.character(df$symbol)
      vals[is.na(vals) | !nzchar(trimws(vals))] <- if ("ASV" %in% colnames(df)) as.character(df$ASV[is.na(vals) | !nzchar(trimws(vals))]) else paste0("Feature_", seq_len(nrow(df)))[is.na(vals) | !nzchar(trimws(vals))]
      return(vals)
    }
    build_taxonomy_plot_label(df)
  }

  detect_tool <- function(filename1) {
    if (grepl("condadist", filename1)) return("condadist")
    if (grepl("deseq2", filename1)) return("deseq2")
    if (grepl("aldex2", filename1)) return("aldex2")
    if (grepl("maaslin2", filename1)) return("maaslin2")
    if (grepl("maaslin", filename1)) return("maaslin2")
    if (grepl("corncob", filename1)) return("corncob")
    if (grepl("ancombc2", filename1)) return("ancom2")
    if (grepl("ancom2", filename1)) return("ancom2")
    NA_character_
  }

  first_existing <- function(df, candidates) {
    hits <- candidates[candidates %in% colnames(df)]
    if (length(hits) == 0) {
      return(NA_character_)
    }
    hits[1]
  }

  resolve_aldex_columns <- function(df) {
    est_col <- if ("diff.btw" %in% colnames(df)) {
      "diff.btw"
    } else {
      est_hits <- grep("\\.Est$", colnames(df), value = TRUE)
      est_hits <- est_hits[!grepl("kendall|spearman", est_hits, ignore.case = TRUE)]
      first_or_na(est_hits)
    }

    p_col <- if ("wi.ep" %in% colnames(df)) {
      "wi.ep"
    } else {
      p_hits <- grep("\\.pval$", colnames(df), value = TRUE)
      p_hits <- p_hits[!grepl("kendall|spearman", p_hits, ignore.case = TRUE)]
      first_or_na(p_hits)
    }

    q_col <- first_existing(df, c("wi.eBH", grep("\\.pval\\.padj$", colnames(df), value = TRUE), grep("\\.pval\\.holm$", colnames(df), value = TRUE)))
    list(effect = est_col, p = p_col, q = q_col)
  }

  build_plot_df <- function(df, tool) {
    baseline_col <- if ("baseline" %in% colnames(df)) "baseline" else if ("basline" %in% colnames(df)) "basline" else NULL
    basline <- if (!is.null(baseline_col)) unique(as.character(df[[baseline_col]]))[1] else "group1"
    smvar <- if ("smvar" %in% colnames(df)) unique(as.character(df$smvar))[1] else "group2"
    mvar <- if ("mvar" %in% colnames(df)) unique(as.character(df$mvar))[1] else "Group"

    if (tool == "deseq2") {
      effect_col <- "log2FoldChange"
      p_col <- "pvalue"
      q_col <- "padj"
      title_tool <- "deseq2"
      file_tool <- "deseq2"
      out_subdir <- "DA_plot"
      score_label <- "Signed significance score (-log10 adj. p)"
    } else if (tool == "aldex2") {
      cols <- resolve_aldex_columns(df)
      effect_col <- cols$effect
      p_col <- cols$p
      q_col <- cols$q
      title_tool <- "aldex2"
      file_tool <- "aldex2"
      out_subdir <- "DA_plot"
      score_label <- "Signed significance score (-log10 adj. p)"
    } else if (tool == "ancom2") {
      effect_col <- "lfc_ancombc"
      p_col <- "pvalue_ancombc"
      q_col <- "qvalue_ancombc"
      title_tool <- "ancom2"
      file_tool <- "ancom2"
      out_subdir <- "DA_plot"
      score_label <- "Signed significance score (-log10 adj. p)"
    } else if (tool == "corncob") {
      effect_col <- "corncob_coef"
      p_col <- "corncob_pvalue"
      q_col <- "corncob_qvalue"
      title_tool <- "corncob"
      file_tool <- "corncob"
      out_subdir <- "DA_plot"
      score_label <- "Signed significance score (-log10 adj. p)"
    } else if (tool == "maaslin2") {
      effect_col <- "maaslin2_coef"
      p_col <- "maaslin2_pvalue"
      q_col <- "maaslin2_qvalue"
      title_tool <- "maaslin2"
      file_tool <- "maaslin2"
      out_subdir <- "DA_plot"
      score_label <- "Signed significance score (-log10 adj. p)"
    } else if (tool == "condadist") {
      effect_col <- "median_effect_size"
      p_col <- "fisher_combined_p"
      q_col <- "fisher_combined_q"
      title_tool <- "ConDA-dist"
      file_tool <- "ConDA-dist"
      out_subdir <- "ConDa_plot"
      score_label <- "Signed priority score"
    } else {
      return(NULL)
    }

    if (!effect_col %in% colnames(df)) {
      return(NULL)
    }

    feature_label <- build_feature_label(df)
    effect <- suppressWarnings(as.numeric(df[[effect_col]]))
    pval <- if (p_col %in% colnames(df)) suppressWarnings(as.numeric(df[[p_col]])) else rep(NA_real_, nrow(df))
    qval <- if (!is.na(q_col) && q_col %in% colnames(df)) suppressWarnings(as.numeric(df[[q_col]])) else rep(NA_real_, nrow(df))

    base_score <- if (tool == "condadist" && "priority_score" %in% colnames(df)) {
      suppressWarnings(as.numeric(df$priority_score))
    } else {
      use_q <- is.finite(qval) & !is.na(qval)
      score_in <- ifelse(use_q, qval, pval)
      -log10(pmax(score_in, 1e-300))
    }

    signed_score <- sign(effect) * base_score
    dir_label <- ifelse(
      !is.finite(effect) | is.na(effect),
      "NS",
      ifelse(effect > 0, smvar, ifelse(effect < 0, basline, "NS"))
    )

    out <- data.frame(
      feature_label = feature_label,
      effect = effect,
      pval = pval,
      qval = qval,
      base_score = base_score,
      signed_score = signed_score,
      direction = dir_label,
      mvar = mvar,
      basline = basline,
      smvar = smvar,
      out_subdir = out_subdir,
      title_tool = title_tool,
      file_tool = file_tool,
      score_label = score_label,
      stringsAsFactors = FALSE
    )

    out <- out[is.finite(out$signed_score) & !is.na(out$signed_score), , drop = FALSE]
    if (nrow(out) == 0) {
      return(NULL)
    }
    out <- out[order(abs(out$signed_score), decreasing = TRUE), , drop = FALSE]
    out <- utils::head(out, top_n)
    out$feature_label <- factor(out$feature_label, levels = rev(out$feature_label))
    out
  }

  if (!is.null(dev.list())) dev.off()

  out_root <- file.path(sprintf("%s_%s", project, format(Sys.Date(), "%y%m%d")))
  if (!dir.exists(out_root)) dir.create(out_root)
  out_pdf <- file.path(sprintf("%s_%s/pdf", project, format(Sys.Date(), "%y%m%d")))
  if (!dir.exists(out_pdf)) dir.create(out_pdf)

  if (!is.null(result)) {
    if (!is.character(result) || !dir.exists(result)) {
      stop("result must be a directory path returned by a Go DA-family tool or ConDA bridge directory.")
    }
    tool_files <- list.files(result, pattern = "\\.csv$", full.names = FALSE)
    if (length(tool_files) == 0) stop(sprintf("No CSV files found in: %s", result))
    file_list <- data.frame(path = result, file = tool_files, stringsAsFactors = FALSE)
  } else if (!is.null(file_path) && !is.null(files)) {
    filenames <- list.files(file_path, pattern = files)
    file_list <- data.frame(path = file_path, file = filenames, stringsAsFactors = FALSE)
  } else {
    tool_dirs <- c("deseq2", "aldex2", "ancom2", "corncob", "maaslin2")
    file_list <- do.call(rbind, lapply(tool_dirs, function(tool) {
      tool_path <- sprintf("%s_%s/table/%s", project, format(Sys.Date(), "%y%m%d"), tool)
      if (!dir.exists(tool_path)) return(NULL)
      tool_files <- list.files(tool_path, pattern = "\\.csv$", full.names = FALSE)
      if (length(tool_files) == 0) return(NULL)
      data.frame(path = tool_path, file = tool_files, stringsAsFactors = FALSE)
    }))
  }

  if (is.null(file_list) || nrow(file_list) == 0) {
    message("No DA result files found.")
    return(invisible(NULL))
  }

  print(file_list)

  for (fn in seq_len(nrow(file_list))) {
    filename1 <- file.path(file_list$path[fn], file_list$file[fn])
    df <- utils::read.csv(filename1, row.names = NULL, check.names = FALSE, stringsAsFactors = FALSE)
    colnames(df) <- gsub('^"|"$', "", colnames(df))

    tool <- detect_tool(filename1)
    if (is.na(tool)) {
      next
    }
    print(tool)

    plot_df <- build_plot_df(df, tool)
    if (is.null(plot_df) || nrow(plot_df) == 0) {
      next
    }

    if (!is.null(mycols) && length(mycols) >= 2) {
      dircolors <- c(mycols[1], "grey70", mycols[2])
    } else {
      dircolors <- c("#f8766d", "grey70", "#7cae00")
    }
    bas.count <- if ("bas.count" %in% colnames(df)) unique(df$bas.count)[1] else NA
    smvar.count <- if ("smvar.count" %in% colnames(df)) unique(df$smvar.count)[1] else NA
    names(dircolors) <- c(plot_df$basline[1], "NS", plot_df$smvar[1])
    legend.labs <- c(
      if (is.na(bas.count)) plot_df$basline[1] else paste0(plot_df$basline[1], " (n=", bas.count, ")"),
      "NS",
      if (is.na(smvar.count)) plot_df$smvar[1] else paste0(plot_df$smvar[1], " (n=", smvar.count, ")")
    )
    names(legend.labs) <- c(plot_df$basline[1], "NS", plot_df$smvar[1])

    out_dir <- file.path(sprintf("%s_%s/pdf/%s", project, format(Sys.Date(), "%y%m%d"), plot_df$out_subdir[1]))
    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

    comparison_token <- NULL
    if ("comparison_token" %in% colnames(df)) {
      vals <- unique(as.character(df$comparison_token))
      vals <- vals[!is.na(vals) & vals != "" & vals != "NA"]
      if (length(vals) > 0) comparison_token <- vals[1]
    }
    if (is.null(comparison_token)) {
      name_token <- if (is.null(name)) NULL else as.character(name)
      comparison_token <- sprintf("%s.vs.%s%s",
                                  plot_df$basline[1],
                                  plot_df$smvar[1],
                                  if (is.null(name_token)) "" else paste(".", name_token, sep = ""))
    }

    p <- ggplot2::ggplot(
      plot_df,
      ggplot2::aes(x = signed_score, y = feature_label, color = direction)
    ) +
      ggplot2::geom_segment(
        ggplot2::aes(x = 0, xend = signed_score, y = feature_label, yend = feature_label),
        linewidth = 0.7,
        alpha = 0.8
      ) +
      ggplot2::geom_point(size = 3) +
      ggplot2::geom_vline(xintercept = 0, linetype = "dotted", linewidth = 0.7, color = "grey40") +
      ggplot2::scale_color_manual(values = dircolors, labels = legend.labs, drop = FALSE) +
      ggplot2::labs(
        title = sprintf("%s, %s (top=%s)", plot_df$mvar[1], plot_df$title_tool[1], top_n),
        x = plot_df$score_label[1],
        y = NULL,
        color = NULL
      ) +
      ggplot2::theme_bw() +
      ggplot2::theme(
        text = ggplot2::element_text(size = font + 6),
        plot.title = ggplot2::element_text(size = font + 8),
        axis.text.y = ggplot2::element_text(face = "italic"),
        panel.grid.major.y = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank(),
        legend.position = "bottom",
        legend.justification = c(0, 0),
        legend.box.just = "left",
        legend.key = ggplot2::element_blank(),
        aspect.ratio = 1/1.35
      )

    pdf_file <- sprintf(
      "%s/%s.lollipop.%s.(%s).%s.(top=%s).%s.pdf",
      out_dir,
      plot_df$file_tool[1],
      plot_df$mvar[1],
      comparison_token,
      project,
      top_n,
      format(Sys.Date(), "%y%m%d")
    )

    ggplot2::ggsave(filename = pdf_file, plot = p, width = width, height = height, device = "pdf")
    print(pdf_file)
  }

  invisible(NULL)
}
