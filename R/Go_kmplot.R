#' Go_kmplot
#'
#' Kaplan–Meier plot by quantile groups of a continuous feature, with log-rank p
#' and automatic dated output folders.
#'
#' @description
#' Given a data frame with time, event, and a continuous feature, this helper:
#' \itemize{
#'   \item derives event status (1/0) from \code{event_col} using numeric values
#'         or string matches to \code{event_positive},
#'   \item bins the feature into \code{n_group} quantile groups (labels optional),
#'   \item fits \code{survival::survfit} and computes a log-rank p-value,
#'   \item draws a KM curve via \code{survminer::ggsurvplot} and saves a PDF to
#'         \verb{<project_YYMMDD>/pdf/}.
#' }
#'
#' @param df A data.frame containing the columns named by \code{time_col},
#'   \code{event_col}, and \code{feature_col}.
#' @param project Character; project prefix for output directory/file naming.
#'   Default \code{"KMplot"}.
#' @param event_col Character; column name for the event indicator (numeric or
#'   categorical). Numeric \code{>0} is treated as event. For non-numeric, any
#'   value in \code{event_positive} is treated as event. Default \code{"caselabel"}.
#' @param time_col Character; column name for follow-up time (must be positive
#'   and finite). Default \code{"time"}.
#' @param feature_col Character; column name of the continuous feature used to
#'   define groups. Default \code{"feature"}.
#' @param name Optional character suffix to distinguish outputs (inserted before
#'   the date in the PDF filename). Default \code{NULL}.
#' @param event_positive Character vector of labels that indicate an event when
#'   \code{event_col} is non-numeric (case-insensitive). Default
#'   \code{c("Case","ACR","Yes","1")}.
#' @param n_group Integer; number of quantile groups to form from \code{feature}.
#'   Default \code{3}.
#' @param group_labels Optional character vector of length \code{n_group} to use
#'   as group labels. If \code{NULL}, labels \code{"G1"}, \code{"G2"}, … are used.
#' @param height Numeric; PDF height (inches). Default \code{5}.
#' @param width Numeric; PDF width (inches). Default \code{6}.
#'
#' @return (Invisibly) a list with:
#' \itemize{
#'   \item \code{fit}: a \code{survfit} object,
#'   \item \code{data}: preprocessed data used for the fit,
#'   \item \code{plot}: \code{ggsurvplot} return (list with \code{$plot}, etc.),
#'   \item \code{pdf}: output PDF filepath,
#'   \item \code{p_value}: log-rank p-value (numeric).
#' }
#'
#' @details
#' Groups are created with \code{ggplot2::cut_number} (equal-count bins) over the
#' non-missing \code{feature}. The log-rank p is computed from
#' \code{survival::survdiff}. The function always creates dated output folders:
#' \verb{<project>_<YYMMDD>/} and \verb{<project>_<YYMMDD>/pdf/}, then writes
#' \code{Kaplan.Meier.plot.<feature_col>.<project>.<name?>.<YYMMDD>.pdf}.
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' df <- data.frame(
#'   time = rexp(120, 1/200),
#'   caselabel = sample(c("Case","Control"), 120, replace=TRUE),
#'   feature = rnorm(120)
#' )
#' Go_kmplot(
#'   df,
#'   project = "DemoKM",
#'   event_col = "caselabel",
#'   time_col = "time",
#'   feature_col = "feature",
#'   n_group = 3,
#'   group_labels = c("Low","Mid","High"),
#'   name = "FeatureTertiles"
#' )
#' }
#'
#' @importFrom dplyr transmute filter mutate
#' @importFrom survival Surv survfit survdiff
#' @importFrom survminer ggsurvplot
#' @importFrom ggplot2 cut_number
#' @importFrom grDevices pdf dev.off
#' @export

Go_kmplot <- function(df,
                      project      = "KMplot",
                      event_col    = "caselabel",
                      time_col     = "time",
                      feature_col  = "feature",
                      name         = NULL,
                      event_positive = c("Case","ACR","Yes","1"),
                      n_group      = 3,
                      group_labels = NULL,
                      height       = 5,
                      width        = 6) {

  # 필요한 패키지 자동 설치 & 로드
  pkgs <- c("dplyr", "survival", "survminer", "ggplot2")
  for (p in pkgs) {
    if (!requireNamespace(p, quietly = TRUE)) {
      install.packages(p, dependencies = TRUE)
    }
    suppressPackageStartupMessages(library(p, character.only = TRUE))
  }

  # 출력 경로 (항상 생성)
  out <- file.path(sprintf("%s_%s", project, format(Sys.Date(), "%y%m%d")))
  if (!file_test("-d", out)) dir.create(out, recursive = TRUE)
  out_path <- file.path(sprintf("%s_%s/pdf", project, format(Sys.Date(), "%y%m%d")))
  if (!file_test("-d", out_path)) dir.create(out_path, recursive = TRUE)

  # 데이터 준비
  df <- df %>%
    dplyr::transmute(
      time    = .data[[time_col]],
      feature = .data[[feature_col]],
      raw_ev  = .data[[event_col]]
    ) #%>%
    #dplyr::filter(is.finite(time), is.finite(feature), time > 0)

  if (is.numeric(df[[time_col]])) {
    df <- df %>%
      dplyr::filter(
        is.finite(.data[[time_col]]),
        is.finite(.data[[feature_col]]),
        .data[[time_col]] > 0
      )
  } else {
    warning(glue::glue("`{time_col}` is not numeric, skipping numeric filtering."))
  }

  # timepoint변환 V1, V2, V3 to 1,2,3
  df <- df %>%
    dplyr::mutate(
      time = if (all(grepl("^V[0-9]+$", time))) {
        as.numeric(gsub("V", "", time))
      } else {
        time
      }
    )



  # 이벤트 → status (1=이벤트, 0=검열)
  if (is.numeric(df$raw_ev)) {
    df$status <- ifelse(df$raw_ev > 0, 1L, 0L)
  } else {
    df$status <- ifelse(toupper(trimws(df$raw_ev)) %in% toupper(event_positive), 1L, 0L)
  }

  # feature 기반 그룹 (분위수 n_group개)
  if (is.null(group_labels)) group_labels <- paste0("G", seq_len(n_group))
  df <- df %>%
    dplyr::mutate(feature_group = ggplot2::cut_number(feature, n_group, labels = group_labels)) %>%
    droplevels()

  # KM 적합
  fit <- survfit(Surv(time, status) ~ feature_group, data = df)

  # Log-rank p 계산 및 라벨 포맷
  lr  <- survdiff(Surv(time, status) ~ feature_group, data = df)
  df_chi <- length(lr$n) - 1
  p_raw  <- tryCatch(1 - pchisq(lr$chisq, df = df_chi), error = function(e) NA_real_)
  p_txt  <- if (is.na(p_raw)) {
    "Log-rank p = NA"
  } else if (p_raw < 1e-4) {
    "Log-rank p < 0.0001"
  } else {
    paste0("Log-rank p = ", signif(p_raw, 3))
  }

  # 팔레트(그룹 수에 맞춰 확장)
  base_pal <- c("#2ca02c", "#1f77b4", "#ff7f0e", "#9467bd", "#8c564b", "#e377c2")
  pal <- base_pal[seq_len(min(n_group, length(base_pal)))]

  # KM 플롯
  p <- ggsurvplot(
    fit, data = df,
    conf.int      = TRUE,
    pval          = p_txt,                 # <- 여기!
    risk.table    = FALSE,
    xlab          = "Days",
    ylab          = sprintf("Event-Free Probability (%s)", event_col),
    legend.title  = sprintf("%s groups", feature_col),
    legend.labs   = paste0(feature_col, "=", levels(df$feature_group)),
    break.time.by = 50,
    xlim          = c(0, max(df$time, na.rm = TRUE)),
    palette       = pal
  )

  # 파일명: feature, project, (name), 날짜
  pdf_file <- sprintf(
    "%s/Kaplan.Meier.plot.%s.%s.%s%s.pdf",
    out_path,
    feature_col,
    project,
    ifelse(is.null(name), "", paste0(name, ".")),
    format(Sys.Date(), "%y%m%d")
  )

  # 항상 PDF 저장
  grDevices::pdf(pdf_file, height = height, width = width)
  print(p$plot)
  grDevices::dev.off()
  message("PDF saved: ", pdf_file)

  # 화면 출력(원하면 유지)
  print(p)

  invisible(list(fit = fit, data = df, plot = p, pdf = pdf_file, p_value = p_raw))
}
