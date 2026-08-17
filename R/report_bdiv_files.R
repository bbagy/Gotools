# Locate beta-diversity PDFs by parsing the metric field instead of relying on
# one regular expression. Go_bdivPM joins one or more metrics with "+".
.gotools_find_bdiv_pdfs <- function(wd, metric, token, subgroup = NULL) {
  candidates <- list.files(
    path = wd,
    pattern = "^ordi\\.PCoA\\..*\\.pdf$",
    full.names = TRUE
  )
  if (length(candidates) == 0) return(character(0))

  candidate_names <- basename(candidates)
  metric_fields <- sub(
    "^ordi\\.PCoA\\.(PM\\.)?([^.]+)\\..*$",
    "\\2",
    candidate_names
  )
  has_metric <- vapply(
    strsplit(metric_fields, "+", fixed = TRUE),
    function(x) metric %in% x,
    logical(1)
  )

  report_token <- as.character(token)[1]
  if (!grepl("^\\(.*\\)$", report_token)) report_token <- paste0("(", report_token, ")")
  has_token <- grepl(report_token, candidate_names, fixed = TRUE)
  has_subgroup <- if (is.null(subgroup)) rep(TRUE, length(candidates)) else {
    grepl(as.character(subgroup)[1], candidate_names, fixed = TRUE)
  }

  candidates[has_metric & has_token & has_subgroup]
}
