#' Perform Paired Permutational Multivariate Analysis of Variance (PERMANOVA)
#'
#' This function performs pairwise PERMANOVA tests to assess differences in microbial community composition
#' between groups while optionally accounting for paired or repeated measures using stratified permutations.
#'
#' @param psIN A `phyloseq` object containing the microbiome data.
#' @param cate.vars A character vector specifying the categorical variables to test.
#' @param project A string indicating the project name (used for directory and file naming).
#' @param distance_metrics A character vector of distance metrics to use for PERMANOVA
#'   (e.g., `"bray"`, `"unifrac"`, `"jaccard"`, etc.).
#' @param cate.conf (Optional) A character vector of confounding variables to adjust for in the model.
#' @param des (Optional) A short description or label for the PERMANOVA run (used in messages and filenames).
#' @param name (Optional) A name string for labeling and output file naming.
#' @param strata_by (Optional) A string specifying a column name in the metadata that defines
#'   paired or repeated-measure structures. When provided, permutations are restricted within each
#'   level of this variable via the `strata` argument in `adonis2()`.
#'
#' @details
#' This function performs **pairwise PERMANOVA** tests between all combinations of group levels
#' within each categorical variable listed in `cate.vars`.
#' The distance matrices are computed via `Go_dist()` for each metric in `distance_metrics`,
#' and then pairwise `adonis2()` tests are applied to all two-group comparisons.
#'
#' - If `cate.conf` is specified, those variables are included as covariates in the model formula.
#' - If `strata_by` is provided, permutations are stratified within that variable (e.g., subject ID),
#'   allowing paired or repeated-measures analysis.
#' - Adjusted p-values (Bonferroni correction) are computed and included in the output table.
#'
#' Each test produces a data frame summarizing:
#' \itemize{
#'   \item `pairs` – the compared group pair (e.g., "Control vs Treatment")
#'   \item `Df`, `SumsOfSqs`, `R2`, `F.Model`, `p.value`, and `padj` (Bonferroni-adjusted)
#'   \item `distance_metric`, `mvar`, and `adjusted` (confounder info)
#' }
#'
#' Results are saved as CSV files under:
#' \code{<project>_<YYMMDD>/table/perm/pair_permanova.<project>.<des>.<name>.csv}
#'
#' @return
#' A `data.frame` containing pairwise PERMANOVA results across all tested variables and distance metrics.
#' The same results are also saved as a CSV file in the project output directory.
#'
#' @importFrom vegan adonis2
#' @importFrom phyloseq sample_data prune_samples
#' @importFrom utils write.csv
#'
#' @examples
#' \dontrun{
#' Go_pairedperm(
#'   psIN = ps,
#'   cate.vars = c("TreatmentGroup"),
#'   project = "MicrobiomeStudy",
#'   distance_metrics = c("bray", "unifrac"),
#'   cate.conf = c("Age", "BMI"),
#'   des = "TreatmentComparison",
#'   name = "TreatmentAnalysis",
#'   strata_by = "SubjectID"   # Restrict permutations within subjects
#' )
#' }
#'
#' @export

Go_pairedperm <- function(psIN, cate.vars, project, distance_metrics, cate.conf=NULL,
                          des=NULL, name=NULL, strata_by=NULL) {

  # Output directory setup
  out <- file.path(sprintf("%s_%s", project, format(Sys.Date(), "%y%m%d")))
  if (!file_test("-d", out)) dir.create(out)
  out_path <- file.path(sprintf("%s/table", out))
  if (!file_test("-d", out_path)) dir.create(out_path)
  out_perm <- file.path(sprintf("%s/perm", out_path))
  if (!file_test("-d", out_perm)) dir.create(out_perm)

  # Run header
  if (!is.null(des)) {
    print(sprintf("#--- Running Paired-PERMANOVA (%s) ---#", des))
  } else {
    print("#--- Running Paired-PERMANOVA ---#")
  }

  set.seed(1)
  mapping.sel <- data.frame(sample_data(psIN))
  res.pair <- data.frame()

  #--------------------------------------------
  # Main loop
  #--------------------------------------------
  for (mvar in cate.vars) {
    mapping.sel.na <- mapping.sel[!is.na(mapping.sel[, mvar]), ]

    if (length(unique(mapping.sel.na[, mvar])) == 1) {
      cat(sprintf("There is no group to compare in '%s' (only one level: %s)\n",
                  mvar, unique(mapping.sel[, mvar])))
      next
    }

    for (distance_metric in distance_metrics) {

      psIN.sel <- prune_samples(rownames(mapping.sel.na), psIN)
      mapping.sel.na[, mvar] <- factor(mapping.sel.na[, mvar])

      # Calculate distance
      distance <- Go_dist(psIN = psIN.sel,
                          project = project,
                          cate.vars = mvar,
                          distance_metrics = distance_metric)

      # Prepare for pairwise comparisons
      x <- as.dist(distance[[distance_metric]])
      map <- mapping.sel.na
      factors <- map[, mvar]
      co <- combn(unique(as.character(factors)), 2)

      R2 <- F.Model <- p.value <- pairs <- SumsOfSqs <- Df <- c()

      #--------------------------------------------
      # Pairwise comparisons
      #--------------------------------------------
      for (elem in 1:ncol(co)) {
        group1 <- as.character(co[1, elem])
        group2 <- as.character(co[2, elem])

        map.pair <- subset(map, map[, mvar] %in% c(group1, group2))
        x1 <- as.matrix(x)[factors %in% c(group1, group2),
                           factors %in% c(group1, group2)]

        # Skip if too few samples per group
        if (table(map.pair[, mvar])[group1] <= 2 | table(map.pair[, mvar])[group2] <= 2) {
          next
        }

        # Formula
        if (!is.null(cate.conf)) {
          form <- as.formula(sprintf("x1 ~ %s + %s", mvar,
                                     paste(setdiff(cate.conf, "SampleType"), collapse = "+")))
        } else {
          form <- as.formula(sprintf("x1 ~ %s", mvar))
        }

        print(form)

        #--------------------------------------------
        # PERMANOVA (stratified if strata_by is given)
        #--------------------------------------------
        if (!is.null(strata_by)) {
          if (strata_by %in% colnames(map.pair)) {
            permanova_result <- adonis2(form,
                                        data = map.pair,
                                        permutations = 999,
                                        by = "terms",
                                        strata = map.pair[[strata_by]])
          } else {
            warning(sprintf("strata_by = '%s' not found in map. Running without strata.", strata_by))
            permanova_result <- adonis2(form, data = map.pair,
                                        permutations = 999, by = "terms")
          }
        } else {
          permanova_result <- adonis2(form, data = map.pair,
                                      permutations = 999, by = "terms")
        }

        # Collect results
        pairs <- c(pairs, paste(group1, "vs", group2))
        Df <- c(Df, permanova_result[1, 1])
        SumsOfSqs <- c(SumsOfSqs, permanova_result[1, 2])
        R2 <- c(R2, permanova_result[1, 3])
        F.Model <- c(F.Model, permanova_result[1, 4])
        p.value <- c(p.value, permanova_result[1, 5])
      }

      pairw.res <- data.frame(
        pairs, Df, SumsOfSqs, R2, F.Model, p.value,
        distance_metric = distance_metric,
        mvar = mvar,
        adjusted = ifelse(is.null(cate.conf), "", paste(setdiff(cate.conf, "SampleType"), collapse = "+"))
      )

      res.pair <- rbind(res.pair, pairw.res)
    }
  }

  #--------------------------------------------
  # Multiple testing correction and output
  #--------------------------------------------
  res.pair$padj <- p.adjust(res.pair$p.value, method = "bonferroni")
  res.pair <- res.pair[, c("pairs", "Df", "SumsOfSqs", "R2", "F.Model",
                           "p.value", "padj", "distance_metric", "mvar", "adjusted")]

  # Output file
  write.csv(
    res.pair,
    quote = FALSE,
    col.names = NA,
    file = sprintf("%s/pair_permanova.%s.%s%s%s%s.csv",
                   out_perm,
                   project,
                   ifelse(is.null(cate.conf), "", paste(cate.conf, "adjusted.", sep = "")),
                   ifelse(is.null(des), "", paste(des, ".", sep = "")),
                   ifelse(is.null(name), "", paste(name, ".", sep = "")),
                   format(Sys.Date(), "%y%m%d"))
  )

  return(res.pair)
}
