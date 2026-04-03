#' Go_OR_fit
#'
#' Fit taxon-wise adjusted logistic regression models from a phyloseq object.
#'
#' @description
#' For a binary outcome defined in \code{sample_data(psIN)}, this function fits
#' one logistic regression model per microbial feature after optional taxonomic
#' agglomeration and abundance transformation. Each model estimates the adjusted
#' odds ratio for one taxon while controlling for user-specified covariates.
#'
#' @param psIN A \code{phyloseq} object.
#' @param outcome Character scalar; outcome column in \code{sample_data(psIN)}.
#' @param pos_class Character or numeric scalar defining the positive class.
#' @param taxrank Character taxonomic rank. Defaults to \code{"Species"}.
#'   Use \code{"ASV"} to keep the original feature level.
#' @param covars Optional character vector of adjustment covariates.
#' @param transform One of \code{"log_rel"}, \code{"relative"}, or
#'   \code{"count"}.
#' @param top_n Optional integer; keep the top N taxa by mean transformed
#'   abundance before fitting models. Defaults to \code{15}.
#' @param prevalence_min Optional prevalence filter on relative abundance
#'   presence rate. Defaults to \code{0}.
#' @param abundance_min Optional mean relative abundance filter. Defaults to
#'   \code{0}.
#' @param project Optional project name used to save the fitted object.
#' @param name Optional file tag for saving the fitted object. When both
#'   \code{project} and \code{name} are provided, the fitted object is saved to
#'   \code{<project>_YYMMDD/table/OR/}.
#' @param seed Integer seed stored for reproducibility and applied before model
#'   fitting. Defaults to \code{1234}.
#' @param na_action Currently only \code{"complete"} is supported.
#' @param verbose Logical; print progress messages.
#'
#' @return A list of class \code{"Go_OR_fit"}.
#'
#' @export
Go_OR_fit <- function(psIN,
                      outcome,
                      pos_class,
                      taxrank = "Species",
                      covars = NULL,
                      transform = c("log_rel", "relative", "count"),
                      top_n = 15,
                      prevalence_min = 0,
                      abundance_min = 0,
                      project = NULL,
                      name = NULL,
                      seed = 1234,
                      na_action = c("complete"),
                      verbose = TRUE) {

  `%||%` <- function(x, y) if (is.null(x) || length(x) == 0) y else x
  log_msg <- function(...) if (isTRUE(verbose)) message(...)
  clean_tag <- function(x) {
    x <- gsub("[^A-Za-z0-9._-]+", "_", x)
    x <- gsub("_+", "_", x)
    x <- gsub("^_|_$", "", x)
    x
  }

  transform <- match.arg(transform)
  na_action <- match.arg(na_action)
  if (!is.null(seed)) seed <- as.integer(seed)[1]

  if (!inherits(psIN, "phyloseq")) stop("`psIN` must be a phyloseq object.")
  if (!requireNamespace("phyloseq", quietly = TRUE)) stop("Package `phyloseq` is required.")
  if (!is.null(seed)) {
    set.seed(seed)
    log_msg(sprintf("[Go_OR_fit] using seed=%d", seed))
  }

  meta <- data.frame(phyloseq::sample_data(psIN), check.names = FALSE, stringsAsFactors = FALSE)
  if (!is.character(outcome) || length(outcome) != 1 || !outcome %in% names(meta)) {
    stop("`outcome` must be a sample_data column name in `psIN`.")
  }
  covars <- intersect(covars %||% character(0), names(meta))
  log_msg(sprintf("[Go_OR_fit] start: outcome='%s' pos_class='%s' taxrank='%s' n_samples=%d",
                  outcome, as.character(pos_class), taxrank, nrow(meta)))

  ps_work <- psIN
  if (!identical(taxrank, "ASV")) {
    tt <- phyloseq::tax_table(ps_work, errorIfNULL = FALSE)
    if (is.null(tt) || !(taxrank %in% colnames(tt))) {
      stop(sprintf("taxrank '%s' not found in tax_table(psIN).", taxrank))
    }
    ps_work <- phyloseq::tax_glom(ps_work, taxrank = taxrank, NArm = FALSE)
  }

  make_tax_labels <- function(ps_obj, taxrank_used) {
    taxa_ids <- phyloseq::taxa_names(ps_obj)
    tt <- phyloseq::tax_table(ps_obj, errorIfNULL = FALSE)
    if (is.null(tt)) {
      return(stats::setNames(taxa_ids, taxa_ids))
    }
    tt <- as.data.frame(tt, stringsAsFactors = FALSE)
    if (identical(taxrank_used, "ASV")) {
      genus_vals <- if ("Genus" %in% colnames(tt)) tt$Genus else rep(NA_character_, nrow(tt))
      species_vals <- if ("Species" %in% colnames(tt)) tt$Species else rep(NA_character_, nrow(tt))
      family_vals <- if ("Family" %in% colnames(tt)) tt$Family else rep(NA_character_, nrow(tt))
      lab <- ifelse(
        !is.na(species_vals) & species_vals != "",
        paste(genus_vals, species_vals),
        ifelse(!is.na(genus_vals) & genus_vals != "",
               paste(genus_vals, "sp."),
               ifelse(!is.na(family_vals) & family_vals != "",
                      paste(family_vals, "taxon"),
                      taxa_ids))
      )
      lab[is.na(lab) | lab == ""] <- taxa_ids[is.na(lab) | lab == ""]
      lab <- make.unique(as.character(lab))
      return(stats::setNames(lab, taxa_ids))
    }
    if ("Species" %in% colnames(tt) && "Genus" %in% colnames(tt)) {
      species_vals <- tt$Species
      genus_vals <- tt$Genus
      lab <- ifelse(!is.na(species_vals) & species_vals != "",
                    paste(genus_vals, species_vals),
                    paste(genus_vals, "sp."))
      lab[is.na(lab) | lab == ""] <- taxa_ids[is.na(lab) | lab == ""]
      lab <- make.unique(lab)
      return(stats::setNames(lab, taxa_ids))
    }
    vals <- tt[[taxrank_used]]
    vals[is.na(vals) | vals == ""] <- taxa_ids[is.na(vals) | vals == ""]
    vals <- make.unique(as.character(vals))
    stats::setNames(vals, taxa_ids)
  }

  otu_mat <- as(phyloseq::otu_table(ps_work), "matrix")
  if (phyloseq::taxa_are_rows(ps_work)) otu_mat <- t(otu_mat)
  if (!nrow(otu_mat) || !ncol(otu_mat)) stop("No OTU data found in `psIN`.")
  log_msg(sprintf("[Go_OR_fit] OTU matrix: n_samples=%d n_features=%d", nrow(otu_mat), ncol(otu_mat)))

  relab <- sweep(otu_mat, 1, pmax(1e-12, rowSums(otu_mat)), "/")
  keep_taxa <- rep(TRUE, ncol(relab))
  if (prevalence_min > 0) keep_taxa <- keep_taxa & (colMeans(relab > 0) >= prevalence_min)
  if (abundance_min > 0) keep_taxa <- keep_taxa & (colMeans(relab) >= abundance_min)
  relab <- relab[, keep_taxa, drop = FALSE]
  otu_mat <- otu_mat[, keep_taxa, drop = FALSE]
  if (!ncol(relab)) stop("No taxa remain after prevalence/abundance filtering.")
  log_msg(sprintf("[Go_OR_fit] after prevalence/abundance filter: n_features=%d", ncol(relab)))

  feat_micro <- switch(
    transform,
    log_rel = log(relab * 100 + 1),
    relative = relab,
    count = otu_mat
  )

  if (!is.null(top_n)) {
    top_n <- min(as.integer(top_n), ncol(feat_micro))
    keep_top <- order(colMeans(feat_micro), decreasing = TRUE)[seq_len(top_n)]
    feat_micro <- feat_micro[, keep_top, drop = FALSE]
  }
  log_msg(sprintf("[Go_OR_fit] after top_n: n_features=%d", ncol(feat_micro)))

  label_map <- make_tax_labels(ps_work, taxrank)
  feature_names_raw <- colnames(feat_micro)
  feature_names_lab <- unname(label_map[feature_names_raw])
  feature_names_lab[is.na(feature_names_lab) | feature_names_lab == ""] <- feature_names_raw[is.na(feature_names_lab) | feature_names_lab == ""]
  feature_names_lab <- make.unique(feature_names_lab)
  feature_names_safe <- make.names(feature_names_lab, unique = TRUE)
  label_df <- data.frame(
    feature_raw = feature_names_raw,
    feature = feature_names_safe,
    feature_label = feature_names_lab,
    stringsAsFactors = FALSE
  )

  feat_df <- data.frame(feat_micro, check.names = FALSE, stringsAsFactors = FALSE)
  colnames(feat_df) <- feature_names_safe
  feat_df$.SampleID <- rownames(feat_df)

  meta_use <- meta[, unique(c(outcome, covars)), drop = FALSE]
  meta_use$.SampleID <- rownames(meta_use)
  work_df <- merge(meta_use, feat_df, by = ".SampleID", all = FALSE, sort = FALSE)
  rownames(work_df) <- work_df$.SampleID

  if (na_action == "complete") work_df <- stats::na.omit(work_df)
  if (!nrow(work_df)) stop("No complete rows available after NA filtering.")
  log_msg(sprintf("[Go_OR_fit] after NA filter: n_samples=%d", nrow(work_df)))

  y_raw <- work_df[[outcome]]
  outcome_vals <- as.character(y_raw)
  y_model <- as.integer(outcome_vals == as.character(pos_class))
  if (length(unique(y_model)) < 2) {
    stop("Both positive and negative classes must be present in the analysis data.")
  }
  neg_levels <- sort(unique(outcome_vals[outcome_vals != as.character(pos_class)]))
  negative_class <- if (length(neg_levels)) paste(neg_levels, collapse = ", ") else "others"

  safe_covars <- intersect(covars, colnames(work_df))
  safe_features <- feature_names_safe[feature_names_safe %in% colnames(work_df)]

  calc_or_row <- function(feature_nm) {
    d <- data.frame(
      y = y_model,
      taxon = work_df[[feature_nm]],
      work_df[, safe_covars, drop = FALSE],
      check.names = FALSE
    )
    rhs <- c("taxon", safe_covars)
    fm <- stats::as.formula(paste("y ~", paste(rhs, collapse = " + ")))
    fit <- tryCatch(
      stats::glm(fm, data = d, family = stats::binomial()),
      error = function(e) NULL
    )
    if (is.null(fit)) return(NULL)
    coef_tab <- summary(fit)$coefficients
    if (!"taxon" %in% rownames(coef_tab)) return(NULL)
    est <- unname(coef_tab["taxon", "Estimate"])
    se <- unname(coef_tab["taxon", "Std. Error"])
    pval <- unname(coef_tab["taxon", "Pr(>|z|)"])
    if (any(is.na(c(est, se, pval)))) return(NULL)
    data.frame(
      feature = feature_nm,
      estimate = est,
      OR = exp(est),
      LCI = exp(est - 1.96 * se),
      UCI = exp(est + 1.96 * se),
      p.value = pval,
      stringsAsFactors = FALSE
    )
  }

  results <- do.call(
    rbind,
    lapply(safe_features, calc_or_row)
  )
  if (is.null(results) || !nrow(results)) stop("No valid logistic regression results were produced.")

  results <- merge(results, label_df[, c("feature", "feature_label")], by = "feature", all.x = TRUE, sort = FALSE)
  results$feature_label[is.na(results$feature_label) | results$feature_label == ""] <- results$feature[is.na(results$feature_label) | results$feature_label == ""]
  results$p.adj <- stats::p.adjust(results$p.value, method = "BH")
  results$sig_fdr     <- results$p.adj   < 0.05
  results$sig_nominal <- results$p.value < 0.05 & !results$sig_fdr
  results$sig         <- results$sig_fdr | results$sig_nominal
  results <- results[order(results$OR, decreasing = FALSE), , drop = FALSE]
  rownames(results) <- NULL

  subtitle <- paste0(
    "Model: adjusted OR | taxrank = ", taxrank,
    " | transform = ", transform,
    " | outcome = ", as.character(pos_class), " vs ", negative_class,
    " | covars = ", if (length(safe_covars)) paste(safe_covars, collapse = ",") else "none",
    " | features = ", nrow(results),
    " | n = ", nrow(work_df)
  )

  out <- list(
    results = results,
    model_info = list(
      engine = "logistic_or",
      family = "binomial",
      outcome = outcome,
      positive_class = as.character(pos_class),
      negative_class = negative_class,
      covars = safe_covars,
      taxrank = taxrank,
      transform = transform,
      top_n = top_n,
      n = nrow(work_df),
      seed = seed
    ),
    subtitle = subtitle,
    data_used = work_df
  )

  class(out) <- c("Go_OR_fit", class(out))
  if (!is.null(project) && !is.null(name)) {
    out_dirs <- Go_path(project = project, pdf = "no", table = "yes", path = NULL)
    out_dir <- file.path(out_dirs$tab, "OR")
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    save_path <- file.path(out_dir, paste0(clean_tag(name), ".rds"))
    saveRDS(out, save_path)
    out$saved_rds <- save_path
    log_msg(sprintf("[Go_OR_fit] saved RDS: %s", save_path))
  }

  out
}
