
#' Microbiome Regression-Based Kernel Association Test (MiRKAT / MiRKAT-LMM)
#'
#' @param psIN Phyloseq object containing microbiome data.
#' @param project Name of the project, used for output file naming.
#' @param cate.vars List of categorical variables for MiRKAT analysis.
#' @param cate.conf Optional list of confounding variables.
#' @param orders A vector specifying the order of levels in the categorical variable.
#' @param name Optional name for the output file.
#' @param strata_var Column name to use as subject ID for MiRKAT-LMM (repeated measures).
#'   When provided and valid (each block contains >1 group level), MiRKAT_LMM() is used.
#'   When strata_var perfectly confounds mvar (cross-sectional design), falls back to
#'   standard MiRKAT() with a warning — identical logic to Go_bdivPM().
#'
#' @details
#' This function performs MiRKAT or MiRKAT-LMM on microbiome data.
#' - No strata_var / confounded strata → MiRKAT() (standard, cross-sectional)
#' - Valid strata_var (repeated measures) → MiRKAT_LMM() (accounts for within-subject correlation)
#'
#' The output CSV includes a 'method' column indicating which test was applied per comparison.
#'
#' @return CSV table of p-values saved to {project}_{date}/table/mirkat/.
#'
#' @examples
#' # Standard (cross-sectional)
#' Go_mirkat(psIN = ps, project = "Proj", cate.vars = "Group",
#'           orders = c("Control", "Treatment"))
#'
#' # With repeated measures (MiRKAT-LMM)
#' Go_mirkat(psIN = ps, project = "Proj", cate.vars = "Group",
#'           orders = c("Neg", "CD"), strata_var = "PatientID")
#'
#' @export

Go_mirkat <- function(psIN, project, cate.vars, cate.conf = NULL, orders,
                      name = NULL, strata_var = NULL) {

  # --- packages -----------------------------------------------------------
  packages <- c("data.table", "CompQuadForm", "devtools", "ecodist",
                "GUniFrac", "lme4", "MASS", "Matrix", "MiRKAT", "permute")
  for (pkg in packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
    library(pkg, character.only = TRUE)
  }

  # --- output directories -------------------------------------------------
  out       <- file.path(sprintf("%s_%s", project, format(Sys.Date(), "%y%m%d")))
  out_table <- file.path(sprintf("%s_%s/table/mirkat", project, format(Sys.Date(), "%y%m%d")))
  for (d in c(out, out_table)) if (!dir.exists(d)) dir.create(d, recursive = TRUE)

  # --- strata confounding check (same logic as Go_bdivPM) -----------------
  .is_strata_confounded <- function(id_vec, mvar_vec) {
    if (is.null(id_vec) || is.null(mvar_vec)) return(FALSE)
    id_vec   <- as.character(id_vec)
    mvar_vec <- as.character(mvar_vec)
    if (length(unique(na.omit(id_vec))) < 2) return(FALSE)
    block_n <- tapply(mvar_vec, id_vec, function(x) length(unique(na.omit(x))))
    all(block_n <= 1)
  }

  # --- main loop ----------------------------------------------------------
  mapping <- data.frame(sample_data(psIN))
  if (!is.null(cate.conf)) {
    for (cv in cate.conf) mapping[, cv] <- factor(mapping[, cv])
  }
  sample_data(psIN) <- mapping

  res <- NULL

  for (mvar in cate.vars) {
    if (length(unique(mapping[, mvar])) == 1) {
      message(sprintf("%s has only 1 variation, skipping.", mvar))
      next
    }

    mapping[, mvar] <- factor(mapping[, mvar], levels = orders)
    mapping[, mvar] <- factor(mapping[, mvar])

    group.cbn        <- combn(levels(mapping[, mvar]), m = 2)
    group_comparisons <- lapply(seq_len(ncol(group.cbn)), function(i) group.cbn[, i])

    set.seed(1)

    for (i in seq_along(group_comparisons)) {
      group.combination <- group_comparisons[[i]]
      baseline <- group.combination[1]
      smvar    <- group.combination[2]

      mapping.cbn <- subset(mapping, mapping[, mvar] %in% c(baseline, smvar))
      psIN.cbn    <- psIN
      sample_data(psIN.cbn) <- mapping.cbn
      mapping.cbn <- data.frame(sample_data(psIN.cbn))

      # --- NA removal (mvar + strata_var + cate.conf) ---------------------
      keep <- !is.na(mapping.cbn[, mvar])
      if (!is.null(strata_var) && strata_var %in% names(mapping.cbn))
        keep <- keep & !is.na(mapping.cbn[, strata_var])
      if (!is.null(cate.conf)) {
        for (cv in setdiff(cate.conf, mvar))
          if (cv %in% names(mapping.cbn)) keep <- keep & !is.na(mapping.cbn[, cv])
      }
      mapping.cbn <- mapping.cbn[keep, , drop = FALSE]
      psIN.cbn    <- prune_samples(rownames(mapping.cbn), psIN.cbn)

      if (nrow(mapping.cbn) < 3 || length(unique(mapping.cbn[, mvar])) < 2) {
        message(sprintf("Skipping %s vs %s: insufficient samples.", baseline, smvar))
        next
      }

      # --- outcome & covariate matrix -------------------------------------
      mapping.cbn[, mvar] <- factor(mapping.cbn[, mvar])
      y <- as.numeric(mapping.cbn[, mvar] == levels(mapping.cbn[, mvar])[1])

      X_mat <- NULL
      if (!is.null(cate.conf) && length(cate.conf) > 0) {
        conf_vars <- setdiff(cate.conf, mvar)
        conf_vars <- intersect(conf_vars, names(mapping.cbn))
        if (length(conf_vars) > 0) {
          X_mat <- do.call(cbind, lapply(conf_vars, function(cv)
            as.numeric(factor(mapping.cbn[, cv]))))
          colnames(X_mat) <- conf_vars
        }
      }

      # --- kernel matrices ------------------------------------------------
      otu.cbn  <- data.frame(otu_table(psIN.cbn))
      tree.cbn <- phy_tree(psIN.cbn)

      unifracs     <- GUniFrac(otu.cbn, tree.cbn, alpha = c(0, 0.5, 1))$unifracs
      D.weighted   <- unifracs[,, "d_1"]
      D.unweighted <- unifracs[,, "d_UW"]
      D.generalized <- unifracs[,, "d_0.5"]
      D.BC         <- as.matrix(vegdist(otu.cbn, method = "bray"))

      Ks <- list(
        K.weighted    = D2K(D.weighted),
        K.unweighted  = D2K(D.unweighted),
        K.generalized = D2K(D.generalized),
        K.BC          = D2K(D.BC)
      )

      # --- strata check (same as Go_bdivPM) -------------------------------
      use_lmm       <- FALSE
      strata_status <- "MiRKAT"

      if (!is.null(strata_var) && strata_var %in% names(mapping.cbn)) {
        confounded <- .is_strata_confounded(mapping.cbn[[strata_var]],
                                            mapping.cbn[, mvar])
        if (confounded) {
          warning(sprintf(
            "[%s vs %s] strata_var '%s' perfectly confounds mvar ",
            baseline, smvar, strata_var
          ), "(each block contains only one group level). ",
          "Falling back to standard MiRKAT [unconstrained].")
          strata_status <- sprintf("MiRKAT [fallback: unconstrained] (strata=%s)", strata_var)
        } else {
          use_lmm       <- TRUE
          strata_status <- sprintf("MiRKAT-LMM (strata=%s)", strata_var)

          # partial confounding note
          block_n   <- tapply(as.character(mapping.cbn[, mvar]),
                              as.character(mapping.cbn[[strata_var]]),
                              function(x) length(unique(na.omit(x))))
          n_pure <- sum(block_n <= 1)
          if (n_pure > 0)
            message(sprintf(
              "Note [%s vs %s]: %d of %d block(s) contain only one mvar level.",
              baseline, smvar, n_pure, length(block_n)
            ))
        }
      }

      # --- run test -------------------------------------------------------
      if (use_lmm) {
        tt <- try(
          result <- MiRKAT::MiRKAT_LMM(
            y      = y,
            X      = X_mat,
            id     = as.character(mapping.cbn[[strata_var]]),
            Ks     = Ks,
            method = "davies"
          ), silent = FALSE
        )
      } else {
        tt <- try(
          result <- MiRKAT::MiRKAT(
            y          = y,
            Ks         = Ks,
            X          = X_mat,
            out_type   = "D",
            method     = "davies",
            omnibus    = "permutation",
            returnKRV  = FALSE,
            returnR2   = FALSE
          ), silent = FALSE
        )
      }

      if (inherits(tt, "try-error")) {
        message(sprintf("Failed: %s vs %s", baseline, smvar))
        next
      }

      # --- format result --------------------------------------------------
      per.df   <- data.frame(t(unlist(result)), stringsAsFactors = FALSE)
      rownames(per.df) <- paste(baseline, "vs", smvar)
      per.df$mvar   <- mvar
      per.df$method <- strata_status
      if (!is.null(cate.conf) && length(cate.conf) > 0)
        per.df$cate.conf <- paste(setdiff(cate.conf, "SampleType"), collapse = "+")

      res <- rbind(res, per.df)
      print(res)
    }
  }

  out_file <- sprintf("%s/mirkat.%s.%s%s.csv",
                      out_table, project,
                      ifelse(is.null(name), "", paste0(name, ".")),
                      format(Sys.Date(), "%y%m%d"))
  write.csv(res, file = out_file, quote = FALSE, row.names = TRUE)
  message("Saved: ", out_file)
}
