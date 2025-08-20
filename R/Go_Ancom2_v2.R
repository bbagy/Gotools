#' Perform ANCOM-II Analysis on Phyloseq Data
#'
#' This function runs ANCOM-II on a phyloseq object. It supports categorical and
#' continuous variables and allows adjusting for confounders. For categorical
#' outcomes, it can iterate over level combinations and save per-contrast results.
#'
#' @param psIN A \code{phyloseq} object containing count data and sample metadata.
#' @param project Character. Project or analysis name (used in output paths).
#' @param taxanames Character or NULL. Taxonomic rank to report (e.g., \code{"Genus"}).
#' @param rand.eff Character or NULL. Random-effect column name in \code{sample_data(psIN)}
#'   (e.g., subject ID) to be modeled as a random intercept. Use \code{NULL} for no random effect.
#' @param data_type Character. One of \code{c("taxonomy","kegg","pathway","RNAseq")}.
#' @param cate.outs Character vector. Categorical outcomes to analyze.
#' @param cate.conf Character vector or NULL. Categorical confounders to adjust for.
#' @param cont.conf Character vector or NULL. Continuous confounders to adjust for.
#' @param orders Character vector or NULL. Custom factor level order for outcomes/confounders.
#' @param name Character or NULL. Extra tag appended to output filenames.
#'
#' @return Invisibly returns a list of data.frames (one per analysis), and writes CSV
#'   files containing ANCOM-II results (W-statistic, p-values, q-values/fdr, effect sizes)
#'   merged with taxonomy where applicable.
#'
#' @details
#' The function prepares the design from \code{sample_data(psIN)}, optionally sets a random
#' intercept via \code{rand.eff}, runs ANCOM-II per outcome (and per level comparison when
#' applicable), and saves tidy results under \code{<project>_<YYMMDD>/table/}.
#'
#' @examples
#' \donttest{
#' Go_Ancom2(
#'   psIN      = psIN,
#'   project   = "MyProject",
#'   data_type = "taxonomy",
#'   cate.outs = c("Treatment","Condition"),
#'   cate.conf = "Gender",
#'   cont.conf = NULL,
#'   orders    = NULL,
#'   name      = "Analysis1"
#' )
#' }
#'
#' @export

Go_Ancom2 <- function(psIN, project,
                      data_type = "other",
                      taxanames = NULL,
                      cate.outs,
                      cate.conf = NULL,
                      cont.conf = NULL,
                      rand.eff = NULL,
                      orders = NULL,
                      name = NULL){

  ########################################################
  # 0. 사전 설정: 출력 제한 해제 + 병렬 비활성(필요시)
  ########################################################
  options(max.print = 100000, width = 10000)

  # 필요하다면 병렬 비활성 (BiocParallel)
  # library(BiocParallel)
  # register(SerialParam())

  ########################################################
  # 1. 출력 디렉토리 생성
  ########################################################
  out <- file.path(sprintf("%s_%s", project, format(Sys.Date(), "%y%m%d")))
  if(!file_test("-d", out)) dir.create(out)

  out_path <- file.path(sprintf("%s_%s/table", project, format(Sys.Date(), "%y%m%d")))
  if(!file_test("-d", out_path)) dir.create(out_path)

  out_DA <- file.path(sprintf("%s_%s/table/Ancom2", project, format(Sys.Date(), "%y%m%d")))
  if(!file_test("-d", out_DA)) dir.create(out_DA)

  mapping <- data.frame(sample_data(psIN))

  ########################################################
  # 2. relative vs. absolute 자동 판단 후 변환
  ########################################################
  detect_abundance_type <- function(physeq) {
    lib_sizes <- sample_sums(physeq)
    mean_lib <- mean(lib_sizes)
    # 100 ± 0.1 → relative, 1000 초과 → absolute, 그 외 unknown
    if (abs(mean_lib - 100) < 0.1) {
      return("relative")
    } else if (mean_lib > 1000) {
      return("absolute")
    } else {
      return("unknown")
    }
  }

  abundance_type <- detect_abundance_type(psIN)
  if (abundance_type == "relative") {
    total_reads <- median(sample_sums(psIN))
    otu_table(psIN) <- otu_table(psIN) * total_reads
    message(" [INFO] The table was relative. Converted to count by median read (QMP (Quantitative Microbial Profiling).")
  } else if (abundance_type == "unknown") {
    message(" [WARN] Abundance type unknown. Proceeding as is.")
  } else {
    message(" [INFO] The table is recognized as absolute count.")
  }

  ########################################################
  # 3. 너무 긴 taxa 이름 및 중복 처리
  ########################################################
  # (1) taxa_names(psIN) 잘라내기
  short_names <- substr(taxa_names(psIN), 1, 100)
  # (2) 중복 해결
  if (any(duplicated(short_names))) {
    short_names <- make.unique(short_names)
  }
  taxa_names(psIN) <- short_names

  # tax_table() 내용도 너무 길면 잘라내기
  tax_table(psIN) <- apply(tax_table(psIN), 2, function(x) substr(x, 1, 50))

  ########################################################
  # 4. 데이터 타입 인식 (taxonomy, KEGG, RNA 등)
  ########################################################
  taxtab.col <- colnames(data.frame(tax_table(psIN)))

  if (any(grepl("Species", taxtab.col))){
    type <- "taxonomy"
    message(" [INFO] Data type = taxonomy")
  } else if (any(grepl("KO", taxtab.col))){
    type <- "kegg"
    message(" [INFO] Data type = KEGG")
  } else if (any(grepl("pathway", taxtab.col))){
    type <- "pathway"
    message(" [INFO] Data type = Pathway")
  } else if (any(grepl("symbol", taxtab.col))){
    type <- "RNAseq"
    message(" [INFO] Data type = RNAseq")
  } else {
    type <- "other"
    message(" [INFO] Data type = other (default)")
  }

  ########################################################
  # 5. Zero variance ASV 제거
  ########################################################
  detect_zero_variance <- function(physeq) {
    otu_mat <- as.matrix(otu_table(physeq))
    zero_var_taxa <- apply(otu_mat, 1, function(x) var(x) == 0)
    return(names(zero_var_taxa[zero_var_taxa]))
  }

  remove_zero_variance <- function(physeq) {
    zero_var_taxa <- detect_zero_variance(physeq)
    if (length(zero_var_taxa) > 0) {
      cat("Removing zero-variance taxa:\n", paste(zero_var_taxa, collapse="\n"), "\n")
      physeq <- prune_taxa(!(taxa_names(physeq) %in% zero_var_taxa), physeq)
    }
    return(physeq)
  }

  psIN <- remove_zero_variance(psIN)

  ########################################################
  # 6. Go_filter : 저빈도 ASV 필터(함수 미정 - 기존 코드에 있다고 가정)
  ########################################################
  # Go_filter 함수가 이미 존재한다고 가정.
  # 이 함수는 cutoff에 따라 저빈도 ASV를 제거하는 것으로 추정.

  ########################################################
  # 7. ANCOM-BC2 실행 함수 정의
  ########################################################
  run_ancombc2 <- function(data, fixed_formula, mvar, rand_formula, taxanames) {
    ancombc2(
      data = data,
      p_adj_method = "BH",
      lib_cut = 1000,  # 필요시 500 등으로 조정 가능
      fix_formula = fixed_formula,
      rand_formula = rand_formula,
      group = mvar,
      tax_level = taxanames,
      struc_zero = TRUE,
      neg_lb = TRUE,
      alpha = 0.05,
      global = TRUE,
      em_control = list(tol = 1e-5, max_iter = 100),
      verbose = FALSE  # 불필요한 출력 줄이기
    )
  }

  ########################################################
  # 8. ANCOM-BC2 분석 루프 시작
  ########################################################
  # cate.outs: 분석할 주요 범주형 변수 리스트
  for (mvar in cate.outs) {
    if (length(unique(mapping[, mvar])) == 1) {
      next
    }

    # NA 제거
    mapping.sel <- data.frame(sample_data(psIN))
    mapping.sel[mapping.sel == ""] <- "NA"
    mapping.sel.na <- mapping.sel[!is.na(mapping.sel[, mvar]), ]
    psIN.na <- prune_samples(rownames(mapping.sel[!is.na(mapping.sel[, mvar]), ]), psIN)
    mapping.sel.na.rem <- data.frame(sample_data(psIN.na))

    if (length(unique(mapping.sel.na.rem[, mvar])) == 1) next

    message(sprintf("##-- %s (total without NA: %s/%s) --##",
                    mvar, dim(mapping.sel.na.rem)[1], dim(mapping.sel)[1]))

    if (length(mapping.sel.na.rem[, mvar]) < 4) {
      message(sprintf(" [WARN] %s is removed because length(%s) < 4", mvar, length(mapping.sel.na.rem[,mvar])))
      next
    }

    # factor 변환
    if (class(mapping.sel.na.rem[, mvar]) == "character" ||
        class(mapping.sel.na.rem[, mvar]) == "integer" ||
        class(mapping.sel.na.rem[, mvar]) == "numeric") {
      mapping.sel.na.rem[, mvar] <- factor(mapping.sel.na.rem[, mvar])
      sample_data(psIN.na) <- mapping.sel.na.rem
    }

    # orders 적용
    if(!is.null(orders)) {
      mapping.sel[, mvar] <- factor(mapping.sel[, mvar], levels = intersect(orders, mapping.sel[, mvar]))
    } else {
      mapping.sel[, mvar] <- factor(mapping.sel[, mvar])
    }

    # 2개씩 비교 (combn)
    cbn <- combn(x = levels(mapping.sel[,mvar]), m = 2)
    my_comparisons <- lapply(seq_len(ncol(cbn)), function(i) cbn[, i])

    # 각 쌍별 분석
    for(i in seq_along(my_comparisons)) {
      print(my_comparisons[i])
      combination <- unlist(my_comparisons[i])
      basline <- combination[1]
      smvar <- combination[2]

      mapping.sel.cb <- subset(mapping.sel, mapping.sel[[mvar]] %in% c(basline, smvar))
      psIN.cb <- psIN.na
      sample_data(psIN.cb) <- mapping.sel.cb


      # 5 이하라면 스킵
      # 그룹별 샘플 수 계산
      bas.count <- sum(mapping.sel.cb[, mvar] == basline)
      smvar.count <- sum(mapping.sel.cb[, mvar] == smvar)

      if (bas.count <= 5 || smvar.count <= 5) {
        message(sprintf("[WARN] The sample size is unbalanced: '%s' (n=%d) vs '%s' (n=%d). Analysis may fail or be unreliable. Skipping...",
                        basline, bas.count, smvar, smvar.count))
        next
      }

      ### categorical/continuous confounder
      if (length(cate.conf) >= 1) {
        for(cate in cate.conf){
          mapping.sel.cb[, cate] <- as.factor(mapping.sel.cb[, cate])
          sample_data(psIN.cb) <- mapping.sel.cb
        }
      }
      if (length(cont.conf) >= 1) {
        for(cont in cont.conf){
          mapping.sel.cb[, cont] <- as.numeric(mapping.sel.cb[, cont])
          sample_data(psIN.cb) <- mapping.sel.cb
        }
      }

      # fixed_formula
      if (!is.null(cate.conf) | !is.null(cont.conf)) {
        confounder <- c(cate.conf, cont.conf)
        fixed_formula <- sprintf("%s + %s", mvar, paste(setdiff(confounder, "SampleType"), collapse = " + "))
      } else {
        confounder <- NULL
        fixed_formula <- mvar
      }

      # rand_formula
      if (!is.null(rand.eff)) {
        rand_formula <- formula(sprintf("~ (1 | %s)", rand.eff))
      } else {
        rand_formula <- NULL
      }

      # ancom 실행
      tt <- try(ancom.out <- run_ancombc2(psIN.cb, fixed_formula, mvar, rand_formula, taxanames), silent = TRUE)

      if (class(tt) == "try-error") {
        # 0 sum 샘플 제거
        psIN.cb1 <- prune_samples(sample_sums(psIN.cb) > 0, psIN.cb)

        cutoff <- 0.001
        increment <- 0.0005
        final_cutoff <- 0.01

        while (cutoff <= final_cutoff) {
          message(" [INFO] Trying cutoff value: ", cutoff)
          psIN.cb2 <- Go_filter(psIN.cb1, cutoff = cutoff)

          test2 <- try(run_ancombc2(psIN.cb2, fixed_formula, mvar, rand_formula, taxanames), silent = TRUE)
          if (!inherits(test2, "try-error")) {
            ancom.out <- test2
            message(" [SUCCESS] Analysis succeeded with cutoff: ", cutoff)
            break
          } else {
            message(" [FAIL] Analysis failed with cutoff: ", cutoff, " - Increasing cutoff")
            cutoff <- cutoff + increment
          }
        }

        if (cutoff > final_cutoff) {
          message(" [FAIL] Analysis failed with all attempted cutoff values up to ", final_cutoff)
          next
        } else {
          message(" [INFO] Final successful cutoff: ", cutoff)
        }
      }

      # 에러 없이 ancom.out 만들었다면 결과 처리
      res.ancom <- ancom.out$res

      # (1) lfc_* 중 Intercept를 빼고 실제 효과 열을 자동으로 하나 잡음
      lfc_cols <- grep("^lfc_", colnames(res.ancom), value = TRUE)
      lfc_cols <- setdiff(lfc_cols, "lfc_(Intercept)")
      suffix   <- sub("^lfc_", "", lfc_cols)[1]   # 예: "caselabelControl"

      # (2) suffix를 재사용해서 필요한 열을 정확히 선택
      ancom_df <- res.ancom %>%
        dplyr::select(
          taxon,
          lfc    = !!rlang::sym(paste0("lfc_", suffix)),
          se     = !!rlang::sym(paste0("se_", suffix)),
          W      = !!rlang::sym(paste0("W_",  suffix)),
          pvalue = !!rlang::sym(paste0("p_",  suffix)),
          qvalue = !!rlang::sym(paste0("q_",  suffix)),
          diff   = !!rlang::sym(paste0("diff_", suffix)),
          pass   = !!rlang::sym(paste0("passed_ss_", suffix))
        )
      rownames(ancom_df) <- ancom_df$taxon
      ancom_df$taxon <- NULL

      colnames(ancom_df) <- c("lfc_ancombc", "se_ancombc", "W_ancombc",
                              "pvalue_ancombc", "qvalue_ancombc", "diff_abn", "pass")

      ancom_df$mvar <- mvar
      ancom_df$basline <- basline
      ancom_df$bas.count <- sum(mapping.sel.cb[,mvar] == basline)
      ancom_df$smvar <- smvar
      ancom_df$smvar.count <- sum(mapping.sel.cb[,mvar] == smvar)
      ancom_df$ancom2.FDR <- ifelse(ancom_df$qvalue_ancombc < 0.05,
                                    ifelse(sign(ancom_df$lfc_ancombc) == 1, "up", "down"), "NS")
      ancom_df$ancom2.P <- ifelse(ancom_df$pvalue_ancombc < 0.05,
                                  ifelse(sign(ancom_df$lfc_ancombc) == 1, "up", "down"), "NS")

      # tax_table merge
      tax_table_df <- data.frame(tax_table(psIN.cb))
      merged_results <- merge(ancom_df, tax_table_df, by = 0, all.x = TRUE)



      merged_results[merged_results == "NA NA"] <- NA
      colnames(merged_results)[colnames(merged_results) == "Row.names"] <- "ASV"

      # pvalue 기준 정렬
      merged_results <- merged_results %>% dplyr::arrange(pvalue_ancombc)
      significant_bacteria <- merged_results[merged_results$qvalue_ancombc < 0.05, ]
      num_significant <- nrow(significant_bacteria)

      # 파일명 생성
      confounder_part <- ifelse(!is.null(confounder), ".with_confounder", "")
      rand.eff_part <- ifelse(!is.null(rand.eff), paste0(".with_random_effect_", rand.eff),"")
      name_part <- ifelse(!is.null(name), paste0(".", name), "")
      filename <- sprintf("ancom2.(%s.vs.%s).Sig%s.%s%s%s%s.%s.csv",
                          basline, smvar, num_significant, mvar,
                          confounder_part, rand.eff_part, name_part, project)

      output_path <- file.path(out_DA, filename)
      write.csv(merged_results, quote = FALSE, col.names = NA, file = output_path)

      message(sprintf(" [DONE] Saved: %s", output_path))
    }
  }
}
