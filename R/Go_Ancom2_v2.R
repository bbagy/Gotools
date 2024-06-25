#' Perform ANCOM-II Analysis on Phyloseq Data
#'
#' This function conducts Analysis of Compositions of Microbiomes (ANCOM-II) on data from a Phyloseq object.
#' It supports various types of variables including categorical and continuous, and allows for adjusting for confounders.
#' The function can handle multiple categories and generates output for each pairwise comparison.
#'
#' @param psIN Phyloseq object containing the data for analysis.
#' @param project Name of the project or analysis.
#' @param taxanames Rank of taxa.
#' @param rand.eff
#' @param data_type Type of data in the Phyloseq object ("taxonomy", "kegg", "pathway", "RNAseq").
#' @param cate.outs Categorical outcomes to be analyzed.
#' @param cate.conf Categorical confounding variables to adjust for in the analysis.
#' @param cont.conf Continuous confounding variables to adjust for in the analysis.
#' @param orders Custom order of factors in the analysis, if applicable.
#' @param name Optional name for the analysis.
#'
#' @return The function generates CSV files containing the results of the ANCOM-II analysis.
#' These include statistics like log fold change, standard error, W statistic, p-values, and adjusted p-values,
#' along with taxonomic information for each feature. Files are saved in a specified directory
#' with a naming convention that includes key details of the analysis.
#'
#' @details
#' The function first checks the type of data in the Phyloseq object.
#' It then processes categorical and continuous variables, adjusting for confounders if specified.
#' ANCOM-II analysis is performed for each variable or combination of variables.
#' The results are merged with taxonomy data and saved as CSV files.
#'
#' For categorical variables, the function considers combinations of levels for analysis.
#' The analysis takes into account the structure of microbiome data and adjusts for compositionality.
#'
#' @examples
#' # psIN is a Phyloseq object
#' # Example usage:
#' Go_Ancom2(psIN = psIN,
#'           project = "MyProject",
#'           data_type = "taxonomy",
#'           cate.outs = c("Treatment", "Condition"),
#'           cate.conf = c("Gender"),
#'           cont.conf = NULL,
#'           orders = NULL,
#'           name = "Analysis1")
#'
#' @export

Go_Ancom2 <- function(psIN,  project,
                      data_type = "other",
                      taxanames = NULL,
                      cate.outs,
                      cate.conf=NULL,
                      cont.conf=NULL,
                      rand.eff=NULL,
                      orders=NULL,
                      name=NULL){

  # out dir
  out <- file.path(sprintf("%s_%s",project, format(Sys.Date(), "%y%m%d")))
  if(!file_test("-d", out)) dir.create(out)
  out_path <- file.path(sprintf("%s_%s/table",project, format(Sys.Date(), "%y%m%d")))
  if(!file_test("-d", out_path)) dir.create(out_path)
  out_DA <- file.path(sprintf("%s_%s/table/Ancom2",project, format(Sys.Date(), "%y%m%d")))
  if(!file_test("-d", out_DA)) dir.create(out_DA)

  #out_DA.Tab <- file.path(sprintf("%s_%s/table/Ancom2/tab",project, format(Sys.Date(), "%y%m%d")))
  #if(!file_test("-d", out_DA.Tab)) dir.create(out_DA.Tab)


  mapping <- data.frame(sample_data(psIN))


  # get data tyep
  print("Check the data type")
  taxtab.col <- colnames(data.frame((tax_table(psIN))))

  if (any(grepl("Species", taxtab.col))){
    taxaTab <- data.frame(tax_table(psIN)[,"Species"])
    type <- "taxonomy"
    print(type)
  }else if(any(grepl("KO", taxtab.col))){
    taxaTab <- data.frame(tax_table(psIN)[,"KO"])
    type <- "kegg"
    print(type)
  }else if(any(grepl("pathway", taxtab.col))){
    taxaTab <- data.frame(tax_table(psIN)[,"pathway"])
    type <- "pathway"
    print(type)
  }else if(any(grepl("symbol", taxtab.col))){
    taxaTab <- data.frame(tax_table(psIN)[,"symbol"])
    type <- "RNAseq"
    print(type)
  }



  # start
  res <- {}
  for (mvar in cate.outs) {
    if (length(unique(mapping[, mvar])) == 1) {
      next
    }

    #na remove
    mapping.sel <- data.frame(sample_data(psIN))
    mapping.sel[mapping.sel==""] <- "NA"
    mapping.sel.na <- mapping.sel[!is.na(mapping.sel[,mvar]), ]
    na.count <- length(mapping.sel.na)
    psIN.na <- prune_samples(rownames(mapping.sel[!is.na(mapping.sel[,mvar]), ]), psIN)
    mapping.sel.na.rem <- data.frame(sample_data(psIN.na ))

    if (length(unique(mapping.sel.na.rem[,mvar])) == 1 )
      next

   print(sprintf("##-- %s (total without NA: %s/%s) --##",
                    mvar, dim(mapping.sel.na.rem)[1], dim(mapping.sel)[1]))

    if (length(mapping.sel.na.rem[,mvar]) < 4){
      next
      print(sprintf("%s is removed because length(%s) less than 4", mvar, length(mapping.sel.na.rem[,mvar])))
    }




    # integer control
    if (class(mapping.sel.na.rem[,mvar]) == "character"){
      mapping.sel.na.rem[,mvar] <- factor(mapping.sel.na.rem[,mvar])
      sample_data(psIN.na) <- mapping.sel.na.rem
    }
    if (class(mapping.sel.na.rem[,mvar]) == "integer" | class(mapping.sel.na.rem[,mvar]) == "numeric"){
      mapping.sel.na.rem[,mvar] <- factor(mapping.sel.na.rem[,mvar])
      sample_data(psIN.na) <- mapping.sel.na.rem
    }



    # for changing "-" to character
    mapping.sel[,mvar] <- gsub("V-","Vn",mapping.sel[,mvar])

    # combination
    if(!is.null(orders)){
     mapping.sel[,mvar] <- factor(mapping.sel[,mvar], levels = intersect(orders, mapping.sel[,mvar]))
    }else{
     mapping.sel[,mvar] <- factor(mapping.sel[,mvar])
    }

    # mapping.sel[,mvar] <- factor(mapping.sel[,mvar])
    cbn <- combn(x = levels(mapping.sel[,mvar]), m = 2)

    my_comparisons <- {}
    for(i in 1:ncol(cbn)){
      x <- cbn[,i]
      my_comparisons[[i]] <- x
    };my_comparisons

    # subset sample by combination
    for(i in 1:length(my_comparisons)){
    print(my_comparisons[i])
    combination <- unlist(my_comparisons[i]);combination
    basline <- combination[1]
    smvar <- combination[2]

    mapping.sel.cb <- subset(mapping.sel, mapping.sel[[mvar]] %in% c(basline, smvar)) # phyloseq subset은 작동을 안한다.

    psIN.cb <- psIN.na

    sample_data(psIN.cb) <- mapping.sel.cb


    ### categorical and continuous confounder control
    if (length(cate.conf) >= 1) {
      for(cate in cate.conf){
        mapping.sel.cb[,cate] <- as.factor(mapping.sel.cb[,cate])
        sample_data(psIN.cb) <- mapping.sel.cb
      }
    }

    if (length(cont.conf) >= 1) {
      for(cont in cont.conf){
        mapping.sel.cb[,cont] <- as.numeric(mapping.sel.cb[,cont])
        sample_data(psIN.cb) <- mapping.sel.cb
      }
    }

    #-- ANCOM-bc for phyloseq --#

    # confounder 설정 및 fixed_formula 생성
    if (!is.null(cate.conf) | !is.null(cont.conf)) {
      confounder <- c(cate.conf, cont.conf)
      fixed_formula <- sprintf("%s + %s", mvar, paste(setdiff(confounder, "SampleType"), collapse = " + "))
      print(fixed_formula)
    } else {
      confounder <- NULL
      fixed_formula <- mvar
    }

    # rand_formula 설정
    if (!is.null(rand.eff)) {
      rand_formula <- formula(sprintf("~ (1 | %s)", rand.eff))
      print(rand_formula)
    } else {
      rand_formula <- NULL
    }

    # tax_level 변수 설정
    taxanames <- NULL

    # ancombc2 실행 함수
    run_ancombc2 <- function(data, fixed_formula, mvar, rand_formula, taxanames) {
      ancombc2(
        data = data,
        p_adj_method = "holm",
        lib_cut = 1000,
        fix_formula = fixed_formula,
        rand_formula = rand_formula,
        group = mvar,
        tax_level = taxanames,
        struc_zero = TRUE,
        neg_lb = TRUE,
        alpha = 0.05,
        global = TRUE,
        em_control = list(tol = 1e-5, max_iter = 100)
      )
    }

    # ancombc2 실행 및 에러 처리
    tt <- try(ancom.out <- run_ancombc2(psIN.cb, fixed_formula, mvar, rand_formula, taxanames), TRUE)

    if (class(tt) == "try-error") {
      # 샘플 0인 ASV 제거
      psIN.cb1 <- prune_samples(sample_sums(psIN.cb) > 0, psIN.cb)

      # 초기 cutoff 설정
      cutoff <- 0.001
      increment <- 0.0005
      final_cutoff <- 0.01

      # cutoff를 점진적으로 증가시키면서 ancombc2 실행
      while (cutoff <= final_cutoff) {
        cat("Trying cutoff value:", cutoff, "\n")
        psIN.cb2 <- Go_filter(psIN.cb1, cutoff = cutoff)

        ancom.out <- try(run_ancombc2(psIN.cb2, fixed_formula, mvar, rand_formula, taxanames), silent = TRUE)

        if (!inherits(ancom.out, "try-error")) {
          cat("Analysis succeeded with cutoff:", cutoff, "\n")
          break
        } else {
          cat("Analysis failed with cutoff:", cutoff, " - Increasing cutoff\n")
          cutoff <- cutoff + increment
        }
      }

      if (cutoff > final_cutoff) {
        cat("Analysis failed with all attempted cutoff values up to", final_cutoff, "\n")
      } else {
        cat("Final successful cutoff:", cutoff, "\n")
        # 추가 분석 진행
      }
    }



    res.ancom = ancom.out$res
    ancom_df = res.ancom %>%
      dplyr::select(taxon, contains(mvar))

    View(res.ancom)

    rownames(ancom_df) <- ancom_df$taxon; ancom_df$taxon <- NULL

    print(dim(ancom_df))

    #rownames(df_mvar) <- df_mvar$taxon
    #names(df_mvar)[length(names(df_mvar))]<-"diff_abn"
    names(ancom_df)<- c("lfc_ancombc", "se_ancombc", "W_ancombc", "pvalue_ancombc", "qvalue_ancombc",  "diff_abn", "pass")

    ancom_df$mvar <- mvar
    ancom_df$basline<-basline
    ancom_df$bas.count <-  sum(with(mapping.sel.cb, mapping.sel.cb[,mvar] == basline))
    ancom_df$smvar <- smvar
    ancom_df$smvar.count <-  sum(with(mapping.sel.cb, mapping.sel.cb[,mvar] == smvar))
    ancom_df$ancom2.FDR <- ifelse(ancom_df$qvalue_ancombc < 0.05, ifelse(sign(ancom_df$lfc_ancombc)==1, "up", "down"), "NS")
    ancom_df$ancom2.P <- ifelse(ancom_df$pvalue_ancombc < 0.05, ifelse(sign(ancom_df$lfc_ancombc)==1, "up", "down"), "NS")

    # Extract the taxonomy table from the phyloseq object
    tax_table <- as.data.frame(tax_table(psIN.cb))
    ancom_df <- as.data.frame(ancom_df)

    merged_results <- merge(ancom_df, tax_table, by = 0, all.x = TRUE)
    merged_results[merged_results == "NA NA"] <- NA

    # Fill missing values with the last available non-missing value in each row
    tt <- try(filled_data <- apply(merged_results[, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")], 1, function(x) na.locf(x, na.rm = FALSE)), T)


    if(inherits(tt, "try-error")){
      filled_data <- apply(merged_results[, c("Rank1", "Phylum", "Class", "Order", "Family", "Genus", "Species")], 1, function(x) na.locf(x, na.rm = FALSE))
    }

    filled_data <- t(as.data.frame(filled_data))

    tt <- try(merged_results[, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")] <- filled_data, T)


    if(inherits(tt, "try-error")){
      merged_results[, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")] <- t(filled_data)
    }
    colnames(merged_results)[colnames(merged_results) == "Row.names"] <- "ASV"

    colnames(merged_results)

    merged_results <- arrange(merged_results, merged_results$pvalue_ancombc)
    significant_bacteria <- merged_results[merged_results$qvalue_ancombc < 0.05, ]
    # Counting the number of significant bacteria
    num_significant <- nrow(significant_bacteria)

    confounder_part <- ifelse(!is.null(confounder), ".with_confounder", "")
    rand.eff_part <- ifelse(!is.null(rand.eff), paste0(".with_random_effect_",rand.eff),"")
    name_part <- ifelse(!is.null(name), paste0(".", name), "")
    filename <- sprintf("ancom2.(%s.vs.%s).Sig%s.%s%s%s%s.%s.csv",
                        basline, smvar, num_significant, mvar,
                        confounder_part,rand.eff_part, name_part, project)

    output_path <- file.path(out_DA, filename)
    write.csv(merged_results, quote = FALSE, col.names = NA, file = output_path)

    }
  }
}

