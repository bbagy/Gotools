#' Perform Differential Abundance Analysis Using ALDEx2 on Phyloseq Data
#'
#' This function conducts differential abundance analysis using the ALDEx2 method on data from a Phyloseq object. It supports various types of analyses including categorical and continuous variables, and allows for adjusting for confounders. The function can handle multiple categories, and generate output for each pairwise comparison.
#'
#' @param psIN Phyloseq object containing the data for analysis.
#' @param project Name of the project or analysis.
#' @param cate.outs Categorical outcomes to be analyzed.
#' @param cont.outs Continuous outcomes to be analyzed.
#' @param cate.conf Categorical confounding variables to adjust for in the analysis.
#' @param cont.conf Continuous confounding variables to adjust for in the analysis.
#' @param orders Custom order of factors in the analysis, if applicable.
#' @param name Optional name for the analysis.
#'
#' @return The function generates CSV files containing the results of the ALDEx2 differential abundance analysis. These include statistics like effect sizes, p-values, and adjusted p-values, along with taxonomic information for each feature. Files are saved in a specified directory with a naming convention that includes key details of the analysis.
#'
#' @details
#' The function first checks the type of data (taxonomy, kegg, pathway, RNAseq) in the Phyloseq object. It then processes categorical and continuous outcomes, adjusting for confounders if specified. ALDEx2 analysis is performed for each variable or combination of variables. The results are merged with taxonomy data and saved as CSV files.
#'
#' For categorical variables, the function considers combinations of levels for analysis. For continuous variables, it performs correlation analysis. The function allows for multiple comparisons and adjusts the results using methods like GLM or t-test based on the type of variable.
#'
#' @examples
#' # psIN is a Phyloseq object
#' # Example usage:
#' Go_Aldex2(psIN = psIN,
#'           project = "MyProject",
#'           cate.outs = c("Treatment", "Condition"),
#'           cont.outs = c("Age", "BMI"),
#'           cate.conf = c("Gender"),
#'           cont.conf = NULL,
#'           orders = NULL,
#'           name = "Analysis1")
#'
#' @export

Go_Aldex2 <- function(psIN,  project,
                      cate.outs = NULL,
                      cont.outs = NULL,
                      cate.conf = NULL,
                      cont.conf = NULL,
                      orders = NULL,
                      name = NULL){

  # out dir
  out <- file.path(sprintf("%s_%s",project, format(Sys.Date(), "%y%m%d")))
  if(!file_test("-d", out)) dir.create(out)
  out_path <- file.path(sprintf("%s_%s/table",project, format(Sys.Date(), "%y%m%d")))
  if(!file_test("-d", out_path)) dir.create(out_path)
  out_DA <- file.path(sprintf("%s_%s/table/Aldex2",project, format(Sys.Date(), "%y%m%d")))
  if(!file_test("-d", out_DA)) dir.create(out_DA)

  #out_DA.Tab <- file.path(sprintf("%s_%s/table/Aldex2/tab",project, format(Sys.Date(), "%y%m%d")))
  #if(!file_test("-d", out_DA.Tab)) dir.create(out_DA.Tab)


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



  # fix outcome column types
  map <- data.frame(sample_data(psIN))
  if(!is.null(cate.outs)){
    for (cate.out in  cate.outs) {
      map[,cate.out] <- factor(map[,cate.out], levels = intersect(orders, map[,cate.out]))
    }
    outcomes <- c(cate.outs)
  }

  if(!is.null(cont.outs)){
    for (cont in  cont.outs) {
      # NA 제거

      map <- na.omit(map[!is.na(map[,cont]), ])
      map[,cont] <- as.numeric(map[[cont]])
    }
    outcomes <- c(cont.outs)
  }


  if(!is.null(cate.outs) || !is.null(cont.outs)){
    outcomes <- unique(c(cate.outs, cont.outs))
    print(outcomes)
  }




  # fix confounders column types
  if(!is.null(cate.conf)){
    for (cvar in cate.conf) {
      map[,cvar] <- factor(map[,cvar])
    }
    confounders <- c(cate.conf)
  }else{
    confounders <- NULL
  }

  if(!is.null(cont.conf)){
    for (cont in  cont.conf) {
      # NA 제거
      map[,cont] <- as.character(map[[cont]]);map[,cont]
      map[,cont][map[,cont]==""] <- "NA";map[,cont]
      map[,cont] <- as.numeric(map[[cont]])
    }
    confounders <- c(cont.outs)
  }else{
    confounders <- NULL
  }


  if(!is.null(cate.conf) || !is.null(cont.conf)){
    confounders <- unique(c(cate.conf, cont.conf))
    print(confounders)
  } else {
    confounders <- NULL
    print(confounders)
  }


  #=========================#
  #  Start data controlling #
  #=========================#
  mapping <- map

  for (mvar in outcomes) {
    holm_col_name <- NULL
    est_col_name <- NULL
    pvalue_col_name <- NULL
    aldex_df <- NULL
    if (class(mapping[,mvar]) == "character" | class(mapping[,mvar]) == "factor" ){
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

      # for changing "-" to character
      mapping.sel.na.rem[,mvar] <- gsub("V-","Vn",mapping.sel.na.rem[,mvar])

      # combination

      if(!is.null(orders)){
        mapping.sel.na.rem[,mvar] <- factor(mapping.sel.na.rem[,mvar], levels = intersect(orders, mapping.sel.na.rem[,mvar]))
      }else{
        mapping.sel.na.rem[,mvar] <- factor(mapping.sel.na.rem[,mvar])
      }



      # mapping.sel[,mvar] <- factor(mapping.sel[,mvar])
      cbn <- combn(x = levels(mapping.sel.na.rem[,mvar]), m = 2)

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

        mapping.sel.cb <- subset(mapping.sel.na.rem, mapping.sel.na.rem[[mvar]] %in% c(basline, smvar)) # phyloseq subset은 작동을 안한다.
        mapping.sel.cb[,mvar] <- factor(mapping.sel.cb[,mvar], levels = intersect(orders, mapping.sel.cb[,mvar]))
        psIN.cb <- psIN.na
        sample_data(psIN.cb) <- mapping.sel.cb

        #========#
        # Aldex2 #
        #========#

        asv_table <- t(otu_table(psIN.cb))

        # Convert to matrix if it's not already
        if(!is.matrix(asv_table)) {
          asv_matrix <- as.matrix(asv_table)
        } else {
          asv_matrix <- asv_table
        }


        # Get the unique levels of your condition
        conds1 <- as.character(mapping.sel.cb[, mvar])
        conds1[conds1 == basline] <- "0"
        conds1[conds1 == smvar] <- "1"
        conds <- factor(conds1)

        # Convert `conds` to character if it's a factor
        if(is.factor(conds)) {
          conds <- as.character(conds)
        }




        # Run ALDEx2 analysis using the 'conditions' argument

        if(!is.null(confounders)){ #========= confounding variation
          print("GLM")
          model <- "GLM"
          form <-as.formula(sprintf("~ %s + %s", mvar, paste(setdiff(confounders, "SampleType"), collapse="+")))
          print(form)

          for (i in c(mvar, confounders)){
            mapping.sel.cb <- mapping.sel.cb[!is.na(mapping.sel.cb[,i]),]
          }

          mod_matrix <- model.matrix(form, data = mapping.sel.cb);dim(mod_matrix)
          sample_names <- rownames(mapping.sel.cb)

          common_samples <- intersect(colnames(asv_matrix), sample_names)
          tt <- try(asv_matrix <- asv_matrix[, common_samples],T)

          if(class(tt)== "try-error"){
            asv_matrix <- t(asv_matrix)
            tt <- try(asv_matrix <- asv_matrix[, common_samples],T)
          }

          dim(asv_matrix)

          set.seed(1)
          tt <- try(aldex_results <- aldex(asv_matrix, mod_matrix, test = "glm"), T)

          if(inherits(tt, "try-error")){
            aldex_results <- aldex(t(asv_matrix), mod_matrix, test = "glm") # test = "t" , "kw","glm", "corr"
          }
        }else{
          print("t-test")
          model <- "t-test"
          # If you suspect the matrix has been read as containing factors, convert it:
          # asv_matrix <- data.frame(lapply(asv_matrix, function(x) if(is.factor(x)) as.numeric(levels(x))[x] else x))

          tt <- try(aldex_results <- aldex(asv_matrix, conds, test = "t"), T)
          if(inherits(tt, "try-error")){
            aldex_results <- aldex(t(asv_matrix), conds, test = "t") # test = "t" , "kw","glm", "corr"
          }
        }

        aldex_df <- as.data.frame(aldex_results)
        aldex_df$model <- model
        aldex_df$mvar <- mvar
        aldex_df$basline<-basline
        aldex_df$bas.count <-  sum(with(mapping.sel.cb, mapping.sel.cb[,mvar] == basline))
        aldex_df$smvar <- smvar
        aldex_df$smvar.count <-  sum(with(mapping.sel.cb, mapping.sel.cb[,mvar] == smvar))



        if (model == "t-test"){
          aldex_df$aldex2.FDR <- ifelse(aldex_df$we.eBH < 0.05, ifelse(sign(aldex_df$effect)==1, "up", "down"), "NS")
          aldex_df$aldex2.P <- ifelse(aldex_df$wi.ep < 0.05, ifelse(sign(aldex_df$effect)==1, "up", "down"), "NS")
          aldex_df <- aldex_df[order(aldex_df$wi.ep, decreasing = F), ]

        }else if(model == "GLM"){
          tt<- try(aldex_df$aldex2.FDR <- ifelse(aldex_df[, sprintf("%s%s.pval.holm", mvar,smvar)] < 0.05,
                                                 ifelse(sign(aldex_df[, sprintf("%s%s.Est", mvar,smvar)]) == 1, "up", "down"),
                                                 "NS"),T)


          if (class(tt) == "try-error"){
            # Create the pval column name
            # holm_col_name <- sprintf("%s%s.pval.holm", mvar, smvar)
            holm_col_name <- sprintf("%s%s.pval.padj", mvar, smvar)
            holm_col_name <- gsub("\\+", ".", holm_col_name)

            # Create the Est column name
            est_col_name <- sprintf("%s%s.Est", mvar, smvar)
            est_col_name <- gsub("\\+", ".", est_col_name)

            pvalue_col_name <- sprintf("%s%s.pval", mvar,smvar)
            pvalue_col_name <- gsub("\\+", ".", pvalue_col_name)

            # Use the modified column names in the ifelse function
            aldex_df$aldex2.FDR <- ifelse(aldex_df[, holm_col_name] < 0.05,
                                          ifelse(sign(aldex_df[, est_col_name]) == 1, "up", "down"), "NS")

            aldex_df$aldex2.P <- ifelse(aldex_df[, pvalue_col_name] < 0.05,
                                          ifelse(sign(aldex_df[, est_col_name]) == 1, "up", "down"), "NS")


            aldex_df <- aldex_df[order(aldex_df[, pvalue_col_name], decreasing = F),]

          }else(
            aldex_df$aldex2.P <- ifelse(aldex_df[, sprintf("%s%s.pval", mvar,smvar)] < 0.05,
                                        ifelse(sign(aldex_df[, sprintf("%s%s.Est", mvar,smvar)]) == 1, "up", "down"),
                                        "NS")
          )

        }


        # Extract the taxonomy table from the phyloseq object
        tax_table <- as.data.frame(tax_table(psIN.cb))
        merged_results <- merge(aldex_df, tax_table, by = 0, all.x = TRUE)
        merged_results[merged_results == "NA NA"] <- NA

        # Fill missing values with the last available non-missing value in each row

        tt <- try(filled_data <- apply(merged_results[, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")], 1, function(x) na.locf(x, na.rm = FALSE)), T)


        if(inherits(tt, "try-error")){
          filled_data <- apply(merged_results[, c("Rank1", "Phylum", "Class", "Order", "Family", "Genus", "Species")], 1, function(x) na.locf(x, na.rm = FALSE))
        }



        filled_data <- t(as.data.frame(filled_data))
        set.seed(1)
        tt <- try(merged_results[, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")] <- filled_data, T)


        if(inherits(tt, "try-error")){
          merged_results[, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")] <- t(filled_data)
        }
        colnames(merged_results)[colnames(merged_results) == "Row.names"] <- "ASV"


        # Filtering based on adjusted p-value threshold
        significant_bacteria <- merged_results[merged_results$we.eBH < 0.05, ]

        # Counting the number of significant bacteria
        num_significant <- nrow(significant_bacteria)

        if (model == "t-test"){
          merged_results <- merged_results[order(merged_results$wi.ep, decreasing = F), ]

        }else if(model == "GLM"){

          pvalue <- sprintf("%s%s.pval", mvar,smvar)
          pvalue <- gsub("\\+", ".", pvalue)
          merged_results <- merged_results[order(merged_results[, pvalue], decreasing = F), ]
        }


        confounder_prefix <- ifelse(is.null(confounders), "", "with_confounder")
        name_prefix <- ifelse(is.null(name), "", paste0(".",name))


        filename <- sprintf("%s/aldex2.(%s.vs.%s).Sig%s.%s.%s.%s%s.%s.csv",
                            out_DA, basline, smvar, num_significant, mvar, model, confounder_prefix, name_prefix, project)
        write.csv(merged_results, quote = F, col.names = NA, file = filename)
      }

    }else if(class(mapping[,mvar]) == "numeric"){ #====== aldex continuous value
      model <- "corr"
      # Na control
      psIN.na <- psIN
      mapping.na <- na.omit(mapping[!is.na(mapping[,mvar]), ])
      sample_data(psIN.na) <- mapping.na


      asv_table <- t(otu_table(psIN.na))
      # Convert to matrix if it's not already
      if(!is.matrix(asv_table)) {
        asv_matrix <- as.matrix(asv_table)
      } else {
        asv_matrix <- asv_table
      }

      cont_variable <- mapping.na[,mvar]
      sample_names <- rownames(mapping.na)

      common_samples <- intersect(colnames(asv_matrix), sample_names)
      asv_matrix <- asv_matrix[, common_samples]

      dim(asv_matrix)

      set.seed(1)

      # aldex.corr
      clr_transformed <- aldex.clr(asv_matrix)
      aldex_results <- aldex.corr(clr_transformed, cont_variable)


      # Extract the taxonomy table from the phyloseq object
      tax_table <- as.data.frame(tax_table(psIN.na))
      aldex_df <- as.data.frame(aldex_results)
      aldex_df$model <- model
      aldex_df$mvar <- mvar


      aldex_df$aldex2.kFDR <- ifelse(aldex_df$kendall.eBH < 0.05, ifelse(sign(aldex_df$kendall.etau)==1, "up", "down"), "NS")
      aldex_df$aldex2.kP <- ifelse(aldex_df$kendall.ep < 0.05, ifelse(sign(aldex_df$kendall.etau)==1, "up", "down"), "NS")

      aldex_df$aldex2.sFDR <- ifelse(aldex_df$spearman.eBH < 0.05, ifelse(sign(aldex_df$spearman.erho)==1, "up", "down"), "NS")
      aldex_df$aldex2.sP <- ifelse(aldex_df$spearman.ep < 0.05, ifelse(sign(aldex_df$spearman.erho)==1, "up", "down"), "NS")


      merged_results <- merge(aldex_df, tax_table, by = 0, all.x = TRUE)
      merged_results[merged_results == "NA NA"] <- NA

      # Fill missing values with the last available non-missing value in each row
      filled_data <- apply(merged_results[, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")], 1, function(x) na.locf(x, na.rm = FALSE))

      filled_data <- t(as.data.frame(filled_data))

      tt <- try(merged_results[, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")] <- filled_data, T)


      if(inherits(tt, "try-error")){
        merged_results[, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")] <- t(filled_data)
      }
      colnames(merged_results)[colnames(merged_results) == "Row.names"] <- "ASV"


      # Filtering based on adjusted p-value threshold
      significant_bacteria <- merged_results[merged_results$kendall.eBH < 0.05, ]

      # Counting the number of significant bacteria
      num_significant <- nrow(significant_bacteria)


      merged_results <- merged_results[order(merged_results$spearman.ep, decreasing = F), ]

      name_prefix <- ifelse(is.null(name), "", paste0( ".",name))
      filename <- sprintf("%s/aldex2.continuous.Sig%s.%s.%s.%s%s.csv",
                          out_DA, num_significant,  model, mvar, name_prefix, project)
      write.csv(merged_results, quote = FALSE, col.names = NA, file = filename)

    }
  }
}


