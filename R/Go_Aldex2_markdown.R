#' Differential Abundance Analysis Using ALDEx2
#'
#' This function performs differential abundance analysis using the ALDEx2 tool on a given phyloseq object.
#' It allows for the analysis of both categorical and continuous outcomes, with optional inclusion of confounding factors.
#'
#' @param psIN The phyloseq object containing OTU/ASV counts and associated sample data.
#' @param cate.outs A vector of column names in the sample data of 'psIN' representing categorical outcome variables.
#' @param cont.outs A vector of column names in the sample data of 'psIN' representing continuous outcome variables.
#' @param cate.conf A vector of column names in the sample data of 'psIN' representing categorical confounding variables.
#' @param cont.conf A vector of column names in the sample data of 'psIN' representing continuous confounding variables.
#' @param orders Optional ordering for the levels of the categorical variables.
#'
#' @return A list of data frames, each representing the results of the differential abundance analysis for one outcome variable.
#'         Each data frame includes statistics like effect size, p-values, and adjusted p-values, along with taxonomic information for each feature.
#'
#' @examples
#' # Assuming 'ps' is a phyloseq object with appropriate data
#' results <- Go_Aldex2_markdown(ps,
#'                              cate.outs = c("Group"),
#'                              cont.outs = c("Age"),
#'                              cate.conf = c("Gender"),
#'                              cont.conf = NULL,
#'                              orders = c("Control", "Treatment"))
#'
#' @export
#' @importFrom ALDEx2 aldex aldex.clr aldex.corr
#' @importFrom phyloseq phyloseq tax_table sample_data otu_table
#' @importFrom dplyr filter
#' @importFrom zoo na.locf
#' @importFrom stats aov kruskal.test

Go_Aldex2_markdown <- function(psIN,
                      cate.outs = NULL,
                      cont.outs = NULL,
                      cate.conf = NULL, 
                      cont.conf = NULL, 
                      orders = NULL){
  
  # get data tyep
  taxtab.col <- colnames(data.frame((tax_table(psIN))))
  
  if (any(grepl("Species", taxtab.col))){
    taxaTab <- data.frame(tax_table(psIN)[,"Species"])
    type <- "taxonomy"
  }else if(any(grepl("KO", taxtab.col))){
    taxaTab <- data.frame(tax_table(psIN)[,"KO"])
    type <- "kegg"
  }else if(any(grepl("pathway", taxtab.col))){
    taxaTab <- data.frame(tax_table(psIN)[,"pathway"])
    type <- "pathway"
  }else if(any(grepl("symbol", taxtab.col))){
    taxaTab <- data.frame(tax_table(psIN)[,"symbol"])
    type <- "RNAseq"
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
  } else {
    confounders <- NULL
  }
  
  
  #=========================#
  #  Start data controlling #
  #=========================#
  mapping <- map
  df_list <- list()
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
              
              
              
              # Run ALDEx2 analysis using the 'conditions' argument
              
              if(!is.null(confounders)){ #========= confounding variation
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
                  model <- "t-test"
                  tt <- suppressMessages(suppressWarnings(try(aldex_results <- aldex(asv_matrix, conds, test = "t"), T)))
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
                      holm_col_name <- sprintf("%s%s.pval.holm", mvar, smvar)
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
          
          
          
      }
      df_list[[i]] <- merged_results
  }
  return(df_list)
}


