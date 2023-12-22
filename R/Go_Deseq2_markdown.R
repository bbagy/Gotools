#' Differential Abundance Analysis Using DESeq2
#'
#' This function performs differential abundance analysis using DESeq2 on data contained within a phyloseq object.
#' It can handle categorical outcomes and optionally accommodate confounding variables. The function also
#' allows for aggregation at different taxonomic levels.
#'
#' @param psIN Phyloseq object containing OTU/ASV counts and associated sample data.
#' @param project Name or identifier for the project or dataset being analyzed.
#' @param taxanames Optional argument to aggregate taxa at a specified taxonomic level (e.g., "Genus", "Family").
#' @param cate.outs Vector of column names in the sample data of 'psIN' representing categorical outcome variables.
#' @param cate.conf Optional vector of column names in the sample data of 'psIN' representing categorical confounding variables.
#' @param cont.conf Optional vector of column names in the sample data of 'psIN' representing continuous confounding variables.
#' @param orders Optional ordering for the levels of the categorical variables.
#' @param name Optional name to assign to the results (useful for identification if storing multiple result sets).
#'
#' @return A list of data frames, where each data frame contains the results of the differential abundance analysis
#'         for one outcome variable. Includes statistics such as log2 fold changes, p-values, and adjusted p-values,
#'         along with taxonomic information for each feature.
#'
#' @examples
#' # Assuming 'ps' is a phyloseq object with appropriate data
#' results <- Go_Deseq2_markdown(ps, "my_project",
#'                              cate.outs = c("Group"),
#'                              cate.conf = c("Gender"),
#'                              cont.conf = NULL,
#'                              orders = NULL,
#'                              name = "Analysis1")
#'
#' @export
#' @importFrom DESeq2 DESeq phyloseq_to_deseq2 results
#' @importFrom phyloseq phyloseq tax_table sample_data otu_table aggregate_taxa
#' @importFrom dplyr filter arrange
#' @importFrom stats p.adjust

Go_Deseq2_markdown <- function(psIN,  project,
                      taxanames=NULL, 
                      cate.outs,  
                      cate.conf=NULL, 
                      cont.conf=NULL, 
                      orders=NULL,
                      name=NULL){
  library("DESeq2")
  # taxa aggregate
  if(!is.null(taxanames)){
    psIN <- aggregate_taxa(psIN, taxanames)
  }else{
    psIN <- psIN
  }
  mapping <- data.frame(sample_data(psIN))
  
  
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
  
  
  
  # start
  df_list <- list()
  res <- {}
  for (mvar in cate.outs ) {
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
      
      
      if (length(mapping.sel.na.rem[,mvar]) < 4){
          next
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
          combination <- unlist(my_comparisons[i]);combination
          basline <- combination[1]
          smvar <- combination[2]
          
          mapping.sel.cb <- subset(mapping.sel, mapping.sel[[mvar]] %in% c(basline, smvar)) # phyloseq subset은 작동을 안한다.
          
          psIN.cb <- psIN.na
          sample_data(psIN.cb) <- mapping.sel.cb
          
          
          #-- DESeq2 for phyloseq --#
          gm_mean = function(x, na.rm=TRUE){
              exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
          }
          
          ### categorical and continuous confounder control
          if (length(cate.conf) >= 1) {
              for(cate in cate.conf){
                  mapping.sel.cb[,cate] <- as.factor(mapping.sel.cb[,cate])
                  mapping.sel.cb <- mapping.sel.cb[!is.na(mapping.sel.cb[,cate]), ]
                  sample_data(psIN.cb) <- mapping.sel.cb
              }
          }
          
          if (length(cont.conf) >= 1) {
              for(cont in cont.conf){
                  mapping.sel.cb[,cont] <- as.numeric(mapping.sel.cb[,cont])
                  mapping.sel.cb <- na.omit(mapping.sel.cb[!is.na(mapping.sel.cb[,cont]), ])
                  sample_data(psIN.cb) <- mapping.sel.cb
              }
          }
          
          
          if (!is.null(cate.conf) | !is.null(cont.conf)) {
              confounder <- c(cate.conf,cont.conf)
              
              form <-as.formula(sprintf("~ %s + %s", mvar, paste(setdiff(confounder, "SampleType"), collapse="+")))
              #print(form)
              
              dds = suppressMessages(suppressWarnings(phyloseq:::phyloseq_to_deseq2(psIN.cb, form)))
          }    else {
              confounder <- NULL
              dds = suppressMessages(suppressWarnings(phyloseq:::phyloseq_to_deseq2(psIN.cb, as.formula(sprintf("~ %s", mvar)))))
              #print(sprintf("~ %s", mvar))
          }
          
          geoMeans = apply(counts(dds), 1, gm_mean)
          dds = suppressMessages(suppressWarnings(estimateSizeFactors(dds, geoMeans = geoMeans)))
          dds = suppressMessages(suppressWarnings(estimateDispersions(dds)))
          vst = getVarianceStabilizedData(dds)
          dds = suppressMessages(suppressWarnings(DESeq(dds, fitType="local")))
          resultsNames(dds)
          
          # calculation
          tmp <- as.data.frame(results(dds, contrast = c(mvar, smvar, basline)))
          tmp$taxa <- unlist(lapply(rownames(tmp), function(x) {
              tmp <- unlist(strsplit(x, ";"))
              tmp[length(tmp)]
          }))
          
          
          tmp$mvar <- mvar
          tmp$basline<-basline
          tmp$bas.count <-  sum(with(mapping.sel.cb, mapping.sel.cb[,mvar] == basline))
          tmp$smvar <- smvar
          tmp$smvar.count <-  sum(with(mapping.sel.cb, mapping.sel.cb[,mvar] == smvar))
          
          
          #-- give taxa name --#
          # Extract the taxonomy table from the phyloseq object
          tax_table <- as.data.frame(tax_table(psIN.cb))
          res <- merge(tmp, tax_table, by = 0, all.x = TRUE)
          res[res == "NA NA"] <- NA
          
          tt <- try(filled_data <- apply(res[, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")], 1, function(x) na.locf(x, na.rm = FALSE)), T)
          
          
          if(inherits(tt, "try-error")){
              filled_data <- apply(res[, c("Rank1", "Phylum", "Class", "Order", "Family", "Genus", "Species")], 1, function(x) na.locf(x, na.rm = FALSE))
          }
          
          filled_data <- t(as.data.frame(filled_data))
          tt <- try(res[, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")] <- filled_data, T)
          
          if(inherits(tt, "try-error")){
              res[, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")] <- t(filled_data)
          }
          
          
          
          if (type == "taxonomy"){
              taxaRanks <- c("Kingdon","Phylum","Class","Order","Family","Genus","Species")
              for(t in 2:length(taxaRanks)){
                  
                  if (!is.null(taxanames)) {
                      if (taxanames == taxaRanks[t-1]){
                          break
                      }
                  }
                  
                  res[,taxaRanks[t]] == "NA"
                  res[,taxaRanks[t]]<- as.character(res[,taxaRanks[t]])
                  res[,taxaRanks[t]][is.na(res[,taxaRanks[t]])] <- "__"
                  
                  mask <- res[,taxaRanks[t]] %in% c("s__", "g__", "f__", "o__", "c__", "p__", "__")
                  res[mask, taxaRanks[t]] <- ""
              }
              
              
              res$TaxaName <- paste(res$Phylum,"",res$Class,"",res$Order,"",res$Family,"",res$Genus,"",res$Species)
              
              res$Species[res$Species=="NA NA"] <- "  "
              
              
              #-- create table --#
              res <- as.data.frame(res)
              res$padj <- p.adjust(res$pvalue, method="fdr")
              res$deseq2.FDR <- ifelse(res$padj < 0.05, ifelse(sign(res$log2FoldChange)==1, "up", "down"), "NS")
              res$deseq2.P <- ifelse(res$pvalue < 0.05, ifelse(sign(res$log2FoldChange)==1, "up", "down"), "NS")
          }else{
              headers <- rownames(res)
              res$taxa <- headers
              #-- create table --#
              res <- as.data.frame(res)
              res$padj <- p.adjust(res$pvalue, method="fdr")
              res$deseq2.FDR <- ifelse(res$padj < 0.05, ifelse(sign(res$log2FoldChange)==1, "up", "down"), "NS")
              res$deseq2.P <- ifelse(res$pvalue < 0.05, ifelse(sign(res$log2FoldChange)==1, "up", "down"), "NS")
          }
          
          
          
          # get ps objectonly significant taxa
          #res.sel <- subset(res, res$ancom  == T & !(res$deseq2  == "NS"));dim(res.sel)[1]
          res.sel <- subset(res, res$deseq2.FDR  %in% c("up","down"));dim(res.sel)[1]
          #taxa_sig <- res.sel$Row.names[1:dim(res.sel)[1]]; summary(taxa_sig)
          
          #if(dim(res.sel)[1] == 0){
          #  ps.taxa.sig <- psIN.cb
          #}else{
          #  ps.taxa.sig <- prune_taxa(taxa_sig, psIN.cb)
          #  print(ps.taxa.sig)
          #}
          
          # "name definition
          if (class(name) == "function"){
              name <- NULL
          }
          
          # for changing "n" to "-"
          res$basline <- gsub("Vn","V-",res$basline)
          res$smvar <- gsub("Vn","V-",res$smvar)
          
          res <- arrange(res, res$padj)
          df_list[[i]] <- res
      }
      return(df_list)
  }
}

