#' Perform Differential Expression Analysis Using DESeq2 on Phyloseq Data
#'
#' This function conducts differential expression analysis using DESeq2 on data from a Phyloseq object.
#' It allows for the analysis of categorical outcomes, with the option to control for categorical
#' and continuous confounding variables. The function can handle data aggregated at different taxonomic levels.
#'
#' @param psIN Phyloseq object containing the data for analysis.
#' @param project Name of the project or analysis.
#' @param taxanames Taxonomic names to aggregate the data. NULL if no aggregation is needed.
#' @param cate.outs Categorical outcomes for differential expression analysis.
#' @param cate.conf Categorical confounding variables to adjust for in the analysis.
#' @param cont.conf Continuous confounding variables to adjust for in the analysis.
#' @param orders Custom order of factors in the analysis, if applicable.
#' @param name Optional name for the analysis.
#'
#' @return Generates CSV files containing the results of the DESeq2 analysis,
#' including statistics like log fold change, p-values, and adjusted p-values,
#' along with taxonomic information for each feature. Files are saved in a specified directory.
#'
#' @details
#' The function preprocesses the Phyloseq object, potentially aggregating data at a specified taxonomic level.
#' It then performs DESeq2 analysis for each categorical outcome, adjusting for confounders if specified.
#' The results are merged with taxonomy data and saved as CSV files. The function supports pairwise comparisons
#' within categorical outcomes and handles data with or without specified confounders.
#'
#' @examples
#' # psIN is a Phyloseq object
#' # Example usage:
#' Go_Deseq2(psIN = psIN,
#'           project = "MyProject",
#'           taxanames = "Genus",
#'           cate.outs = c("Treatment", "Condition"),
#'           cate.conf = c("Gender"),
#'           cont.conf = NULL,
#'           orders = NULL,
#'           name = "Analysis1")
#'
#' @export

Go_Deseq2 <- function(psIN,  project,
                      taxanames=NULL,
                      cate.outs,
                      cate.conf=NULL,
                      cont.conf=NULL,
                      orders=NULL,
                      name=NULL){

  # out dir
  out <- file.path(sprintf("%s_%s",project, format(Sys.Date(), "%y%m%d")))
  if(!file_test("-d", out)) dir.create(out)
  out_path <- file.path(sprintf("%s_%s/table",project, format(Sys.Date(), "%y%m%d")))
  if(!file_test("-d", out_path)) dir.create(out_path)
  out_DA <- file.path(sprintf("%s_%s/table/Deseq2",project, format(Sys.Date(), "%y%m%d")))
  if(!file_test("-d", out_DA)) dir.create(out_DA)



  # taxa aggregate
  if(!is.null(taxanames)){
    psIN <- aggregate_taxa(psIN, taxanames)
  }else{
    psIN <- psIN
  }
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
      print(form)

      dds = phyloseq:::phyloseq_to_deseq2(psIN.cb, form)
    }    else {
      confounder <- NULL
      dds = phyloseq:::phyloseq_to_deseq2(psIN.cb, as.formula(sprintf("~ %s", mvar)))
      print(sprintf("~ %s", mvar))
    }

    geoMeans = apply(counts(dds), 1, gm_mean)
    dds = estimateSizeFactors(dds, geoMeans = geoMeans)
    dds = estimateDispersions(dds)
    vst = getVarianceStabilizedData(dds)
    dds = DESeq(dds, fitType="local")
    resultsNames(dds)

    # calculation
      print("pass2")
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
        filled_data <- apply(res[, c("locus_tag", "symbol")], 1, function(x) na.locf(x, na.rm = FALSE))
      }


      filled_data <- t(as.data.frame(filled_data))
      tt <- try(res[, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")] <- filled_data, T)

      if(inherits(tt, "try-error")){
        res[, c("locus_tag", "symbol")] <- t(filled_data)
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


        print("pass4")
        res$TaxaName <- paste(res$Phylum,"",res$Class,"",res$Order,"",res$Family,"",res$Genus,"",res$Species)

        res$Species[res$Species=="NA NA"] <- "  "



        print("pass5")
        #-- create table --#
        res <- as.data.frame(res)
        res$padj <- p.adjust(res$pvalue, method="fdr")
        res$deseq2.FDR <- ifelse(res$padj < 0.05, ifelse(sign(res$log2FoldChange)==1, "up", "down"), "NS")
        res$deseq2.P <- ifelse(res$pvalue < 0.05, ifelse(sign(res$log2FoldChange)==1, "up", "down"), "NS")
      }else{
        headers <- rownames(res)
        res$taxa <- headers
        print("pass5")
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


      # Construct the common part of the file name
      file_prefix <- sprintf("(%s.vs.%s).Sig%s.%s", basline, smvar, dim(res.sel)[1], mvar)

      # Optional parts
      taxa_suffix <- ifelse(is.null(taxanames), "", paste0(taxanames, "."))
      confounder_suffix <- ifelse(is.null(confounder), "", paste0("with_confounder", "."))
      name_suffix <- ifelse(is.null(name), "", paste0(name, "."))

      # Complete file names
      csv_filename <- sprintf("%s/deseq2.%s.%s%s%s%s.csv", out_DA, file_prefix, taxa_suffix, confounder_suffix, name_suffix, project)
      # rds_filename <- sprintf("%s/deseq2.%s.%s%s%s%s.rds", out_DA.ps, file_prefix, taxa_suffix, confounder_suffix, name_suffix, project)

      # Write the csv and rds files
      write.csv(res, quote = FALSE, col.names = NA, file = csv_filename)
      #saveRDS(ps.taxa.sig, rds_filename)



    }
  }
}

