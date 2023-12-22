
#' Perform Regression Differential Abundance Analysis
#'
#' This function conducts regression differential abundance analysis on microbiome data using linear mixed-effects models or generalized linear models.
#'
#' @param psIN Phyloseq object containing the data for analysis.
#' @param cate.outs Categorical outcomes for regression analysis.
#' @param con.outs Continuous outcomes for regression analysis.
#' @param cate.conf Confounding categorical variables to adjust in the model.
#' @param con.conf Confounding continuous variables to adjust in the model.
#' @param StudyID Column name identifying individual study participants or samples.
#' @param project Name of the project or analysis.
#' @param pval P-value threshold for significance in the analysis.
#' @param orders Order of levels in the categorical variables.
#' @param nsamps_threshold Threshold for minimum sample abundance.
#' @param filt_threshold Threshold for filtering taxa based on sample prevalence.
#' @param taxanames Taxonomic rank to aggregate (e.g., "Genus", "Species").
#' @param name Optional name for the analysis.
#'
#' @return Generates tables and RDS files for differential abundance results and saves them in specified directories.
#'
#' @details
#' The function applies linear mixed-effects models (LMEMs) or generalized linear models (GLMs) to identify taxa that are differentially abundant across the specified outcomes while adjusting for confounders.
#'
#' @examples
#' # Example usage:
#' Go_regDA(psIN = phyloseq_object,
#'          cate.outs = c("Group1", "Group2"),
#'          con.outs = c("Age", "BMI"),
#'          cate.conf = c("Sex", "Diet"),
#'          con.conf = c("Time", "Dosage"),
#'          StudyID = "SubjectID",
#'          project = "MyProject",
#'          pval = 0.05,
#'          orders = c("Control", "Treatment"),
#'          nsamps_threshold = 0.01,
#'          filt_threshold = 0.1,
#'          taxanames = "Genus",
#'          name = "Analysis1")
#'
#' @export

Go_regDA <- function(psIN,
                    cate.outs=NULL, 
                    con.outs=NULL,
                    cate.conf=NULL, 
                    con.conf = NULL,
                    StudyID, 
                    project, 
                    pval=0.05,
                    orders, 
                    nsamps_threshold, 
                    filt_threshold, 
                    taxanames, 
                    name){

  # out dir
  out <- file.path(sprintf("%s_%s",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out)) dir.create(out)
  out_path <- file.path(sprintf("%s_%s/table",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out_path)) dir.create(out_path)
  out_table <- file.path(sprintf("%s_%s/table/lmem",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out_table)) dir.create(out_table)

  out_lmem.Tab <- file.path(sprintf("%s_%s/table/lmem/tab",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out_lmem.Tab)) dir.create(out_lmem.Tab)
  
  out_lmem.ps <- file.path(sprintf("%s_%s/table/lmem/ps",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out_lmem.ps)) dir.create(out_lmem.ps)
  
  # Check if taxanames is NULL or empty before proceeding
  if (!is.null(taxanames) && length(taxanames) > 0) {
    if (taxanames == "ASVs") {
      taxanames <- NULL
    }
  } 
  
  if(!is.null(taxanames)){
    psIN.agg <- aggregate_taxa(psIN, taxanames);psIN.agg
  }else{
    psIN.agg <- psIN
  }

  map <- data.frame(sample_data(psIN.agg))

  # fix outcome column types
  if(!is.null(cate.outs)){
    for (cate.out in  cate.outs) {
      map[,cate.out] <- factor(map[,cate.out], levels = intersect(orders, map[,cate.out]))
    }
    outcomes <- c(cate.outs)
  } 
  
  if(!is.null(con.outs)){
    for (con.out in  con.outs) {
      # NA 제거
      map[,con.out] <- as.character(map[[con.out]]);map[,con.out]
      map[,con.out][map[,con.out]==""] <- "NA";map[,con.out]
      map[,con.out] <- as.numeric(map[[con.out]])
    }
    outcomes <- c(con.outs)
  }
  
  if(!is.null(cate.outs) & !is.null(con.outs)){
    outcomes <- unique(c(cate.outs, con.outs))
  }
  
  
  # fix confounders column types
  if(!is.null(cate.conf)){
    for (cvar in cate.conf) {
      map[, cvar] <- map[df$SampleID, cvar]
    }
    confounders <- c(cate.conf)
  }else{
    confounders <- NULL
  }

  if(!is.null(con.conf)){
    for (con.out in  con.conf) {
      # NA 제거
      map[,con.conf] <- as.character(map[[con.conf]]);map[,con.conf]
      map[,con.conf][map[,con.conf]==""] <- "NA";map[,con.conf]
      map[,con.conf] <- as.numeric(map[[con.conf]])
    }
    confounders <- c(con.outs)
  }else{
    confounders <- NULL
  }
  
  if(!is.null(cate.conf) & !is.null(con.conf)){
    confounders <- unique(c(cate.outs, con.outs))
  }else{
    confounders <- NULL
  }
  
  
  
  #===============#
  # RUN the Model #
  #===============#
  for (mvar in outcomes){
    
    if (!is.null(taxanames) && length(taxanames) > 0) {
      if (taxanames == "ASVs") {
        taxanames <- NULL
      }
    } 
    
    if(!is.null(taxanames)){
      otu.filt <- as.data.frame(t(otu_table(psIN.cb)))
      tt <- try(otu.filt[,taxanames]  <- getTaxonomy(otus=rownames(otu.filt), taxRanks = colnames(tax_table(psIN.cb)), tax_tab=tax_table(psIN.cb), level=taxanames),T)
      
      if(class(tt) == "try-error"){
        print("other table")
        otu.filt <- as.data.frame(otu_table(psIN.cb)) 
        otu.filt[,taxanames] <- getTaxonomy(otus=rownames(otu.filt), tax_tab=tax_table(psIN.cb), taxRanks=colnames(tax_table(psIN.cb)),level=taxanames)
      }else{
        otu.filt <- as.data.frame(t(otu_table(psIN.cb)))
        print("DADA2 table")
        otu.filt[,taxanames] <- getTaxonomy(otus=rownames(otu.filt), tax_tab=tax_table(psIN.cb), taxRanks=colnames(tax_table(psIN.cb)),level=taxanames)
      }
    }else{
      taxanames <- "ASVs"
      otu.filt <- as.data.frame(t(otu_table(psIN)))
      otu.filt$ASVs <- rownames(otu.filt)
    }
    
    agg <- aggregate(as.formula(sprintf(". ~ %s" , taxanames)), otu.filt, sum, na.action=na.pass)
    genera <- agg[,taxanames]
    
    agg <- agg[,-1]
    #agg <- normalizeByCols(agg)
    rownames(agg) <- genera
    dim(agg)
    ftk <- names(which(unlist(apply(agg, 1, function(x) length(which(x>=nsamps_threshold)))) > ceiling(filt_threshold*ncol(agg))))
    agg <- agg[intersect(ftk,ftk),]
    # control data set after filter
    if (dim(agg)[1] == 0)
      next
    
    agg[,taxanames] <- rownames(agg)
    
    print("1. Taxa table preparation")
    
    if (class(map[,mvar]) == "factor"){ # ==== catagorical value
      # combination
      mapping.sel <- data.frame(sample_data(psIN.agg))
      mapping.sel[,mvar] <- factor(mapping.sel[,mvar], levels = orders)
      
      mapping.sel[,mvar] <- factor(mapping.sel[,mvar])
      cbn <- combn(x = levels(mapping.sel[,mvar]), m = 2)
      
      my_comparisons <- {}
      for(i in 1:ncol(cbn)){
        x <- cbn[,i]
        my_comparisons[[i]] <- x
      };my_comparisons
      
      
      for(i in 1:length(my_comparisons)){
        print(my_comparisons[i])
        combination <- unlist(my_comparisons[i]);combination
        baseline <-combination[1];baseline
        smvar <- combination[2];smvar
        
        mapping.sel.cb <- subset(mapping.sel, mapping.sel[[mvar]] %in% c(baseline, smvar));dim(mapping.sel.cb) # phyloseq subset은 작동을 안한다.
        psIN.cb <- psIN.agg
        sample_data(psIN.cb) <- mapping.sel.cb
        
#--------------    model    -------------#
        res <- {}
        for (f in agg[,taxanames]) {
          # clean bacteria name
          if (f == "s__" || f == "g__" || f == "f__" || f == "o__" || f == "c__"|| f == "p__"){
            next
          }
          
          df <- melt(agg[f,]); colnames(df) <- c("Genus", "SampleID", "value"); df$SampleID <- as.character(df$SampleID)
          df$StudyID <- mapping.sel.cb[df$SampleID, StudyID]
          
          # add groups
          for (cate in cate.conf) {
            df$Group <- as.character(mapping.sel[df$SampleID, cate])
            df[,cate] <- mapping.sel[df$SampleID, cate]
            
            # order
            if (length(orders) >= 1) {
              df[,mvar] <- factor(df[,cate], levels = orders)
            }
            else {
              df[,mvar] <- factor(df[,cate])
            }
          }
          
          # na remove
          mapping.sel.cb[mapping.sel.cb==""] <- "NA"
          mapping.na <- mapping.sel.cb[!is.na(mapping.sel.cb[,mvar]), ]
          na.count <- length(mapping.na)
          if (length(unique(mapping.na[,mvar])) == 1)
            next
          
          # na count
          print(sprintf("##-- %s (total without NA: %s/%s) --##", mvar, dim(mapping.na)[1], dim(mapping.sel.cb)[1]))
          df[,mvar] <- mapping.na[df$SampleID, mvar]
          print("2. Matching taxa table and mapping file")
          
          #=====================#
          #  Regression method  #
          #=====================#
          if(!is.null(StudyID)){
            reg <- "LMEM"
            form <- as.formula(sprintf("value ~ (1 | StudyID) + %s  %s", mvar, 
                                       ifelse(is.null(confounders), "", paste("+",setdiff(confounders, "SampleType"), collapse=""))))
            print(form)
            tt <- try(mod <- lmer(form, data=df, control=lmerControl(check.nobs.vs.nlev = "ignore",check.nobs.vs.rankZ = "ignore",check.nobs.vs.nRE="ignore")),T)
            if(class(tt) == "try-error"){
              next
            }else{
              mod <- lmer(form, data=df, control=lmerControl(check.nobs.vs.nlev = "ignore",check.nobs.vs.rankZ = "ignore",check.nobs.vs.nRE="ignore"))
              # mod2 <- lmer(form, data=df)
            }
          }else{
            reg <- "GLM"
            form <- as.formula(sprintf("value ~ %s  %s", mvar, 
                                       ifelse(is.null(confounders), "", paste("+",setdiff(confounders, "SampleType"), collapse=""))));form
            print(form)
            #mod <- glm(form, data=df,  family = binomial(link='logit'))
            tt <- try(mod <- glm(form, data=df,  family = quasipoisson(link = "log")),T)
            
            if("try-error" %in% class(tt)){
              mod <- glm(form, data=df, family = gaussian(link = "identity")) # data contained negative value
              reg <- "GLM_gaussian"
            } else {
              mod <- glm(form, data=df, family = quasipoisson(link = "log")) # for overfitting when use poisson
              reg <- "GLM_quasipoisson"
            }
            #mod <- glm(form, data=df,  gaussian())
          }
          
          coef <- summary(mod)$coefficients
          coef <- coef[grep(mvar, rownames(coef)),,drop=F]
          res <- rbind(res, cbind(f, mvar, rownames(coef), coef, baseline))
          dim(res)
        }
        print("check")
        
        
        print("3. Test Model")
        #-- create table --#
        res <- as.data.frame(res)
        #colnames(res) <- c("taxa", "metadata", "coefficient", "Estimate", "SE", "df", "t", "pvalue", "baseline")
        print(5)
        
        if(!is.null(StudyID)){
          reg <- "LMEM"
          res$pvalue <- as.numeric(as.character(res$`Pr(>|t|)`))
          
        }else{
          #reg <- "GLM"
          
          tt <- try(res$pvalue <- as.numeric(as.character(res$`Pr(>|z|)`)),T)
          
          if(class(tt) == "try-error"){
            res$pvalue <- as.numeric(as.character(res$`Pr(>|t|)`))
            res$`Pr(>|t|)` <- NULL
          }else{
            res$pvalue <- as.numeric(as.character(res$`Pr(>|z|)`))
            res$`Pr(>|z|)` <- NULL
          }
        }
        
        res$Estimate <- as.numeric(as.character(res$Estimate))
        res$SE <- as.numeric(as.character(res$`Std. Error`))
        res$padj <- p.adjust(res$pvalue, method="fdr")
        res$Model <- reg
        res <- res[order(res$pvalue),]
        
        
        # Make sure tax_table is in the correct format for the merge function
        tax_table <- as(tax_table(psIN.agg), "matrix")
        rownames(res) <- res$f
        res$f <- NULL
        common_otus <- intersect(rownames(res), rownames(tax_table))
        
        # Subset the aldex_diff_abundance and tax_table based on common_otus
        res <- res[rownames(res) %in% common_otus, ]; head(res)
        tax_table <- tax_table[rownames(tax_table) %in% common_otus, ]
        res <- merge(res, tax_table, by.x = "row.names", by.y = "row.names", all.x = TRUE)
        
        res.sel <- res
        res.sel <- as.data.frame(subset(res, pvalue < pval))
        taxa_sig <- res.sel$taxa[1:dim(res.sel)[1]]; summary(taxa_sig)
        
        res.sel <- res
        res.sel <- as.data.frame(subset(res, pvalue < pval))
        taxa_sig <- res.sel$taxa[1:dim(res.sel)[1]]; summary(taxa_sig)
        
        print("4. Prepare the table")
        
        if(dim(res.sel)[1] == 0){
          next
        }else{
          res.sel$bas.count <-  unique(sum(with(mapping.na, mapping.na[,mvar] == baseline)))
          res.sel$coef.count <-  unique(sum(with(mapping.na, mapping.na[,mvar] == smvar)))
        }
        
        if(dim(res.sel)[1] == 0){
          ps.taxa.sig <- psIN.cb
        }else{
          tt <- try(ps.taxa.sig <- prune_taxa(taxa_sig, psIN.cb),T)
          
          if(class(tt) == "try-error"){
            pathwayTab <- data.frame(otu_table(psIN.cb))
            pathwayRank <- data.frame(tax_table(psIN.cb))
            rownames(pathwayRank) <- pathwayRank[,taxRanks]
            rownames(pathwayTab) <- pathwayRank[,taxRanks]
            pathwayRank <- as.matrix(pathwayRank)
            pathwayTab <- as.matrix(t(pathwayTab))
            psIN.cb <- phyloseq(otu_table(pathwayTab, taxa_are_rows=FALSE), tax_table(pathwayRank));psIN.cb
            ps.taxa.sig <- prune_taxa(taxa_sig, psIN.cb)
            print(ps.taxa.sig)
          }else{
            ps.taxa.sig <- prune_taxa(taxa_sig, psIN.cb)
            print(ps.taxa.sig)
          }
        }
        
        
        res.sel.sel <- subset(res.sel, Genus != "NA")
        
        write.csv(res.sel.sel, quote = FALSE,col.names = NA,file=sprintf("%s/Categorical.outs.(%s.vs.%s).Sig%s.%s.%s%s.%s%s.csv",out_lmem.Tab,
                                                                         baseline, 
                                                                         smvar,
                                                                         dim(res.sel)[1],
                                                                         mvar, 
                                                                         ifelse(is.null(taxanames), "", paste(taxanames, ".", sep = "")), 
                                                                         project,
                                                                         ifelse(is.null(name), "", paste(name, ".", sep = "")), 
                                                                         reg, 
                                                                         sep="/"))
        saveRDS(ps.taxa.sig,sprintf("%s/Categorical.outs.(%s.vs.%s).Sig%s.%s.%s%s.%s%s.rds",out_lmem.ps,
                                    baseline, 
                                    smvar,
                                    dim(res.sel)[1],
                                    mvar, 
                                    ifelse(is.null(taxanames), "", paste(taxanames, ".", sep = "")), 
                                    project,
                                    ifelse(is.null(name), "", paste(name, ".", sep = "")), 
                                    reg))
      }
    } else { #====== continuous outcomes
      print("continuous outcomes")
      #--------------    lmer    -------------#
      res <- {}
      for (f in agg[,taxanames]) {
        # clean bacteria name
        if (f == "s__" || f == "g__" || f == "f__" || f == "o__" || f == "c__"|| f == "p__"){
          next
        }
        
        df <- melt(agg[f,]); colnames(df) <- c("Genus", "SampleID", "value"); df$SampleID <- as.character(df$SampleID)
        df$StudyID <- mapping.sel.cb[df$SampleID, StudyID]
        
        # add groups
        for (cate in cate.conf) {
          df$Group <- as.character(mapping.sel[df$SampleID, cate])
          df[,cate] <- mapping.sel[df$SampleID, cate]
          
          # order
          if (length(orders) >= 1) {
            df[,mvar] <- factor(df[,cate], levels = orders)
          }
          else {
            df[,mvar] <- factor(df[,cate])
          }
        }
        
        # na remove
        mapping.sel.cb[mapping.sel.cb==""] <- "NA"
        mapping.na <- mapping.sel.cb[!is.na(mapping.sel.cb[,mvar]), ]
        na.count <- length(mapping.na)
        if (length(unique(mapping.na[,mvar])) == 1)
          next
        
        # na count
        print(sprintf("##-- %s (total without NA: %s/%s) --##", mvar, dim(mapping.na)[1], dim(mapping.sel.cb)[1]))
        df[,mvar] <- mapping.na[df$SampleID, mvar]
        print("2. Matching taxa table and mapping file")
        
        #=====================#
        #  Regression method  #
        #=====================#
        if(!is.null(StudyID)){
          reg <- "LMEM"
          form <- as.formula(sprintf("value ~ (1 | StudyID) + %s  %s", mvar, 
                                     ifelse(is.null(confounders), "", paste("+",setdiff(confounders, "SampleType"), collapse=""))))
          print(form)
          tt <- try(mod <- lmer(form, data=df, control=lmerControl(check.nobs.vs.nlev = "ignore",check.nobs.vs.rankZ = "ignore",check.nobs.vs.nRE="ignore")),T)
          if(class(tt) == "try-error"){
            next
          }else{
            mod <- lmer(form, data=df, control=lmerControl(check.nobs.vs.nlev = "ignore",check.nobs.vs.rankZ = "ignore",check.nobs.vs.nRE="ignore"))
            # mod2 <- lmer(form, data=df)
          }
        }else{
          reg <- "GLM"
          form <- as.formula(sprintf("value ~ %s  %s", mvar, 
                                     ifelse(is.null(confounders), "", paste("+",setdiff(confounders, "SampleType"), collapse=""))));form
          print(form)
          #mod <- glm(form, data=df,  family = binomial(link='logit'))
          tt <- try( mod <- glm(form, data=df,  family = quasipoisson(link = "log")),T)
          if("try-error" %in% class(tt)){
            mod <- glm(form, data=df, family = gaussian(link = "identity")) # data contained negative value
            reg <- "GLM_gaussian"
          } else {
            mod <- glm(form, data=df, family = quasipoisson(link = "log")) # for overfitting when use poisson
            reg <- "GLM_quasipoisson"
          }
          #mod <- glm(form, data=df,  family = poisson(link = "log")
        }
        
        coef <- summary(mod)$coefficients
        coef <- coef[grep(mvar, rownames(coef)),,drop=F]
        res <- rbind(res, cbind(f, mvar, rownames(coef), coef, baseline))
        dim(res)
      }
      print("3. Test Model")
      #-- create table --#
      res <- as.data.frame(res)
      #colnames(res) <- c("taxa", "metadata", "coefficient", "Estimate", "SE", "df", "t", "pvalue", "baseline")
      print(5)
      
      if(!is.null(StudyID)){
        reg <- "LMEM"
        res$pvalue <- as.numeric(as.character(res$`Pr(>|t|)`))
        
      }else{
        #reg <- "GLM"
        
        tt <- try(res$pvalue <- as.numeric(as.character(res$`Pr(>|z|)`)),T)
        
        if(class(tt) == "try-error"){
          res$pvalue <- as.numeric(as.character(res$`Pr(>|t|)`))
          res$`Pr(>|t|)` <- NULL
        }else{
          res$pvalue <- as.numeric(as.character(res$`Pr(>|z|)`))
          res$`Pr(>|z|)` <- NULL
        }
      }
      
      res$Estimate <- as.numeric(as.character(res$Estimate))
      res$SE <- as.numeric(as.character(res$`Std. Error`))
      res$padj <- p.adjust(res$pvalue, method="fdr")
      res$Model <- reg
      res <- res[order(res$pvalue),]
      
      
      # Make sure tax_table is in the correct format for the merge function
      tax_table <- as(tax_table(psIN.agg), "matrix")
      rownames(res) <- res$f
      res$f <- NULL
      common_otus <- intersect(rownames(res), rownames(tax_table))
      
      # Subset the aldex_diff_abundance and tax_table based on common_otus
      res <- res[rownames(res) %in% common_otus, ]; head(res)
      tax_table <- tax_table[rownames(tax_table) %in% common_otus, ]
      res <- merge(res, tax_table, by.x = "row.names", by.y = "row.names", all.x = TRUE)
      
      res.sel <- res
      res.sel <- as.data.frame(subset(res, pvalue < pval))
      taxa_sig <- res.sel$taxa[1:dim(res.sel)[1]]; summary(taxa_sig)
      
      if(dim(res.sel)[1] == 0){
        next
      }else{
        res.sel$bas.count <-  unique(sum(with(mapping.na, mapping.na[,mvar] == baseline)))
        res.sel$coef.count <-  unique(sum(with(mapping.na, mapping.na[,mvar] == smvar)))
      }
      
      if(dim(res.sel)[1] == 0){
        ps.taxa.sig <- psIN.cb
      }else{
        tt <- try(ps.taxa.sig <- prune_taxa(taxa_sig, psIN.cb),T)
        
        if(class(tt) == "try-error"){
          pathwayTab <- data.frame(otu_table(psIN.cb))
          pathwayRank <- data.frame(tax_table(psIN.cb))
          rownames(pathwayRank) <- pathwayRank[,taxRanks]
          rownames(pathwayTab) <- pathwayRank[,taxRanks]
          pathwayRank <- as.matrix(pathwayRank)
          pathwayTab <- as.matrix(t(pathwayTab))
          psIN.cb <- phyloseq(otu_table(pathwayTab, taxa_are_rows=FALSE), tax_table(pathwayRank));psIN.cb
          ps.taxa.sig <- prune_taxa(taxa_sig, psIN.cb)
          print(ps.taxa.sig)
        }else{
          ps.taxa.sig <- prune_taxa(taxa_sig, psIN.cb)
          print(ps.taxa.sig)
        }
      }
      
      
      res.sel.sel <- subset(res.sel, Genus != "NA")
      
      write.csv(res.sel.sel, quote = FALSE,col.names = NA,file=sprintf("%s/Continuous.out.Sig%s.%s.%s%s.%s%s.csv",out_lmem.Tab,
                                                                       dim(res.sel)[1],
                                                                       mvar, 
                                                                       ifelse(is.null(taxanames), "", paste(taxanames, ".", sep = "")), 
                                                                       project,
                                                                       ifelse(is.null(name), "", paste(name, ".", sep = "")), 
                                                                       reg, 
                                                                       sep="/"))
      saveRDS(ps.taxa.sig,sprintf("%s/Continuous.out.Sig%s.%s.%s%s.%s%s.rds",out_lmem.ps,
                                  dim(res.sel)[1],
                                  mvar, 
                                  ifelse(is.null(taxanames), "", paste(taxanames, ".", sep = "")), 
                                  project,
                                  ifelse(is.null(name), "", paste(name, ".", sep = "")), 
                                  reg))
      
  }
 }
}

