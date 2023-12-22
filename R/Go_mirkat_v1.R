
#' Microbiome Regression-Based Kernel Association Test (MiRKAT)
#'
#' @param psIN Phyloseq object containing microbiome data.
#' @param project Name of the project, used for output file naming.
#' @param cate.vars List of categorical variables for MiRKAT analysis.
#' @param cate.conf Optional list of confounding variables.
#' @param orders A vector specifying the order of levels in the categorical variable.
#' @param name Optional name for the output file.
#'
#' @details
#' This function performs the MiRKAT (Microbiome Regression-Based Kernel Association Test) on microbiome data. MiRKAT is a non-parametric method used for assessing the association between microbial community composition and a categorical outcome variable, accounting for potential confounders.
#'
#' The function handles multiple comparisons by considering different combinations of groups within the categorical variables and accounts for confounding factors if provided. It computes various UniFrac distances and converts them into kernel matrices for the MiRKAT analysis.
#'
#' @return
#' The function generates a CSV file containing the MiRKAT analysis results, including P-values for the association between microbial communities and the categorical outcome. The file is saved in the specified output directory.
#'
#' @examples
#' Go_mirkat(psIN = my_phyloseq_object,
#'          project = "MyMicrobiomeProject",
#'          cate.vars = c("TreatmentGroup"),
#'          cate.conf = c("Age", "BMI"),
#'          orders = c("Control", "Treatment"),
#'          name = "MiRKAT_analysis")
#'
#' @export

Go_mirkat<- function(psIN, project, cate.vars, cate.conf = NULL,  orders,name=NULL){
  # install bioconductor
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  bioconductors <- c("data.table","phyloseq","dirmult","vegan")
  
  for (bioconductor in bioconductors){
    if(!bioconductor %in% installed.packages()){
      library(BiocManager)
      BiocManager::install(bioconductor)
    }else{library(bioconductor, character.only = TRUE)}
  }
  
  

  # install package from install_github
  github <- c("GLMMMiRKAT")
  if(!github %in% installed.packages()){
    install_github("hk1785/GLMM-MiRKAT", force=T)
  }else{library(package, character.only = TRUE)}
  # install package 
  packages <- c("data.table","CompQuadForm","devtools","ecodist","GUniFrac","GLMMMiRKAT","lme4","MASS","Matrix","MiRKAT","permute") 
  for (package in packages){
    if(!package %in% installed.packages()){
      install.packages(package)
    }else{library(package, character.only = TRUE)}
  }
  
  

  #install.packages("data.table", version = "1.13.0")
  
  
  # out dir
  out <- file.path(sprintf("%s_%s",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out)) dir.create(out)
  out_path <- file.path(sprintf("%s_%s/table",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out_path)) dir.create(out_path)
  out_table <- file.path(sprintf("%s_%s/table/mirkat",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out_table)) dir.create(out_table)
  

  
  # check by variables
  mapping <- data.frame(sample_data(psIN))
  for (mvar in  cate.conf) {
   mapping[,mvar] <- factor(mapping[,mvar])
  }
  
  sample_data(psIN) <- mapping
  
  
  res <- {}
  for (mvar in cate.vars) {
    if (length(unique(mapping[,mvar])) == 1){
      print(sprintf("%s has only 1 variation, which wouldn't be able to compare.",mvar))
      next
    }
    #-----------------------#
    # for group combination #
    #-----------------------#
    # order for each variation
    mapping[,mvar] <- factor(mapping[,mvar], levels = orders)
    mapping[,mvar] <- factor(mapping[,mvar])
    
    # list of the combination
    group.cbn <- combn(x = levels(mapping[,mvar]), m = 2)
    
    #print(count(group.cbn))
    
    group_comparisons <- {}
    for(i in 1:ncol(group.cbn)){
      x <- group.cbn[,i]
      group_comparisons[[i]] <- x
    };group_comparisons
    
    set.seed(1)
    
    # looping for using the combination of list 
    for(i in 1:length(group_comparisons)){
      print(group_comparisons[i])
      group.combination <- unlist(group_comparisons[i]);group.combination
      
      basline <- group.combination[1]
      smvar <- group.combination[2]
      mapping.cbn <- subset(mapping, mapping[,mvar] %in% c(basline,smvar)) 
      
      psIN.cbn <- psIN
      sample_data(psIN.cbn) <- mapping.cbn
      dim(psIN.cbn)
      
      # extract phyloseq elements
      otu.cbn <- data.frame(otu_table(psIN.cbn))
      mapping.cbn <- data.frame(sample_data(psIN.cbn))
      tree.cbn <- phy_tree(psIN.cbn)
      
      
      # create empty df for this combination
      df <- data.frame(matrix(ncol = 1, nrow = dim(mapping.cbn)[1]));df[1] <- list(NULL)
      df.covar <- df
      
      df[,mvar] <- as.numeric(mapping.cbn[,mvar]  == unique(mapping.cbn[,mvar] )[1])
      
      # add corvatiate into the df
      
      if (!is.null(cate.conf)){
      for (covar in cate.conf) {
        df.covar[,covar] <- as.numeric(mapping.cbn[,covar])
        if (mvar == covar){
        next
        }
      }
      }

      
      # Create the UniFrac Distances
      unifracs <- GUniFrac(otu.cbn, tree.cbn, alpha = c(0, 0.5, 1))$unifracs
      D.weighted = unifracs[,,"d_1"]
      D.unweighted = unifracs[,,"d_UW"]
      D.generalized = unifracs[,,"d_0.5"]
      D.BC = as.matrix(vegdist(otu.cbn, method="bray"))
      
      # Convert Distance to Kernel Matrices
      K.weighted = D2K(D.weighted)
      K.unweighted = D2K(D.unweighted)
      K.generalized = D2K(D.generalized)
      K.BC = D2K(D.BC)
      Ks = list(K.weighted = K.weighted, K.unweighted = K.unweighted, K.BC = K.BC)
      
      
      if (length(cate.conf) >=1){
        for (covar in cate.conf) {
          # Cauchy
          # cauchy <- MiRKAT(y = df[,mvar], Ks = Ks, X = df.covar, out_type = "D", method = "davies",  
          # omnibus = "cauchy", returnKRV = FALSE, returnR2 = FALSE)
          # print(cauchy)  tt <- try()
          
          
          tt <- try(permutation <- MiRKAT(y = df[,mvar], Ks = Ks, X = df.covar, out_type = "D", method = "davies"
                                          ,omnibus = "permutation", returnKRV = FALSE, returnR2 = FALSE), T)
          
          if (class(tt) == "try-error"){
            print(sprintf("Fail comparison %s vs %s",basline,smvar))
            
            next
          }
        }
        
      }else if(length(cate.conf) ==0){
        
        permutation <- MiRKAT(y = df[,mvar], Ks = Ks, X = NULL, out_type = "D", method = "davies", 
                              omnibus = "permutation", returnKRV = FALSE, returnR2 = FALSE)
      }
      
      # create table
      per.df <- data.frame(unlist(permutation));per.df
      
      colnames(per.df) <- paste(basline,"vs",smvar);per.df
      per.df.t <- as.data.frame(t(per.df));per.df.t
      
      per.df.t$mvar <- mvar
      class(per.df.t)
      
      
      if (length(cate.conf) >=1){
        covars <- cate.conf[mvar != cate.conf]
        per.df.t$cate.conf <- paste(setdiff(cate.conf, "SampleType"), collapse="+")
        per.df.t$covar <- paste(setdiff(df.covar, "SampleType"), collapse="+")
        
      }else if(length(cate.conf) ==0){
        per.df.t <- per.df.t
      }
      
      
      res <- rbind(res, per.df.t);res
      print(res)
    }
  }
  write.csv(res, quote = FALSE,col.names = NA,file=sprintf("%s/mirkat.%s.%s%s.csv",out_table,
                                                          project,
                                                          ifelse(is.null(name), "", paste(name, ".", sep = "")),  
                                                          format(Sys.Date(), "%y%m%d"), sep="/")) 
  
}
