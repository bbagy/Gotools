
Go_zicoseq <- function(psIN,
                       cate.outs,
                       group1,
                       group2,
                       rand.eff =NULL,
                       con.vari = NULL,
                       orders){
  #===== Differential abundant test  (ZicoSeq)
  # install.packages("GUniFrac")
  library(GUniFrac)
  
  map <- data.frame(sample_data(psIN))
  
  
  map.cb <- subset(map, map[[cate.outs]] %in% c(group1, group2)) # phyloseq subset은 작동을 안한다.
  
  psIN.cb <- psIN
  sample_data(psIN.cb) <- map.cb
  
  #ps2.sel <- subset_samples(psIN, map[,cate.outs] %in% c(group1,group2));ps2.sel
  #===== extract mata data
  meta_dat <- data.frame(sample_data(psIN.cb))
  meta_dat[,cate.outs] <- as.factor(meta_dat[,cate.outs] )
  
  meta_dat[,cate.outs]  <- relevel(meta_dat[,cate.outs] , ref = group1)
  
  #===== extract otu and species names
  otu_mat <- as.matrix(t(otu_table(psIN.cb)))
  
  species_names <- tax_table(psIN.cb)[, "Species"]  # Assuming "Species" is the column for species names
  rownames(otu_mat) <- species_names
  #===== check the both created data 
  rownames(meta_dat) <- rownames(t(otu_mat))
  
  #===== Option: If the ZicoSeq() show errors related 0 value 
  non_zero_rows <- rowSums(otu_mat) != 0
  
  # Filter out rows with all zeros
  otu_mat_filtered <- otu_mat[non_zero_rows, ]
  
  # Alternatively, you can check for any variance and remove columns with zero variance
  non_zero_var_rows <- apply(otu_mat, 1, var) != 0
  otu_mat_filtered <- otu_mat[non_zero_var_rows, ]
  
  # Check how many features are left after filtering
  cat("Number of features after filtering:", nrow(otu_mat_filtered), "\n")
  set.seed(123)
  
  zico_results <- ZicoSeq(
    meta.dat = meta_dat,
    feature.dat = otu_mat_filtered,
    grp.name = cate.outs,
    adj.name = con.vari,  # Adjust this if you need to account for covariates
    feature.dat.type = 'count',
    prev.filter = 0.1, # 0.2
    mean.abund.filter = 0,
    max.abund.filter = 0.01, # 0.002
    is.winsor = TRUE,
    outlier.pct = 0.05, # 0.03
    is.post.sample = TRUE,
    post.sample.no = 50,#25
    link.func = list(function(x) sign(x) * (abs(x))^0.5),
    perm.no = 999, #99
    strata = rand.eff,
    ref.pct = 0.5,
    stage.no = 6,
    excl.pct = 0.2,
    p.max = 1000,#500
    is.fwer = TRUE, # FALSE
    verbose = TRUE,
    return.feature.dat = TRUE
  )
  
  return(zico_results)
}



