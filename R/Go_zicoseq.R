#' Go_zicoseq
#'
#' This function performs a differential abundance test using ZicoSeq on microbiome count data provided in a Phyloseq object. It filters data based on a categorical grouping and adjusts for covariates if provided.
#'
#' @param psIN A `phyloseq` object containing species count data along with sample and taxonomic metadata.
#' @param cate.outs A character string specifying the column in the sample metadata that defines the categorical groups for comparison.
#' @param group1 The first group in the comparison.
#' @param group2 The second group in the comparison.
#' @param rand.eff An optional character vector naming random effects variables to be used in the model, default is NULL.
#' @param con.vari An optional character string specifying a column in the sample metadata that should be used as a covariate in the model, default is NULL.
#' @param orders A character vector specifying the order of levels for the factor used in the model. This is particularly useful for setting the reference level in the analysis.
#'
#' @return A list containing the results of the ZicoSeq analysis, including p-values and effect sizes for features that are differentially abundant between the specified groups.
#'
#' @details
#' The function subsets the `phyloseq` object to include only the samples of interest based on `cate.outs`. It then relevels the factor of interest to set the first specified group as the reference. It filters the OTU table to remove features with zero counts across all samples, and optionally adjusts for covariates and random effects if specified.
#'
#' @importFrom phyloseq sample_data otu_table tax_table
#' @importFrom dplyr filter
#' @import GUniFrac
#'
#' @examples
#' \dontrun{
#' data(GlobalPatterns)
#' psIN <- GlobalPatterns
#' result <- Go_zicoseq(psIN, "SampleType", "Feces", "Soil", orders = c("Feces", "Soil"))
#' print(result)
#' }
#' @export
#'
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

  non_zero_var_rows <- apply(otu_mat, 1, function(x) var(x) != 0)
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



