
#' Generate QQ Plots and Histograms for Alpha Diversity Metrics
#'
#' @param psIN Phyloseq object containing microbiome data.
#' @param project Name of the project, used for output file naming.
#' @param alpha_metrics Vector of alpha diversity metrics to analyze (e.g., "Shannon", "Simpson").
#' @param name Optional name to include in the output file's name.
#' @param height Height of the output PDF file.
#' @param width Width of the output PDF file.
#'
#' @details
#' This function creates QQ plots and histograms for various alpha diversity metrics specified in the `alpha_metrics` argument. It computes these metrics using the `estimate_richness` function from the phyloseq package and then generates corresponding plots to assess the distribution and normality of each metric.
#'
#' @return
#' Generates a PDF file containing QQ plots and histograms for each alpha diversity metric. The file is saved in the specified output directory. Additionally, the function prints the Shapiro-Wilk normality test results for each metric to the console.
#'
#' @examples
#' Go_qq(psIN = my_phyloseq_object,
#'       project = "MyMicrobiomeProject",
#'       alpha_metrics = c("Shannon", "Simpson"),
#'       name = "Analysis1",
#'       height = 10,
#'       width = 8)
#'
#' @export

Go_qq <- function(psIN, project, alpha_metrics, name, height, width){
    if(!is.null(dev.list())) dev.off()
  # out dir
  out <- file.path(sprintf("%s_%s",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out)) dir.create(out)
  out_path <- file.path(sprintf("%s_%s/pdf",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out_path)) dir.create(out_path)
  
  
      # logic for out file
  pdf(sprintf("%s/QQplot.%s%s.%s.pdf", out_path, 
              project, 
              ifelse(is.null(name), "", paste(name, ".", sep = "")), 
              format(Sys.Date(), "%y%m%d")), height = height, width = width)
  
  # 1st adiv table
  mapping.sel <- data.frame(sample_data(psIN))
  adiv <- estimate_richness(psIN, measures=alpha_metrics)
  rownames(adiv) <- gsub("^X", "", rownames(adiv))
  adiv$SampleID <- rownames(adiv)
  rownames(adiv) <- rownames(mapping.sel)
  adiv <- merge(adiv, mapping.sel, by="row.names")
  rownames(adiv) <- adiv$SampleID
  adiv$ShannonLn <-log(adiv$Shannon)
  # show last column name
  rev(names(adiv))[1]
  
  #----------- QQ plot and histogram -----------#
  par(mfrow = c(3,2))
  mes <- c(alpha_metrics, rev(names(adiv))[1])
  for (am in mes){
    test <- shapiro.test(adiv[,am])
    hist(adiv[,am], freq=F, xlab= am, main=sprintf("Histogram of %s (%s)", project, am ), cex.main=1) 
    lines(density(adiv[,am])) 
    rug(adiv[,am])
    # remove inf
    adiv.inf <- adiv[!is.infinite(adiv[,am]),]
    
    qqnorm(adiv.inf[,am], main=sprintf("Normal Q-Q Plot (%s p=%.2g)", "shapiro", test$p.value), cex.main=1)
    qqline(adiv.inf[,am])
    
    print(sprintf("%s %s shapiro test (p=%.2g)",am, project, test$p.value))
  }
  
  dev.off()
}
