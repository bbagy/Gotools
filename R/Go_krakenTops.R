
#' Process Kraken MPA Output for Top Taxa and Prepare for RPKM Analysis
#'
#' This function processes Kraken MPA output, focusing on top taxa, and prepares data for RPKM analysis.
#'
#' @param project Name of the project or analysis.
#' @param mpa File path of the Kraken MPA output.
#'
#' @return The function processes the Kraken MPA data, extracts the top taxa, and saves the relevant
#'         information for further RPKM analysis. It saves a Phyloseq object and other required data as .csv and .rds files.
#'
#' @details
#' The function performs a series of steps to process Kraken MPA output. It reads the MPA file,
#' extracts top taxa information, formats the data, and saves it for downstream RPKM analysis.
#' It generates a Phyloseq object and exports relevant tables for species names and OTU.
#'
#' @examples
#' # Example usage:
#' Go_krakenTops(project = "MyProject",
#'               mpa = "path/to/kraken_mpa_output.txt")
#'
#' @export


Go_krakenTops <- function(project,
                          mpa,
                          kingdom = "d__Bacteria"){

  # out dir
  # RPKM_tem <- file.path(sprintf("%s","RPKM_tem"));if(!file_test("-d", RPKM_tem)) dir.create(RPKM_tem)
  # RPKM_input <- file.path(sprintf("%s","RPKM_input"));if(!file_test("-d", RPKM_input)) dir.create(RPKM_input)

  print(sprintf("Currently kingdom is %s",kingdom))
  rds <- file.path(sprintf("%s","2_rds"));if(!file_test("-d", rds)) dir.create(rds)

  # read kraken mpa
  mpatable <- read.table(mpa, header=T, as.is=T, sep="\t", comment.char="", row.names=1, quote="");head(mpatable)

  #L1 <- subset(mpatable, grepl("k__Fungi", rownames(mpatable)))
  L1 <- subset(mpatable, grepl(kingdom, rownames(mpatable)))
  L2 <- subset(L1, grepl("p__", rownames(L1)))
  L3 <- subset(L2, grepl("c__", rownames(L2)))
  L4 <- subset(L3, grepl("o__", rownames(L3)))
  L5 <- subset(L4, grepl("f__", rownames(L4)))
  L6 <- subset(L5, grepl("g__", rownames(L5)))
  L7 <- subset(L6, grepl("s__", rownames(L6)))

  data <- L7[, setdiff(1:ncol(L7), grep(".Bacterial.kraken.1", colnames(L7)))]

  # 데이터 정리 시도
  tt <- try({
    if (any(grepl("_mpa|_out.txt|_metaphlan_bugs_list", colnames(data)))) {
      colnames(data) <- gsub("_mpa", "", colnames(data))
      colnames(data) <- gsub("_out.txt", "", colnames(data))
      colnames(data) <- gsub("_metaphlan_bugs_list", "", colnames(data))
    }

    if (any(grepl("^X", colnames(data)))) {
      colnames(data) <- gsub("^X", "", colnames(data))
    }
  }, silent = TRUE)

  # 오류 발생 시 대체 데이터(L7) 활용
  if (inherits(tt, "try-error")) {
    if (any(grepl("_mpa|_out.txt|_metaphlan_bugs_list", colnames(L7)))) {
      colnames(L7) <- gsub("_mpa", "", colnames(L7))
      colnames(L7) <- gsub("_out.txt", "", colnames(L7))
      colnames(L7) <- gsub("_metaphlan_bugs_list", "", colnames(L7))
    }

    if (any(grepl("^X", colnames(L7)))) {
      colnames(L7) <- gsub("^X", "", colnames(L7))
    }

    rownames(L7) <- gsub("\\|", ";", rownames(L7))
    data <- L7
  } else {
    rownames(data) <- gsub("\\|", ";", rownames(data))
  }

  # 최종 컬럼명 확인
  head(colnames(data))



  #ranklist <- c("Rank1", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  ranklist <- c("Rank1", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  taxlist <- lapply(rownames(data), function(x) parse_taxonomy_qiime(x))

  # 데이터 프레임 대신 matrix를 사용하면 인덱싱 오류가 발생할 수 있으므로 data.frame으로 변환
  taxa <- as.data.frame(matrix(NA, nrow=nrow(data), ncol=length(ranklist)))
  colnames(taxa) <- ranklist

  for (i in seq_along(taxlist)) {
    if (!is.null(taxlist[[i]])) {  # NULL 값 방지
      valid_names <- intersect(names(taxlist[[i]]), ranklist)  # ranklist에 있는 이름만 선택
      if (length(valid_names) > 0) {
        taxa[i, valid_names] <- taxlist[[i]][valid_names]
      }
    }
  }

  # phyloseq을 위해 데이터 프레임을 matrix로 변환
  taxa_matrix <- as.matrix(taxa)

  # phyloseq 객체 생성
  tt <- tax_table(taxa_matrix)
  rownames(tt) <- rownames(data)
  ps <- phyloseq(otu_table(data, taxa_are_rows=T), tt)

  print(ps)
  return(ps)


  #random_tree = rtree(ntaxa(ps), rooted=TRUE, tip.label=taxa_names(ps))
  #ps5 <- merge_phyloseq(ps, random_tree);ps5
  #ps.relative <- transform_sample_counts(ps, function(x) x / sum(x) )
  #ps.rarefied <- rarefy_even_depth(ps, sample.size=40224, rngseed=nsamples(ps))


  #------------------------------#
  #---       Kraken RPKM      ---#
  #------------------------------#
  # make a list species names
  # taxaTab <- data.frame(taxa)
  # taxaTab$Species <-  gsub(" ", "_", taxaTab$Species);head(taxaTab$Species)

  # write.table(taxaTab$Species, quote=F, sep="\t", row.names=F, col.names=F,
  #             file=sprintf("%s/%s.speciesName.txt", RPKM_input, project))

  # write.csv(data.frame(taxa), quote = FALSE,col.names = NA, #row.names = FALSE,
  #           file=sprintf("%s/%s.taxaTable_for_rpkm.%s.csv", RPKM_input, project, format(Sys.Date(), "%y%m%d"),sep="/"))

  # run GoGenomeSize.sh speciesName.txt genomeSize.txt  / 이름을 Gofriend로 할까나.

  # make a list species names
  # otu <- otu_table(data, taxa_are_rows=T)
  # rownames(otu) <- taxaTab$Species
  # otu <- data.frame(otu)

  # saveRDS(ps, sprintf("%s/ps.%s.%s.rds",rds, project, format(Sys.Date(), "%y%m%d")))

  # write.csv(otu, quote = FALSE,col.names = NA, #row.names = FALSE,
  #           file=sprintf("%s/%s.speciesTab_for_rpkm.%s.csv", RPKM_input, project, format(Sys.Date(), "%y%m%d"),sep="/"))


}
