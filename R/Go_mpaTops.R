
#' Process Kraken and metaphlan MPA Output for Top Taxa
#'
#' This function processes Kraken and metaphlan  MPA output, focusing on top taxa.
#'
#' @param project Name of the project or analysis.
#' @param mpa File path of the Kraken and metaphlan MPA output.
#'
#' @return The function processes the Kraken MPA data, extracts the top taxa, and saves the relevant
#'         information for further RPKM analysis. It saves a Phyloseq object and other required data as .csv and .rds files.
#'
#' @details
#' The function performs a series of steps to process Kraken and metaphlan MPA output. It reads the MPA file,
#' extracts top taxa information, formats the data, and saves it for downstream RPKM analysis.
#' It generates a Phyloseq object and exports relevant tables for species names and OTU.
#'
#' @examples
#' # Example usage:
#' Go_krakenTops(project = "MyProject",
#'               mpa = "path/to/kraken_mpa_output.txt")
#'
#' @export


Go_mpaTops <- function(project,
                       mpa,
                       kingdom = "d__Bacteria") {

  print(sprintf("Currently kingdom is %s", kingdom))
  rds <- file.path(sprintf("%s", "2_rds"))
  if (!file_test("-d", rds)) dir.create(rds)

  # read kraken mpa
  mpatable <- read.table(mpa, header=T, as.is=T, sep="\t", comment.char="", row.names=1, quote="")

  # 특정 킹덤만 필터링
  L1 <- subset(mpatable, grepl(kingdom, rownames(mpatable)))
  L2 <- subset(L1, grepl("p__", rownames(L1)))
  L3 <- subset(L2, grepl("c__", rownames(L2)))
  L4 <- subset(L3, grepl("o__", rownames(L3)))
  L5 <- subset(L4, grepl("f__", rownames(L4)))
  L6 <- subset(L5, grepl("g__", rownames(L5)))
  L7 <- subset(L6, grepl("s__", rownames(L6)))  # Species까지 필터링

  # `t__`가 있는 경우 확인
  has_t__ <- grepl("t__", rownames(L7))

  if (any(has_t__)) {
    # `t__`가 하나라도 있으면 기존 방식 유지
    s_names <- gsub(";t__.*", "", rownames(L7[has_t__, , drop=FALSE]))
    needs_t__ <- !has_t__ & rownames(L7) %in% s_names
    rownames(L7)[needs_t__] <- paste0(rownames(L7)[needs_t__], ";t__", gsub("^.*;s__", "", rownames(L7)[needs_t__]))
  } else {
    # `t__`가 전혀 없는 경우: `s__`를 `t__`로 변환
    rownames(L7) <- paste0(rownames(L7), ";t__", gsub("^.*;s__", "", rownames(L7)))
  }

  # `t__`가 있는 데이터만 유지
  L8 <- subset(L7, grepl("t__", rownames(L7)))

  data <- L8[, setdiff(1:ncol(L8), grep(".Bacterial.kraken.1", colnames(L8)))]

  # 🔹 샘플 이름 정리 코드 다시 추가됨 🔹
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

  # 오류 발생 시 대체 데이터(L8) 활용
  if (inherits(tt, "try-error")) {
    if (any(grepl("_mpa|_out.txt|_metaphlan_bugs_list", colnames(L8)))) {
      colnames(L8) <- gsub("_mpa", "", colnames(L8))
      colnames(L8) <- gsub("_out.txt", "", colnames(L8))
      colnames(L8) <- gsub("_metaphlan_bugs_list", "", colnames(L8))
    }
    if (any(grepl("^X", colnames(L8)))) {
      colnames(L8) <- gsub("^X", "", colnames(L8))
    }
    rownames(L8) <- gsub("\\|", ";", rownames(L8))
    data <- L8
  } else {
    rownames(data) <- gsub("\\|", ";", rownames(data))
  }

  # Taxonomic Rank 정리
  ranklist <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Strain")
  taxlist <- lapply(rownames(data), function(x) {
    parsed <- parse_taxonomy_qiime(x)

    # ✅ `t__`가 존재하면 Species 다음에 Strain으로 저장
    if ("Rank8" %in% names(parsed)) {
      names(parsed)[which(names(parsed) == "Rank8")] <- "Strain"
    }

    return(parsed)
  })

  # 데이터 프레임 변환
  taxa <- as.data.frame(matrix(NA, nrow=nrow(data), ncol=length(ranklist)))
  colnames(taxa) <- ranklist

  for (i in seq_along(taxlist)) {
    if (!is.null(taxlist[[i]])) {
      valid_names <- intersect(names(taxlist[[i]]), ranklist)
      if (length(valid_names) > 0) {
        taxa[i, valid_names] <- taxlist[[i]][valid_names]
      }
    }
  }

  # ✅ Rank1을 항상 "Bacteria"로 설정
  taxa$Kingdom <- "Bacteria"

  # ✅ Strain 설정 (올바른 형식 유지)
  for (i in seq_len(nrow(taxa))) {
    # `t__`가 존재하는 경우 → "Species t__"
    if (!is.na(taxa[i, "Strain"])) {
      taxa[i, "Strain"] <- paste(taxa[i, "Species"], taxa[i, "Strain"])
    }
    # `t__`가 없는 경우 → "Species t__na"
    else {
      taxa[i, "Strain"] <- paste(taxa[i, "Species"], "t__na")
    }
  }

  # phyloseq을 위해 데이터 프레임을 matrix로 변환
  taxa_matrix <- as.matrix(taxa)
  tt <- tax_table(taxa_matrix)
  rownames(tt) <- rownames(data)

  # phyloseq 객체 생성
  ps <- phyloseq(otu_table(data, taxa_are_rows=T), tt)

  print(ps)
  return(ps)
}
