
#' Convert Functional Data to Phyloseq Object
#'
#' Converts functional profiling data from PICRUSt2 or HUMAnN into a \code{phyloseq} object.
#' Automatically detects the functional type (KEGG KO, Pathway, EC) from row names.
#'
#' @param tabPath Path to the functional data file (tab-delimited, row names as feature IDs).
#' @param project (Optional) Project name included in the output RDS file name.
#' @param func.type Data source: \code{"PICRUSt"} (or \code{"PICRUSt2"}) or \code{"HUMAnN"} (or \code{"HUMAnN2"}, \code{"HUMAnN3"}). Case-insensitive.
#' @param name (Optional) Additional suffix for the output RDS file name.
#' @param collapse_pairs Logical. For HUMAnN inputs only. When \code{TRUE}, automatically detects
#'   and collapses paired technical read columns (e.g., \code{_R1}/\code{_R2}, \code{_r1}/\code{_r2},
#'   \code{_1}/\code{_2}) into a single sample-level column by summing. The pattern is detected
#'   from column names — no manual suffix specification required. Default is \code{FALSE}.
#'
#' @details
#' For PICRUSt2 input, feature IDs are taken from row names and classified as KEGG KO
#' (\code{^K\\d{5}}), Pathway (\code{^PWY}), or EC number (\code{^\\d+\\.\\d+\\.\\d+\\.\\d+}).
#'
#' For HUMAnN input, rows containing \code{UNGROUPED}, \code{UNMAPPED}, \code{UNINTEGRATED},
#' \code{unclassified}, or \code{ribosomal protein} are removed. Feature IDs and descriptions
#' are split on \code{": "}.
#'
#' The output RDS is saved to \code{2_rds/} in the working directory with the naming convention:
#' \code{ps.<func.type>.<func>.<project>.<name>.<YYMMDD>.rds}.
#'
#' @return A \code{phyloseq} object with \code{otu_table} and \code{tax_table}. Also saves an RDS file.
#'
#' @examples
#' # PICRUSt2
#' Go_function2ps(tabPath = "path/to/pred_metagenome_unstrat.tsv",
#'                project = "MyProject",
#'                func.type = "PICRUSt2")
#'
#' # HUMAnN with automatic R1/R2 collapsing
#' Go_function2ps(tabPath = "path/to/humann_genefamilies.tsv",
#'                project = "MyProject",
#'                func.type = "HUMAnN",
#'                collapse_pairs = TRUE)
#'
#' @export

Go_function2ps <- function(tabPath, project = NULL, func.type, name = NULL, collapse_pairs = FALSE) {
  # Read tab
  func.tab <- try(read.delim(tabPath, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE), silent = TRUE)

  if (inherits(func.tab, "try-error")) {
    print("Error reading the table. Check the file path and format.")
    return(NULL)
  }

  # Identify data type and perform necessary preprocessing
  func.type <- tolower(func.type)

  if (func.type %in% c("picrust", "picrust2", "picrust2")) {
    func.type <- "PICRUSt2"

    # ID (KO, PWY, EC 등) 설정
    if (!is.null(rownames(func.tab))) {
      func.tab$ID <- rownames(func.tab)
    } else if ("function" %in% colnames(func.tab)) {
      func.tab$ID <- func.tab$`function`
    } else {
      stop("Cannot identify function ID column.")
    }

    # Sample matrix와 tax table 분리
    NumOfSample <- ncol(func.tab) - 1
    otu <- as.matrix(func.tab[, 1:NumOfSample])
    tax <- as.matrix(func.tab[, NumOfSample + 1, drop = FALSE])
    rownames(otu) <- func.tab$ID
    rownames(tax) <- func.tab$ID

    # ID가 KO, Pathway, 또는 EC인지 판단
    if (any(grepl("^K\\d{5}$", rownames(tax)))) {
      func <- "KEGG"
      colnames(tax) <- "KO"

    } else if (any(grepl("^PWY", rownames(tax)))) {
      func <- "pathway"
      colnames(tax) <- "Pathway"

    } else if (any(grepl("^\\d+\\.\\d+\\.\\d+\\.\\d+$", rownames(tax)))) {
      func <- "EC"
      colnames(tax) <- "EC"

    } else {
      func <- "Unknown"
      colnames(tax) <- "Function"
    }

  } else if (func.type %in% c("humann", "humann2", "humann3")) {
    func.type <- "Humann3"

    collapse_humann_pairs <- function(df) {
      sample_cols <- setdiff(colnames(df), c("ID", "Description"))

      # R1 컬럼 탐색: _R1, _r1, _1 패턴 (우선순위 순)
      pair_patterns <- list(
        list(r1 = "_(R1)(_|$)", r2_replace = "_R2"),
        list(r1 = "_(r1)(_|$)", r2_replace = "_r2"),
        list(r1 = "_(1)(_|$)",  r2_replace = "_2")
      )

      pairs   <- list()
      paired_cols <- character(0)

      for (pat in pair_patterns) {
        r1_candidates <- grep(pat$r1, sample_cols, value = TRUE, perl = TRUE)
        for (col in r1_candidates) {
          if (col %in% paired_cols) next
          partner <- sub(pat$r1, paste0(pat$r2_replace, "\\2"), col, perl = TRUE)
          if (partner %in% sample_cols && !partner %in% paired_cols) {
            base <- sub(pat$r1, "\\2", col, perl = TRUE)  # suffix 이후 부분 보존
            base <- sub("^_", "", base)                    # 앞 _ 제거
            # base가 비어있으면 col에서 _R1 부분만 제거
            if (nchar(base) == 0) base <- sub(pat$r1, "", col, perl = TRUE)
            pairs[[base]] <- c(col, partner)
            paired_cols <- c(paired_cols, col, partner)
          }
        }
      }

      # 페어 없는 컬럼은 그대로 유지
      for (col in setdiff(sample_cols, paired_cols)) {
        pairs[[col]] <- col
      }

      n_pairs <- sum(lengths(pairs) == 2)
      if (n_pairs > 0) {
        message(sprintf("Go_function2ps(): auto-detected %d R1/R2 paired column(s). Collapsing by sum.", n_pairs))
      } else {
        warning("Go_function2ps(): collapse_pairs = TRUE but no R1/R2 pairs detected in column names.")
      }

      out <- lapply(pairs, function(cols) {
        vals <- as.matrix(df[, cols, drop = FALSE])
        mode(vals) <- "numeric"
        rowSums(vals, na.rm = TRUE)
      })

      collapsed <- as.data.frame(out, check.names = FALSE, stringsAsFactors = FALSE)
      rownames(collapsed) <- rownames(df)
      collapsed
    }

    collapse_humann_ids <- function(df) {
      sample_cols <- setdiff(colnames(df), c("ID", "Description"))

      if (!anyDuplicated(df$ID)) {
        return(df)
      }

      message(sprintf(
        "Go_function2ps(): detected %d duplicated HUMAnN feature ID(s). Collapsing duplicated rows by sum.",
        length(unique(df$ID[duplicated(df$ID)]))
      ))

      summed <- stats::aggregate(
        df[, sample_cols, drop = FALSE],
        by = list(ID = df$ID),
        FUN = function(x) sum(as.numeric(x), na.rm = TRUE)
      )

      descriptions <- stats::aggregate(
        Description ~ ID,
        data = df[, c("ID", "Description"), drop = FALSE],
        FUN = function(x) {
          x <- x[!is.na(x) & x != ""]
          if (length(x) == 0) "Unknown" else x[1]
        }
      )

      merged <- merge(descriptions, summed, by = "ID", sort = FALSE)
      merged <- merged[, c("ID", "Description", sample_cols), drop = FALSE]
      rownames(merged) <- merged$ID
      merged
    }

    # 특정 단어 포함된 행 제거
    func.tab <- func.tab[!grepl("ribosomal protein|UNGROUPED|unclassified|UNMAPPED|UNINTEGRATED", rownames(func.tab)), ]

    # ID와 Description 분리
    # HUMAnN pathway tables use "ID: Description", while KO tables may provide only IDs.
    split_ids <- stringr::str_split_fixed(rownames(func.tab), ": ", n = 2)
    func.tab1 <- data.frame(
      ID = split_ids[, 1],
      Description = split_ids[, 2],
      stringsAsFactors = FALSE
    )
    func.tab1$Description[func.tab1$Description == "" | is.na(func.tab1$Description)] <- "Unknown"

    # 병합 및 재구성
    func.tab <- cbind(func.tab1, func.tab)
    func.tab <- collapse_humann_ids(func.tab)
    rownames(func.tab) <- func.tab$ID

    if (isTRUE(collapse_pairs)) {
      collapsed_otu <- collapse_humann_pairs(func.tab)
      otu <- as.matrix(collapsed_otu)
    } else {
      NumOfSample <- ncol(func.tab)
      otu <- as.matrix(func.tab[, 3:NumOfSample])
    }
    tax <- as.matrix(func.tab[, 1:2])
    rownames(otu) <- tax[, 1]
    rownames(tax) <- tax[, 1]

    # KO, PWY, EC 판별
    if (any(grepl("^K\\d{5}", tax[, 1]))) {
      colnames(tax) <- c("KO", "KO.des")
      func <- "KEGG"

    } else if (any(grepl("^PWY", tax[, 1]))) {
      colnames(tax) <- c("pathway", "path.des")
      func <- "pathway"

    } else if (any(grepl("^\\d+\\.\\d+\\.\\d+\\.\\d+$", tax[, 1]))) {
      colnames(tax) <- c("EC", "EC.des")
      func <- "EC"

    } else {
      colnames(tax) <- c("Function", "Description")
      func <- "Unknown"
    }

  } else {
    stop("Please define the correct data type: PICRUSt2 or Humann3")
  }

  # Merge phyloseq
  ps <- phyloseq::phyloseq(
    phyloseq::otu_table(otu, taxa_are_rows = TRUE),
    phyloseq::tax_table(tax)
  )

  print(ps)

  return(ps)
}
