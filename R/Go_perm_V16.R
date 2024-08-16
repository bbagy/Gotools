
#' Global Permutation Analysis for Microbiome Data
#'
#' This function performs global permutation (PERMANOVA) analysis on microbiome data.
#'
#' @param psIN A phyloseq object containing microbiome data.
#' @param vars A vector of categorical variables for which PERMANOVA will be performed.
#' @param project A string representing the project name, used for file naming and directory creation.
#' @param distance_metrics A vector of distance metrics to be used in the analysis.
#' @param mul.vars A boolean value indicating whether to use multiple variables in the PERMANOVA model.
#' @param strata_var A variable used to define strata for the PERMANOVA analysis, ensuring that permutations are constrained within each stratum (e.g., patient ID).
#' @param name An optional string for naming the output files.
#'
#' @details
#' The function conducts PERMANOVA analysis to understand the impact of various categorical variables on microbial community composition. It supports multiple distance metrics and can handle both single and multiple variables in the model. The function outputs the results in a CSV file and returns them as a data frame.
#'
#' @return
#' Returns a data frame containing the results of the PERMANOVA analysis, including degrees of freedom, sums of squares, R-squared values, F-model statistics, p-values, and adjusted p-values.
#'
#' @examples
#' # Example usage:
#' permanova_results <- Go_perm(psIN = my_phyloseq_object,
#'                              vars = c("TreatmentGroup", "AgeGroup"),
#'                              project = "MyMicrobiomeStudy",
#'                              distance_metrics = c("bray", "unifrac"),
#'                              multi = FALSE,
#'                              strata_var=NULL,
#'                              name = "MyAnalysis")
#'
#' @export

Go_perm <- function(psIN, vars, project, distance_metrics, multi=FALSE, name=NULL, strata_var=NULL){
  # out dir
  out <- file.path(sprintf("%s_%s",project, format(Sys.Date(), "%y%m%d")))
  if(!file_test("-d", out)) dir.create(out)
  out_path <- file.path(sprintf("%s_%s/table",project, format(Sys.Date(), "%y%m%d")))
  if(!file_test("-d", out_path)) dir.create(out_path)
  out_dist <- file.path(sprintf("%s_%s/table/dist",project, format(Sys.Date(), "%y%m%d")))
  if(!file_test("-d", out_dist)) dir.create(out_dist)

  results_df <- data.frame()
  for (dist in distance_metrics) {
    # 파일 경로 설정
    file_path <- sprintf("%s/%s_distance.RData", out_dist,dist)

    # 거리 행렬 계산 또는 로드
    if (file.exists(file_path)) {
      load(file_path)
      print(sprintf("%s loading the distance_matrix", dist))
    } else {
      print(sprintf("%s calculating the distance_matrix", dist))
      distance_matrix <- phyloseq::distance(psIN, method = dist)
      save(distance_matrix, file = file_path)
    }

    # 거리 행렬과 map의 샘플 수 및 순서 일치 확인
    distance_samples <- rownames(as.matrix(distance_matrix))

    map <- data.frame(sample_data(psIN))
    map_samples <- rownames(map)

    # map을 distance_matrix의 순서에 맞추어 정렬
    map <- map[distance_samples, ]

    # NA 값이 있는 행 제거
    complete_cases <- complete.cases(map[, vars])
    map_no_na <- map[complete_cases, ]
    distance_matrix_no_na <- as.matrix(distance_matrix)[complete_cases, complete_cases]

    # 데이터 확인
    if (!all(rownames(distance_matrix_no_na) == rownames(map_no_na))) {
      stop("The sample order of the distance matrix and the data frame does not match.")
    }

    if (multi) {
      # 다변량 PERMANOVA 분석
      formula_str <- paste("distance_matrix_no_na ~", paste(vars, collapse = " + "))
      form <- as.formula(formula_str)

      set.seed(123)
      print(paste("Running multivariate adonis2 with", dist, "distance"))
      if (!is.null(strata_var)) {
        permanova_result <- adonis2(form, data = map_no_na, permutations = 999, strata = map_no_na[[strata_var]])
      } else {
        permanova_result <- adonis2(form, data = map_no_na, permutations = 999)
      }

      # PERMANOVA 결과를 데이터 프레임으로 변환하고 변수와 메소드 정보를 추가
      result_df <- as.data.frame(permanova_result)
      result_df$Variable <- "Multivariate"
      result_df$Method <- dist

      # 결과를 결과 데이터 프레임에 추가
      results_df <- rbind(results_df, result_df)
    } else {
      # 각 변수에 대해 단변량 PERMANOVA 분석
      for (var in vars) {
        formula_str <- paste("distance_matrix_no_na ~", var)
        form <- as.formula(formula_str)

        set.seed(123)
        print(paste("Running univariate adonis2 for", var, "with", dist, "distance"))
        if (!is.null(strata_var)) {
          permanova_result <- adonis2(form, data = map_no_na, permutations = 999, strata = map_no_na[[strata_var]])
        } else {
          permanova_result <- adonis2(form, data = map_no_na, permutations = 999)
        }

        # PERMANOVA 결과를 데이터 프레임으로 변환하고 변수와 메소드 정보를 추가
        result_df <- as.data.frame(permanova_result)
        result_df$Variable <- var
        result_df$Method <- dist

        # 결과를 결과 데이터 프레임에 추가
        results_df <- rbind(results_df, result_df)
      }
    }
  }

  # 결과 출력
  write.csv(results_df, sprintf("%s/permanova.%s.%s.%s%s.csv",
                                out_dist,
                                project,
                                ifelse(multi, "multivariate", "univariate"),
                                ifelse(is.null(strata_var), "", paste("strata=", strata_var, ".", sep = "")),
                                ifelse(is.null(name), "", paste(name, ".", sep = "")),
                                format(Sys.Date(), "%y%m%d")))
}

