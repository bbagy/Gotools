#' Go_network
#'
#' Constructs and visualizes a network based on microbial community data across different samples, integrating multiple data tables.
#'
#' @param tab1_path Path to the CSV file containing the first table of data.
#' @param tab2_path Optional path to the CSV file containing the second table of data.
#' @param tab3_path Optional path to the CSV file containing the third table of data.
#' @param Sampledata Path to the CSV file containing sample metadata.
#' @param mainGroup The name of the main grouping variable in the sample metadata.
#' @param subgroup The name of the subgroup for filtering within the main group.
#' @param tab1_name Descriptive name for the first data table used in the network visualization.
#' @param tab2_name Optional descriptive name for the second data table.
#' @param tab3_name Optional descriptive name for the third data table.
#' @param cutoff Correlation coefficient threshold for including edges in the network.
#' @param pval_threshold Optional threshold for p-values to filter edges based on statistical significance.
#' @param qval_threshold Optional threshold for FDR-adjusted q-values to filter edges based on statistical significance.
#' @param node_font Font size for node labels in the network plot.
#' @param node_names Boolean to indicate whether to display node names or IDs.
#' @param name Optional name for additional naming in output files.
#' @param width Width of the output plot in inches.
#' @param height Height of the output plot in inches.
#'
#' @return Saves network plots as PDF files and optionally outputs centrality measures for nodes in the network.
#'
#' @details
#' This function reads microbial community data from specified CSV files, processes it, and constructs a network based on correlations among community members. It allows for extensive customization of the network visualization, including the ability to filter edges based on correlation coefficients, p-values, and q-values. Additionally, it can merge multiple tables and utilize sample metadata for grouping and subsetting the data.
#'
#' @importFrom dplyr filter select mutate left_join
#' @importFrom igraph graph_from_data_frame degree betweenness closeness evcent
#' @importFrom ggplot2 ggplot aes_string geom_point
#' @importFrom stats p.adjust
#' @importFrom utils read.csv
#' @importFrom stats setNames
#'
#' @examples
#' \dontrun{
#' Go_network(
#'   tab1_path = "path/to/data1.csv",
#'   tab2_path = "path/to/data2.csv",
#'   Sampledata = "path/to/sampledata.csv",
#'   mainGroup = "Treatment",
#'   subgroup = "Control",
#'   tab1_name = "Data Table 1",
#'   tab2_name = "Data Table 2",
#'   cutoff = 0.3,
#'   pval_threshold = 0.05,
#'   node_font = 12,
#'   node_names = TRUE,
#'   name = "NetworkPlot",
#'   width = 12,
#'   height = 12
#' )
#' }
#' @export

Go_network <- function(
    tab1_path,
    tab2_path = NULL,
    tab3_path = NULL,
    Sampledata,

    mainGroup,
    subgroup,

    tab1_name,
    tab2_name = NULL,
    tab3_name = NULL,

    #===== Parameters
    cutoff = 0.3,      # 상관계수 절대값 필터링
    pval_threshold = NULL, # p-value 필터링
    qval_threshold = NULL, # FDR 필터링
    node_font,
    node_names = TRUE,
    name = NULL,
    width = 10,
    height = 10
) {

  library(dplyr)
  library(ggplot2)
  # 패키지 설치 (설치되어 있지 않은 경우)
  if (!requireNamespace("purrr", quietly = TRUE)) install.packages("purrr")
  if (!requireNamespace("igraph", quietly = TRUE)) install.packages("igraph")
  if (!requireNamespace("RColorBrewer", quietly = TRUE)) install.packages("RColorBrewer")  ### [추가]

  library(igraph)
  library(purrr)
  library(dplyr)
  library(reshape2)
  library(RColorBrewer)  ### [추가]


  # 데이터 처리 및 전치(transpose) 함수 정의
  process_table <- function(tab) {
    # 수치형 컬럼만 선택하여 rowSums 계산 후 새로운 컬럼 추가
    tab$RowSum <- rowSums(tab[, sapply(tab, is.numeric), drop = FALSE])

    # Species와 RowSum을 조합하여 새로운 이름 생성
    if ("Species" %in% colnames(tab)) {
      tab$names <- paste(tab$Species, tab$RowSum, sep = "_")
      tab$names <- make.unique(tab$names)  # 중복 방지
    } else {
      tab$names <- rownames(tab)  # Species가 없을 경우 기존 rownames 유지
    }

    # rownames 설정 및 특수문자 제거
    rownames(tab) <- gsub("\\[|\\]", "", tab$names)

    # 특정 taxonomic rank 컬럼 제거 (존재하는 경우만 삭제)
    for (rank in c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "RowSum", "names")) {
      if (rank %in% colnames(tab)) {
        tab[, rank] <- NULL
      }
    }

    # 데이터 전치 후 데이터프레임으로 반환
    return(as.data.frame(t(tab), stringsAsFactors = FALSE))
  }

  # 📌 데이터 로드 함수 (파일 존재 여부 및 공백 처리)
  read_and_process <- function(file_path) {
    if (!is.null(file_path) && file.exists(file_path)) {
      tab <- suppressWarnings(
        read.csv(file_path, row.names = 1, check.names = FALSE, fill = TRUE, strip.white = TRUE)
      )

      # 비어있는 경우 NULL 반환
      if (nrow(tab) == 0 || ncol(tab) == 0) {
        warning(sprintf("⚠️ Warning: %s is empty or malformed. Returning NULL.", file_path))
        return(NULL)
      }

      # 데이터 처리 및 전치
      return(process_table(tab))
    } else {
      return(NULL)
    }
  }

  # 📌 tab1 (필수), tab2, tab3 로드
  tab1 <- read_and_process(tab1_path)
  if (is.null(tab1)) stop("❌ tab1 파일이 존재하지 않습니다.")

  tab2 <- read_and_process(tab2_path)
  tab3 <- read_and_process(tab3_path)

  # 📌 sampledata 로드

  if (class(Sampledata) == "character"){
    sampledata <- read.csv(Sampledata, row.names = 1, check.names = FALSE)
  }else{
    sampledata = Sampledata
  }



  # 📌 공통 샘플 찾기 (sampledata까지 포함)
  common_samples <- rownames(tab1)
  if (!is.null(tab2)) common_samples <- intersect(common_samples, rownames(tab2))
  if (!is.null(tab3)) common_samples <- intersect(common_samples, rownames(tab3))
  common_samples <- intersect(common_samples, rownames(sampledata))  # sampledata 포함

  # 📌 공통 샘플 유지
  tab1 <- tab1[common_samples, , drop = FALSE]
  if (!is.null(tab2)) tab2 <- tab2[common_samples, , drop = FALSE]
  if (!is.null(tab3)) tab3 <- tab3[common_samples, , drop = FALSE]
  sampledata <- sampledata[common_samples, , drop = FALSE]  # sampledata도 동일 적용

  # 📌 병합 수행 (tab1은 필수, tab2, tab3는 있는 경우만 추가)
  merged_table.1 <- tab1

  if (!is.null(tab2)) {
    merged_table.1 <- merge(merged_table.1, tab2, by = "row.names", all = TRUE)
  }

  if (!is.null(tab3)) {
    tab3 <- as.data.frame(tab3)
    tab3$Row.names <- rownames(tab3)
    merged_table.1 <- merge(merged_table.1, tab3, by = "Row.names", all = TRUE)
  }

  # 📌 sampledata 병합
  sampledata$Row.names <- rownames(sampledata)

  # `merged_table.1`의 rownames를 "Row.names" 컬럼으로 변환 (병합을 위해 필요)
  if (!"Row.names" %in% colnames(merged_table.1)) {
    merged_table.1 <- data.frame(Row.names = rownames(merged_table.1), merged_table.1, check.names = FALSE)
  }

  # 📌 최종 병합 (sampledata는 필수)
  final_merged_table <- merge(merged_table.1, sampledata[, c("Row.names", mainGroup)], by = "Row.names", all = TRUE)

  # 📌 최종 공통 샘플 유지
  final_merged_table <- final_merged_table[final_merged_table$Row.names %in% common_samples, , drop = FALSE]


  rownames(final_merged_table) <- final_merged_table$Row.names


  # (3) 유일한 세균 리스트 생성
  all_bacteria <- colnames(tab1)  # tab1은 항상 존재

  if (!is.null(tab2)) {
    all_bacteria <- unique(c(all_bacteria, colnames(tab2)))
  }
  if (!is.null(tab3)) {
    all_bacteria <- unique(c(all_bacteria, colnames(tab3)))
  }

  cat("\n🔍 total bacteria after de-duplicate:", length(all_bacteria), "\n")

  # (4) global_bacteria_map 생성 (중복 없이 고유 ID 할당)
  if (!exists("global_bacteria_map")) {  # 처음 한 번만 실행
    global_bacteria_map <- data.frame(
      Bacteria = all_bacteria,
      ID = as.character(seq_along(all_bacteria))  # 1부터 고유한 ID 할당
    )
  }


  # (6) 최종 결과 확인
  # cat("\n✅ Global Map 생성 완료: 총", nrow(global_bacteria_map), "개의 세균이 포함됨 ✅\n")
  # print(global_bacteria_map)
  global_bacteria_map

  bacteria_map <- global_bacteria_map


  # 데이터 필터링
  data <- final_merged_table %>%
    dplyr::filter(!!sym(mainGroup) == subgroup) %>%  # mainGroup 변수 사용하여 동적으로 필터링
    dplyr::select(-all_of(mainGroup), -Row.names)

  # 문자열을 숫자로 변환
  data <- data %>%
    dplyr::mutate(across(everything(), ~ as.numeric(trimws(.))))

  # (2) 중복된 `.x`, `.y` 컬럼 정리
  data_clean <- data %>%
    dplyr::select(-matches("\\.y$"))  # ".y"로 끝나는 중복 컬럼 제거
  colnames(data_clean) <- gsub("\\.x$", "", colnames(data_clean))  # ".x" 삭제하여 컬럼명 정리





  #===== Run
  #=== Step 2: 상관관계 및 p-value 계산 ===#
  cor_results <- expand.grid(Source = colnames(data), Target = colnames(data)) %>%
    dplyr::filter(Source != Target) %>%
    rowwise() %>%
    dplyr::mutate(
      test_result = list(cor.test(data[[Source]], data[[Target]], method = "kendall", use = "pairwise.complete.obs")),
      Correlation = test_result$estimate,
      p_value = test_result$p.value
    ) %>%
    ungroup() %>%
    dplyr::mutate(q_value = p.adjust(p_value, method = "fdr")) %>%  # FDR 보정
    dplyr::select(Source, Target, Correlation, p_value, q_value)

  #=== Step 3: edge 데이터 생성 (p-value 적용) ===#

  if (!is.null(qval_threshold)) {
    # FDR (q-value) 기준 필터링
    edges <- cor_results %>%
      dplyr::filter(abs(Correlation) > cutoff & q_value < qval_threshold) %>%
      dplyr::select(Source, Target, Correlation, p_value, q_value)
    sig <- "FDR"
    sigval <- "FDR < 0.05"
    print(sig)

  } else if (!is.null(pval_threshold)) {
    # p-value 기준 필터링
    edges <- cor_results %>%
      dplyr::filter(abs(Correlation) > cutoff & p_value < pval_threshold) %>%
      dplyr::select(Source, Target, Correlation, p_value, q_value)
    sig <- "p"
    sigval <- "p < 0.05"
    print(sig)

  } else {
    # 필터링 없이 모든 관계 포함
    edges <- cor_results %>%
      dplyr::filter(abs(Correlation) > cutoff)
    sig <- "p_FDR"
    signame <-"FDR"
    sigval <- "FDR Highlighted"
    print(sig)
  }

  # 중복 엣지 제거
  edges_unique <- edges %>%
    dplyr::mutate(pair = pmap_chr(list(Source, Target), ~ paste(sort(c(.x, .y)), collapse = "-"))) %>%
    dplyr::group_by(pair) %>%
    #slice(1) %>%
    dplyr::ungroup() %>%
    dplyr::select(Source, Target, Correlation, p_value, q_value)  # p_value 포함

  #=== Step 4: node 데이터 생성 ===#
  all_nodes <- unique(c(edges_unique$Source, edges_unique$Target))

  # 데이터프레임 생성 (초기 Type을 NA로 설정)
  nodes <- data.frame(Node = all_nodes, Type = NA)

  # Type 지정 (각 노드가 어느 그룹에 속하는지 체크)
  if (!is.null(tab1)) {
    nodes$Type[nodes$Node %in% colnames(tab1)] <- tab1_name
  }

  if (!is.null(tab2)) {
    nodes$Type[nodes$Node %in% colnames(tab2)] <- tab2_name
  }

  if (!is.null(tab3) && !is.null(tab3_name)) {
    nodes$Type[nodes$Node %in% colnames(tab3)] <- tab3_name
  }

  # 그룹에 속하지 않은 노드 처리 (기본값 설정)
  nodes$Type[is.na(nodes$Type)] <- "Unknown"  # 필요에 따라 변경 가능

  #=== Step 5: 네트워크 생성 ===#
  network <- graph_from_data_frame(d = edges_unique, vertices = nodes, directed = FALSE)

  #=== Replacing Node names ad number
  # (1) edges_unique에 등장하는 세균만 선택하여 bacteria_map 생성
  bacteria_map_filtered <- bacteria_map %>%
    dplyr::filter(Bacteria %in% unique(c(edges_unique$Source, edges_unique$Target)))


  # (2) Source & Target을 ID로 변환
  edges_numeric <- edges_unique %>%
    dplyr::left_join(bacteria_map_filtered, by = c("Source" = "Bacteria")) %>%
    dplyr::rename(Source_ID = ID) %>%
    dplyr::left_join(bacteria_map_filtered, by = c("Target" = "Bacteria")) %>%
    dplyr::rename(Target_ID = ID) %>%
    dplyr::select(Source_ID, Target_ID, Correlation, p_value, q_value)  # 숫자로 변환된 노드 사용

  # (3) igraph 객체 생성 (숫자로 변환된 엣지 사용)
  g <- graph_from_data_frame(edges_numeric, directed = FALSE)

  # (4) 노드 레이블을 숫자로 변경
  V(g)$label <- V(g)$name  # 노드 라벨 = 숫자 ID

  #=== Step 6: 시각화 ===#
  loops <- which_loop(network)
  loop_indices <- which(loops)

  # 루프 삭제
  if (length(loop_indices) > 0) {
    network <- delete_edges(network, loop_indices)
    print("Loops removed from the network.")
  } else {
    print("No loops found in the network.")
  }

  #------------------------------------------------------------------------------------#
  # [추가] tab1만 존재할 때 => Genus별 색상 적용 (sub("[_\\.].*", "", ...) 사용)
  #------------------------------------------------------------------------------------#
  if (!is.null(tab1) && is.null(tab2) && is.null(tab3)) {
    # 1) Genus 이름 추출: 첫 번째 '_' 또는 '.' 전까지를 Genus로 인식
    all_bacteria_tab1 <- colnames(tab1)
    genus_names <- sub("[_\\.].*", "", all_bacteria_tab1)  ### [중요 수정]
    genus_names <- sub("[ _.].*", "", all_bacteria_tab1)
    unique_genera <- unique(genus_names)
    num_genera <- length(unique_genera)

    # 2) Genus별로 색상 생성 (최대 9개 이상이면 colorRampPalette로 확장)
    if (num_genera > 9) {
      genus_colors <- setNames(colorRampPalette(brewer.pal(9, "Set1"))(num_genera),
                               unique_genera)
    } else {
      genus_colors <- setNames(brewer.pal(num_genera, "Set1"), unique_genera)
    }

    # 3) 네트워크 노드별 Genus 매핑
    #    (V(network)$name과 tab1의 colnames를 match)
    V(network)$Genus <- genus_names[match(V(network)$name, all_bacteria_tab1)]
    node_colors <- sapply(V(network)$Genus, function(genus) genus_colors[[genus]] %||% "gray")

    # 4) 범례 설정 (Genus Groups)
    legend_labels <- unique_genera
    legend_colors <- genus_colors[legend_labels]
    legend_title <- "Genus Groups"

  } else {
    #----------------------------------------------------------------------------------#
    # [기존 코드] 노드 타입별 색상 매핑 (유연하게 변경 가능)
    #----------------------------------------------------------------------------------#
    color_mapping <- setNames(
      c("deepskyblue", "yellow", "green", "gray"),  # 색상 목록
      c(tab1_name, tab2_name, tab3_name, "Unknown") # 실제 타입명
    )

    # 노드 색상 할당 (유연한 설정)
    node_colors <- sapply(V(network)$Type, function(type) {
      color_mapping[[type]] %||% "gray"  # 해당 타입이 없으면 기본값(gray) 사용
    })

    # legend 설정 (기존 Node Types)
    used_types <- unique(V(network)$Type)
    legend_labels <- used_types[used_types %in% names(color_mapping)]
    if (length(legend_labels) == 0) legend_labels <- "Unknown"
    legend_colors <- sapply(legend_labels, function(type) color_mapping[[type]] %||% "gray")
    legend_title <- "Node Types"
  }
  #------------------------------------------------------------------------------------#
  # [추가/수정 끝]
  #------------------------------------------------------------------------------------#


  # out dir
  out <- file.path(sprintf("%s_%s",project, format(Sys.Date(), "%y%m%d")))
  if(!dir.exists(out)) dir.create(out)
  out_pdf <- file.path(sprintf("%s_%s/pdf",project, format(Sys.Date(), "%y%m%d")))
  if(!dir.exists(out_pdf)) dir.create(out_pdf)
  out_table <- file.path(sprintf("%s_%s/table",project, format(Sys.Date(), "%y%m%d")))
  if(!dir.exists(out_table)) dir.create(out_table)
  out_network <- file.path(sprintf("%s_%s/table/network",project, format(Sys.Date(), "%y%m%d")))
  if(!dir.exists(out_network)) dir.create(out_network)


  # 네트워크 시각화
  if (!is.null(dev.list())) dev.off()

  pdf(sprintf("%s/network.%s.%s.(%s).%s.%s.%s%s%s%s%s.pdf",  out_pdf, mainGroup,subgroup, cutoff,sig,
              ifelse(node_names, "IDs", "Names"),
              ifelse(is.null(tab1_name), "", paste(tab1_name, ".", sep = "")),
              ifelse(is.null(tab2_name), "", paste(tab2_name, ".", sep = "")),
              ifelse(is.null(tab3_name), "", paste(tab3_name, ".", sep = "")),
              ifelse(is.null(name), "", paste(name, ".", sep = "")),
              format(Sys.Date(), "%y%m%d")), width = width, height = height)

  if(sig=="FDR"){
    edge_styles <- ifelse(E(network)$q_value < 0.05, 1, 2)  # 1 = 실선, 2 = 점선
  }else if(sig=="p"){
    edge_styles <- ifelse(E(network)$p_value < 0.05, 1, 2)  # 1 = 실선, 2 = 점선
  }else{
    edge_styles <- ifelse(E(network)$q_value < 0.05, 1, 2)  # 1 = 실선, 2 = 점선
  }


  set.seed(123)
  # (1) 노드 이름을 숫자로 변환할지 여부 설정
  # V(network)$name에 맞춰 ID를 정확히 매핑
  node_labels <- if (node_names) {
    bacteria_map_filtered$ID[match(V(network)$name, bacteria_map_filtered$Bacteria)]
  } else {
    V(network)$name
  }
  plot(network,
       vertex.color = node_colors,
       vertex.label = node_labels,
       vertex.label.cex = node_font,
       vertex.size = 5 + 10 * (degree(network) / max(degree(network))),
       edge.width = abs(E(network)$Correlation) * 5,
       edge.color = ifelse(E(network)$Correlation > 0, "red", "blue"),
       edge.lty = edge_styles,  # 실선(1) vs 점선(2)
       edge.curved = 0.15,
       main = sprintf("%s-%s Group Network (%s)", mainGroup, subgroup,sigval))


  # 범례 (legend) 동적 생성
  legend("bottomright",
         legend = legend_labels,
         col = legend_colors,
         pch = 19,
         pt.cex = 1.5,
         bty = "n",
         title = legend_title)

  legend("bottomleft",
         legend = c("Positive Correlation", "Negative Correlation"),
         col = c("red", "blue"),
         lty = 1,
         lwd = 3,
         bty = "n",
         title = "Edge Correlation")


  if (node_names){
    legend("topleft", legend = paste(bacteria_map$ID, bacteria_map$Bacteria, sep = " = "),
           cex = 0.6, bty = "n")
    print("replacing node names as number")
  }

  if (sig == "p_FDR"){
    legend("topright",
           legend = c(sprintf("%s < 0.05",signame), sprintf("%s ≥ 0.05",signame)),
           lty = c(1, 2),
           lwd = 3,
           bty = "n",
           title = sprintf("%s Significance",signame))
  }

  dev.off()
  #=== Centrality (중심성) 노드 찾기
  degree_values <- degree(network)
  betweenness_values <- betweenness(network)
  closeness_values <- closeness(network)
  eigenvector_values <- evcent(network)$vector

  # 중심성 지표 계산
  centrality_data <- data.frame(
    Node = V(network)$name,  # 노드 이름
    Degree = degree(network),  # Degree 중심성
    Betweenness = betweenness(network),  # Betweenness 중심성
    Closeness = closeness(network),  # Closeness 중심성
    Eigenvector = eigen_centrality(network)$vector  # Eigenvector 중심성
  )

  # 중심성 값 정렬 (Degree 중심성 기준으로 예시)
  sorted_centrality <- centrality_data[order(-centrality_data$Degree), ]

  # 상위 10개 노드 출력
  head(sorted_centrality, 10)
  write.csv(sorted_centrality,sprintf("%s/%s.%s.sorted_centrality.csv",out_network,mainGroup,sig))

}
