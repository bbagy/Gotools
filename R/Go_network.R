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
#' @import base read.csv setNames
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
  library(igraph)
  library(ggplot2)
  # 패키지 설치 (설치되어 있지 않은 경우)
  if (!requireNamespace("purrr", quietly = TRUE)) install.packages("purrr")

  library(purrr)
  library(dplyr)
  library(reshape2)
  
  # 데이터 처리 및 전치(transpose) 함수 정의
  process_table <- function(tab) {
    # 수치형 컬럼만 선택하여 rowSums 계산 후 새로운 컬럼 추가
    tab$RowSum <- rowSums(tab[, sapply(tab, is.numeric), drop = FALSE])
    
    # Species와 RowSum을 조합하여 새로운 이름 생성
    if ("Species" %in% colnames(tab)) {
      tab$names <- paste(tab$Species, tab$RowSum, sep = "_")
    } else {
      tab$names <- paste("Unknown", tab$RowSum, sep = "_")  # Species가 없을 경우 대비
    }
    
    # rownames을 위에서 생성한 새로운 이름으로 설정
    rownames(tab) <- tab$names
    
    # rownames에서 '['와 ']' 제거
    rownames(tab) <- gsub("\\[", "", gsub("\\]", "", rownames(tab)))
    
    # 특정 taxonomic rank 컬럼 제거 (존재하는 경우에만 삭제)
    for (rank in c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "RowSum", "names")) {
      if (rank %in% colnames(tab)) {
        tab[, rank] <- NULL
      }
    }
    
    # 데이터 전치 (transpose) 후 반환
    return(t(tab))
  }
  
  # 데이터 로드 함수
  read_and_process <- function(file_path) {
    if (!is.null(file_path) && file.exists(file_path)) {
      tab <- read.csv(file_path, row.names = 1, check.names = FALSE)
      return(process_table(tab))  # 전처리 및 전치 포함
    } else {
      return(NULL)
    }
  }
  
  # tab1은 반드시 존재해야 함
  tab1 <- read_and_process(tab1_path)
  sampledata <- read.csv(Sampledata, row.names = 1, check.names = FALSE)
  
  if (is.null(tab1)) stop("tab1 파일이 존재하지 않습니다.")
  
  # tab2와 tab3는 선택적으로 로드
  tab2 <- read_and_process(tab2_path)
  tab3 <- read_and_process(tab3_path)
  
  # 공통 샘플 찾기
  common_samples <- rownames(tab1)
  if (!is.null(tab2)) common_samples <- intersect(common_samples, rownames(tab2))
  if (!is.null(tab3)) common_samples <- intersect(common_samples, rownames(tab3))
  
  # sampledata 존재 여부 확인
  if (exists("sampledata") && !is.null(sampledata)) {
    common_samples <- intersect(common_samples, rownames(sampledata))
    sampledata <- sampledata[common_samples, , drop = FALSE]
  } else {
    stop("sampledata가 존재하지 않습니다.")
  }
  
  # 공통 샘플 유지
  tab1 <- tab1[common_samples, , drop = FALSE]
  if (!is.null(tab2)) tab2 <- tab2[common_samples, , drop = FALSE]
  if (!is.null(tab3)) tab3 <- tab3[common_samples, , drop = FALSE]
  
  # 병합 수행 (tab1은 필수, tab2, tab3는 있는 경우만 추가)
  merged_table.1 <- tab1
  if (!is.null(tab2)) merged_table.1 <- merge(merged_table.1, tab2, by = "row.names", all = TRUE)
  if (!is.null(tab3)) merged_table.1 <- merge(merged_table.1, tab3, by = "Row.names", all = TRUE)
  
  # sampledata 병합
  sampledata$Row.names <- rownames(sampledata)
  #final_merged_table <- merge(merged_table.1, sampledata[, c("Row.names", mainGroup)], by = "Row.names", all = TRUE)
  
  
  
  # `merged_table.1`의 rownames를 "Row.names" 컬럼으로 변환 (병합을 위해 필요)
  if (!"Row.names" %in% colnames(merged_table.1)) {
    merged_table.1 <- data.frame(Row.names = rownames(merged_table.1), merged_table.1, check.names = FALSE)
  }
  
  # `sampledata`가 존재하는 경우에만 병합 진행
  if (exists("sampledata") && !is.null(sampledata)) {
    
    # `sampledata`에 "Row.names" 컬럼이 없으면 rownames로 생성
    if (!"Row.names" %in% colnames(sampledata)) {
      sampledata <- data.frame(Row.names = rownames(sampledata), sampledata, check.names = FALSE)
    }
    
    # 병합 수행 (두 데이터프레임 모두 "Row.names"이 존재함)
    final_merged_table <- merge(merged_table.1, sampledata[, c("Row.names", mainGroup)], 
                                by = "Row.names", all = TRUE)
    
  } else {
    warning("⚠️ sampledata가 존재하지 않거나 'Row.names' 컬럼이 없습니다. 병합을 건너뜁니다.")
    final_merged_table <- merged_table.1  # 병합하지 않고 원본 유지
  }
  
  # 병합 후 rownames 복원
  rownames(final_merged_table) <- final_merged_table$Row.names
  final_merged_table$Row.names <- NULL  # 불필요한 컬럼 제거
  
  
  
  
  
  
  
  
  
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
  
  final_merged_table <- merge(merged_table.1, sampledata[, c("Row.names", mainGroup)], by = "Row.names", all = TRUE)
  
  
  
  
  
  # 데이터 필터링
  data <- final_merged_table %>%
    filter(!!sym(mainGroup) == subgroup) %>%  # mainGroup 변수 사용하여 동적으로 필터링
    select(-all_of(mainGroup), -Row.names)
  
  # 문자열을 숫자로 변환
  data <- data %>%
    mutate(across(everything(), ~ as.numeric(trimws(.))))
  
  # (2) 중복된 `.x`, `.y` 컬럼 정리
  data_clean <- data %>%
    select(-matches("\\.y$"))  # ".y"로 끝나는 중복 컬럼 제거
  colnames(data_clean) <- gsub("\\.x$", "", colnames(data_clean))  # ".x" 삭제하여 컬럼명 정리
  
  
  
  
  
  #===== Run
  #=== Step 2: 상관관계 및 p-value 계산 ===#
  cor_results <- expand.grid(Source = colnames(data), Target = colnames(data)) %>%
    filter(Source != Target) %>%
    rowwise() %>%
    mutate(
      test_result = list(cor.test(data[[Source]], data[[Target]], method = "kendall", use = "pairwise.complete.obs")),
      Correlation = test_result$estimate,
      p_value = test_result$p.value
    ) %>%
    ungroup() %>%
    mutate(q_value = p.adjust(p_value, method = "fdr")) %>%  # FDR 보정
    select(Source, Target, Correlation, p_value, q_value)
  
  #=== Step 3: edge 데이터 생성 (p-value 적용) ===#
  
  if (!is.null(qval_threshold)) {
    # FDR (q-value) 기준 필터링
    edges <- cor_results %>%
      filter(abs(Correlation) > cutoff & q_value < qval_threshold) %>%
      select(Source, Target, Correlation, p_value, q_value)
    sig <- "FDR"
    sigval <- "FDR < 0.05"
    print(sig)
    
  } else if (!is.null(pval_threshold)) {
    # p-value 기준 필터링
    edges <- cor_results %>%
      filter(abs(Correlation) > cutoff & p_value < pval_threshold) %>%
      select(Source, Target, Correlation, p_value, q_value)
    sig <- "p"
    sigval <- "p < 0.05"
    print(sig)
    
  } else {
    # 필터링 없이 모든 관계 포함
    edges <- cor_results %>%
      filter(abs(Correlation) > cutoff) 
    sig <- "p_FDR"
    signame <-"FDR"
    sigval <- "FDR Highlighted"
    print(sig)
  }
  
  # 중복 엣지 제거
  edges_unique <- edges %>%
    mutate(pair = pmap_chr(list(Source, Target), ~ paste(sort(c(.x, .y)), collapse = "-"))) %>%
    group_by(pair) %>%
    #slice(1) %>%
    ungroup() %>%
    select(Source, Target, Correlation, p_value, q_value)  # p_value 포함
  
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
    filter(Bacteria %in% unique(c(edges_unique$Source, edges_unique$Target)))
  
  
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
  
  # 노드 타입별 색상 매핑 (유연하게 변경 가능)
  color_mapping <- setNames(
    c("deepskyblue", "yellow", "green", "gray"),  # 색상 목록
    c(tab1_name, tab2_name, tab3_name, "Unknown") # 실제 타입명
  )
  
  # 노드 색상 할당 (유연한 설정)
  node_colors <- sapply(V(network)$Type, function(type) {
    color_mapping[[type]] %||% "gray"  # 해당 타입이 없으면 기본값(gray) 사용
  })
  
  
  # out dir
  out <- file.path(sprintf("%s_%s",project, format(Sys.Date(), "%y%m%d")))
  if(!file_test("-d", out)) dir.create(out)
  out_pdf <- file.path(sprintf("%s_%s/pdf",project, format(Sys.Date(), "%y%m%d")))
  if(!file_test("-d", out_pdf)) dir.create(out_pdf)
  out_table <- file.path(sprintf("%s_%s/table",project, format(Sys.Date(), "%y%m%d")))
  if(!file_test("-d", out_table)) dir.create(out_table)
  out_network <- file.path(sprintf("%s_%s/table/network",project, format(Sys.Date(), "%y%m%d")))
  if(!file_test("-d", out_network)) dir.create(out_network)
  
  
  # 네트워크 시각화
  if (!is.null(dev.list())) dev.off()
  
  pdf(sprintf("%s/network.%s.%s.(%s).%s.%s.%s%s.pdf",  out_pdf, mainGroup,subgroup, cutoff,sig,
              ifelse(node_names, "IDs", "Names"),  
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
  
  
  used_types <- unique(V(network)$Type)
  legend_labels <- used_types[used_types %in% names(color_mapping)]
  legend_colors <- sapply(legend_labels, function(type) color_mapping[[type]] %||% "gray")
  
  
  # 범례 (legend) 동적 생성
  legend("bottomright", 
         legend = legend_labels, 
         col = legend_colors,
         pch = 19,
         pt.cex = 1.5,
         bty = "n",
         title = "Node Types")
  
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
  
  if (sig == "All"){
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


