Go_Maaslin2 <- function(psIN, 
                        project,
                        fixed_effects = NULL,
                        random_effects = NULL,
                        normalization = "TSS",   # MaAsLin2의 normalization 방법 (TSS, CSS, NONE 등)
                        transform = TRUE,        # 상대 abundance 시 median read 곱해 정수 변환
                        orders = NULL,           # 특정 인자를 factor로 지정 + 레벨 순서
                        out_dir = NULL           # 결과 저장 경로 (NULL이면 자동 생성)
) {
  
  ########################################################
  # 1. 출력 디렉토리 생성 (project_YYYYMMDD/table/MaAsLin2)
  ########################################################
  # 예: "MyStudy_230101/table/MaAsLin2"
  out <- file.path(sprintf("%s_%s", project, format(Sys.Date(), "%y%m%d")))
  if (!file_test("-d", out)) dir.create(out)
  
  out_path <- file.path(sprintf("%s_%s/table", project, format(Sys.Date(), "%y%m%d")))
  if (!file_test("-d", out_path)) dir.create(out_path)
  
  out_DA <- file.path(sprintf("%s_%s/table/MaAsLin2", project, format(Sys.Date(), "%y%m%d")))
  if (!file_test("-d", out_DA)) dir.create(out_DA)
  
  # out_dir 인자가 없으면 자동으로 out_DA 사용
  if (is.null(out_dir)) {
    out_dir <- out_DA
  }
  
  message(sprintf("[INFO] MaAsLin2 result will be saved in: %s", out_dir))
  
  ########################################################
  # 2. phyloseq 데이터 추출 + 샘플 이름 확인
  ########################################################
  # 메타데이터
  metadata_df <- as.data.frame(sample_data(psIN))
  
  # OTU 테이블 (행 = feature, 열 = sample)
  # MaAsLin2는 "feature x sample" 이어야 하고,
  # sample_data는 "sample x variables" 형태
  # -> taxa_are_rows(psIN)에 따라 전치 여부
  if (taxa_are_rows(psIN)) {
    otu_mat <- as.matrix(otu_table(psIN))  # 이미 (taxa x sample)
  } else {
    otu_mat <- t(as.matrix(otu_table(psIN)))  # 전치
  }
  
  # 데이터프레임으로 변환
  otu_mat <- as.data.frame(otu_mat)
  
  # 샘플 이름(열 이름)과 메타데이터 rownames 일치 확인
  sample_names_otu <- colnames(otu_mat)
  sample_names_meta <- rownames(metadata_df)
  
  # 교집합만 추출 (불일치 해결)
  common_samples <- intersect(sample_names_otu, sample_names_meta)
  if (length(common_samples) < 2) {
    stop("[ERROR] Not enough common samples between OTU table and metadata!")
  }
  
  # OTU/Metadata를 교집합 샘플만 남긴다
  otu_mat <- otu_mat[, common_samples, drop=FALSE]  # col이 샘플
  metadata_df <- metadata_df[common_samples, , drop=FALSE]
  
  message(sprintf("[INFO] # of common samples: %d", length(common_samples)))
  
  ########################################################
  # 3. Abundance 타입 확인 + 필요 시 변환
  ########################################################
  detect_abundance_type <- function(physeq) {
    lib_sizes <- sample_sums(physeq)
    mean_lib <- mean(lib_sizes)
    if (abs(mean_lib - 100) < 0.1) {
      return("relative")
    } else if (mean_lib > 1000) {
      return("absolute")
    } else {
      return("unknown")
    }
  }
  
  abundance_type <- detect_abundance_type(psIN)
  if (abundance_type == "relative" && transform) {
    # median read를 곱해서 대략적인 정수 형태로 변환
    total_reads <- median(sample_sums(psIN))
    message("[INFO] Detected relative abundance. Multiplying by median read (", total_reads, ") and rounding.")
    otu_mat <- round(otu_mat * total_reads)
  } else {
    message("[INFO] Either 'absolute' data or transform=FALSE. No integer conversion performed.")
  }
  
  ########################################################
  # 4. orders 적용 (factor level 순서)
  ########################################################
  # 예: orders = list(Timepoint = c("D0","D9","D20"))
  if (!is.null(orders) && length(orders) > 0 && !is.null(fixed_effects)) {
    for (fx in fixed_effects) {
      if (fx %in% names(orders)) {
        metadata_df[, fx] <- factor(metadata_df[, fx], levels = orders[[fx]])
      }
    }
  }
  
  ########################################################
  # 5. fixed_effects / random_effects 유효성 확인
  ########################################################
  if (is.null(fixed_effects) || length(fixed_effects) == 0) {
    stop("[ERROR] No fixed effects provided (fixed_effects=NULL). Please specify at least one variable.")
  }
  # fixed_effects가 메타데이터에 실제 존재하는지
  missing_fx <- setdiff(fixed_effects, colnames(metadata_df))
  if (length(missing_fx) > 0) {
    stop("[ERROR] Some fixed_effects not found in metadata: ", paste(missing_fx, collapse=", "))
  }
  
  # random_effects가 메타데이터에 존재하는지
  if (!is.null(random_effects)) {
    missing_re <- setdiff(random_effects, colnames(metadata_df))
    if (length(missing_re) > 0) {
      message("[WARN] Some random_effects not found in metadata: ", paste(missing_re, collapse=", "))
      message("[WARN] They will be ignored.")
      random_effects <- setdiff(random_effects, missing_re)
      if (length(random_effects) == 0) random_effects <- NULL
    }
  }
  
  ########################################################
  # 6. MaAsLin2 실행
  ########################################################
  suppressMessages(library(Maaslin2))
  
  message("[INFO] Running MaAsLin2...")
  
  # otu_mat : feature x sample
  # metadata_df : sample x metadata
  # -> MaAsLin2 자동 감지
  fit_data <- Maaslin2(
    input_data = otu_mat,
    input_metadata = metadata_df,
    output = out_dir,
    fixed_effects = fixed_effects,
    random_effects = random_effects,
    normalization = normalization,
    transform = "LOG",   # 기본 LOG 변환
    plot_heatmap = FALSE,
    plot_scatter = FALSE
  )
  
  # 실행 후, 결과 파일 확인
  files_in_outdir <- list.files(out_dir, recursive = TRUE)
  if (length(files_in_outdir) == 0) {
    message("[ERROR] MaAsLin2 did not generate any output files. Check logs or data format.")
  } else {
    message(sprintf("[SUCCESS] MaAsLin2 completed. Files in: %s", out_dir))
  }
  
  ########################################################
  # 7. (선택) PDF 보고서 생성 (예시)
  ########################################################
  # 원하시면 R Markdown 등을 통해 PDF 자동 생성 가능.
  # 여기서는 간단한 예시 메시지
  # message("[INFO] PDF report not implemented. You can manually create plots from the results TSV.")
  
  return(invisible(fit_data))
}