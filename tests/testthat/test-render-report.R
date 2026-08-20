## End-to-end smoke tests: run each Download_templates pipeline against a
## small synthetic dataset and confirm Go_renderReport() produces an HTML
## file. These exercise the same code paths (Go_filter, Go_bdivPM,
## Go_ConDaDist, Go_Groupheatmap, Go_renderReport) that broke during manual
## verification, so they should catch regressions in that pipeline.

skip_if_not_installed("ConDAdist")
skip_if_not_installed("phyloseq")
skip_if_not_installed("ape")
skip_on_cran()

suppressPackageStartupMessages(library(ConDAdist))

run_da_and_report <- function(comparison_ps, project, comparison_var, comparison_orders, dir_root, sumInfo_extra = NULL) {
  da_result_dir <- NULL
  for (method in c("ancombc2", "aldex2", "deseq2")) {
    da_result_dir <- Go_ConDaDist(
      psIN = comparison_ps, project = project, group_var = comparison_var,
      group_1 = comparison_orders[1], group_2 = comparison_orders[2],
      covariates = NULL, methods = method, distances = NULL, name = NULL
    )
  }
  Go_volcanoPlot(project = project, result = da_result_dir, fc = 0, mycols = c("#026CCBFF", "#F51E02FF"), name = NULL, font = 3, height = 5, width = 6)
  extract_da_plot_labels(file.path(dir_root, "pdf", "DA_plot"))
}

test_that("16S template renders end to end", {
  root <- file.path(tempdir(), paste0("gt16S_", as.integer(Sys.time())))
  dir.create(root, recursive = TRUE)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  make_16S_fixture(root)

  old_wd <- getwd()
  on.exit(setwd(old_wd), add = TRUE)
  setwd(root)

  project <- "TestFixture16S"
  today <- format(Sys.Date(), "%y%m%d")
  final_asv_file <- sprintf("1_out/%s.final_asvTable.%s.csv", project, today)
  tree <- ape::read.tree("1_out/TREE/tree.nwk")
  sampledata <- read.csv("3_map/MAPPING.csv", row.names = 1, check.names = FALSE)

  ps0 <- phyloseq::merge_phyloseq(Go_tabTops(csv = final_asv_file, project = project),
                                   phyloseq::sample_data(sampledata), phyloseq::phy_tree(tree))
  ps_filtered <- Go_filter(Go_SeqLengths(ps0, from = 397, to = 450), cutoff = 0.00005)

  comparison_var <- "TreatmentGroup"
  comparison_orders <- c("GroupA", "GroupB")
  ps_main <- ps_filtered

  Go_barchart(psIN = ps_main, project = project, cutoff = 0.005, simple = FALSE, relative = TRUE,
              taxanames = "Phylum", cate.vars = comparison_var, facet = NULL, orders = comparison_orders,
              x_label = NULL, mycols = Go_myCols(custumCols = "cols3"), legend = "bottom",
              name = ensure_report_token("Phylum"), height = 4, width = 8)

  adiv_all <- Go_adiv(psIN = ps_main, project = project, alpha_metrics = c("Chao1", "Shannon"), name = ensure_report_token("AllSamples"))
  Go_boxplot(df = adiv_all, project = project, cate.vars = comparison_var, outcomes = c("Chao1", "Shannon"),
             orders = comparison_orders, mycols = Go_myCols(piratepal = "basel"), covariates = NULL, paired = NULL,
             name = ensure_report_token(comparison_var), cutoff = 0.1, plotCols = 2, plotRows = 1)

  Go_bdivPM(psIN = ps_main, cate.vars = comparison_var, project = project, orders = comparison_orders,
            distance_metrics = c("bray"), mycols = Go_myCols(piratepal = "basel"), name = ensure_report_token(comparison_var),
            facet = NULL, strata_var = NULL, height = 4.5, width = 7, plotCols = 1, plotRows = 1)

  dir_root <- sprintf("%s_%s", project, today)
  da_labels <- run_da_and_report(ps_main, project, comparison_var, comparison_orders, dir_root)

  alpha_tab_file <- sprintf("%s/table/adiv/adiv.%s.%s.%s.csv", dir_root, project, ensure_report_token("AllSamples"), today)

  expInfo <- Go_expInfo(Project_name = "Test", Samples_info = "n=10", Sequencing_date = today,
                         Sequencing_platform = "test", kit_number = 6, pos_number = 1, prep_number = 1, spikein_number = 2,
                         authorName1 = "Test Author", authorEmail1 = "test@example.edu")
  dirInfo <- Go_dirInfo(Current_working_dir = getwd(), Image_locations = sprintf("%s/pdf/", dir_root), DA_image_location = sprintf("%s/pdf/DA_plot/", dir_root))
  tabInfo <- Go_tabInfo(ASVs_Tab = final_asv_file, Tract_Tab = "1_out/track.csv", Alpha_divTab = alpha_tab_file)
  imgInfo <- Go_imgInfo(Rarefaction = NA, Barchart = "Phylum", Bac.heatmap = comparison_var, Adivplot = comparison_var, Bdivplot = comparison_var, DAplot = da_labels)
  sumInfo <- Go_sumInfo(Prep_overview = "test")

  out <- Go_renderReport(family = "16S", project = project, expInfo = expInfo, dirInfo = dirInfo,
                          tabInfo = tabInfo, imgInfo = imgInfo, sumInfo = sumInfo, output_dir = file.path(root, "report"))
  expect_true(file.exists(out))
  expect_gt(file.size(out), 1e5)
  expect_length(
    list.files(dirInfo$image.wd, pattern = "^ordi\\.PCoA\\.PM\\.bray.*\\.png$"),
    1
  )
})

test_that("MG template renders end to end", {
  root <- file.path(tempdir(), paste0("gtMG_", as.integer(Sys.time())))
  dir.create(root, recursive = TRUE)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  make_MG_fixture(root)

  old_wd <- getwd()
  on.exit(setwd(old_wd), add = TRUE)
  setwd(root)

  project <- "TestFixtureMG"
  today <- format(Sys.Date(), "%y%m%d")
  sampledata <- read.csv("3_map/MAPPING.csv", row.names = 1, check.names = FALSE)

  ps0 <- phyloseq::merge_phyloseq(Go_kraken2Tops(kraken2 = "1_out/KRAKEN2_MPA.txt", bracken = "1_out/BRACKEN_MPA.txt", kingdom = "Bacteria"),
                                   phyloseq::sample_data(sampledata))
  ps_main <- Go_filter(ps0, cutoff = 0.00001)

  comparison_var <- "TreatmentGroup"
  comparison_orders <- c("GroupA", "GroupB")

  Go_barchart(psIN = ps_main, project = project, cutoff = 0.005, simple = FALSE, relative = TRUE,
              taxanames = "phylum", cate.vars = comparison_var, facet = NULL, legend = "bottom",
              orders = comparison_orders, x_label = NULL, mycols = Go_myCols(custumCols = "cols3"),
              name = "(Phylum)", height = 4, width = 8)

  adiv_all <- Go_adiv(psIN = ps_main, project = project, alpha_metrics = c("Chao1", "Shannon"), name = ensure_report_token("AllSamples"))
  Go_boxplot(df = adiv_all, project = project, cate.vars = comparison_var, outcomes = c("Chao1", "Shannon"),
             orders = comparison_orders, mycols = Go_myCols(piratepal = "basel"), covariates = NULL, paired = NULL,
             name = ensure_report_token(comparison_var), cutoff = 0.1, plotCols = 2, plotRows = 1)

  Go_bdivPM(psIN = ps_main, cate.vars = comparison_var, project = project, orders = comparison_orders,
            distance_metrics = c("jaccard"), mycols = Go_myCols(piratepal = "basel"), name = ensure_report_token(comparison_var),
            facet = NULL, strata_var = NULL, height = 4.5, width = 7, plotCols = 1, plotRows = 1)

  dir_root <- sprintf("%s_%s", project, today)
  da_labels <- run_da_and_report(ps_main, project, comparison_var, comparison_orders, dir_root)
  alpha_tab_file <- sprintf("%s/table/adiv/adiv.%s.%s.%s.csv", dir_root, project, ensure_report_token("AllSamples"), today)

  # HUMAnN3 functional branch
  humann_pathabundance_file <- "1_out/PATHABUNDANCE.tsv"
  humann_genefamilies_file  <- "1_out/GENEFAMILIES.tsv"
  ps2.path <- Go_function2ps(tabPath = humann_pathabundance_file, project = project, func.type = "Humann", name = NULL)
  ps2.path <- phyloseq::merge_phyloseq(ps2.path, phyloseq::sample_data(data.frame(sampledata)))
  Go_extendedBarplot(psIN = ps2.path, project = project, mvar = comparison_var, group1 = comparison_orders[1], group2 = comparison_orders[2], func = "path.des", wilcox.p = 0.1, name = "PATHWAY", height = 3, width = 10)
  Go_psTotab(ps2.path, project = "PATH")
  path_df <- read.csv(sprintf("1_out/PATH.%s.psTotab.asvTable.csv", today), row.names = 1, check.names = FALSE)
  path_df$pathway <- NULL
  path_df$path.des <- NULL
  Go_Groupheatmap(df = path_df, SampleData = sampledata, project = project, Group = comparison_var, orders = comparison_orders, title = NULL, mycols = NULL, name = ensure_report_token("PATHWAY"), top_n = 50, normalization = "log", width = 11, height = 8, x_label = NULL)
  functional_eb_tokens <- sprintf("%s.vs.%s.PATHWAY", comparison_orders[1], comparison_orders[2])
  functional_sh_tokens <- "PATHWAY"

  expInfo <- Go_expInfo(Project_name = "Test", Samples_info = "n=10", Sequencing_date = today,
                         Sequencing_platform = "test", kit_number = 2, pos_number = 1, prep_number = 6, spikein_number = 2,
                         authorName1 = "Test Author", authorEmail1 = "test@example.edu")
  dirInfo <- Go_dirInfo(Current_working_dir = getwd(), Image_locations = sprintf("%s/pdf/", dir_root), DA_image_location = sprintf("%s/pdf/DA_plot/", dir_root))
  tabInfo <- Go_tabInfo(Taxa_Tab = "1_out/BRACKEN_MPA.txt", Func_Tab = humann_pathabundance_file, Tract_Tab = "1_out/QC_summary.csv", Alpha_divTab = alpha_tab_file)
  imgInfo <- Go_imgInfo(Barchart = "phylum", Bac.heatmap = comparison_var, Adivplot = comparison_var, Bdivplot = comparison_var, DAplot = da_labels,
                         EBplot = functional_eb_tokens, SHeatmap = functional_sh_tokens)
  sumInfo <- Go_sumInfo(Prep_overview = "test")

  out <- Go_renderReport(family = "MG", project = project, expInfo = expInfo, dirInfo = dirInfo,
                          tabInfo = tabInfo, imgInfo = imgInfo, sumInfo = sumInfo, output_dir = file.path(root, "report"))
  expect_true(file.exists(out))
  expect_gt(file.size(out), 1e5)
  expect_length(list.files(dirInfo$image.wd, pattern = "^barchart.*\\(Phylum\\).*\\.png$"), 1)
  expect_length(
    list.files(dirInfo$image.wd, pattern = "^ordi\\.PCoA\\.PM\\.jaccard.*\\.png$"),
    1
  )
})

test_that("RNAseq template renders end to end", {
  root <- file.path(tempdir(), paste0("gtRNAseq_", as.integer(Sys.time())))
  dir.create(root, recursive = TRUE)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  make_RNAseq_fixture(root)

  old_wd <- getwd()
  on.exit(setwd(old_wd), add = TRUE)
  setwd(root)

  project <- "TestFixtureRNAseq"
  today <- format(Sys.Date(), "%y%m%d")
  sampledata <- read.csv("3_map/MAPPING.csv", row.names = 1, check.names = FALSE)
  ps.rna <- phyloseq::merge_phyloseq(Go_tabTops(csv = "1_out/RNAseq_expression_table.csv", project = project),
                                      phyloseq::sample_data(sampledata))

  comparison_var <- "TreatmentGroup"
  comparison_orders <- c("GroupA", "GroupB")

  Go_pheatmap(psIN = ps.rna, project = project, Ntax = 20, title = NULL, group1 = comparison_var,
              show_rownames = TRUE, show_colnames = FALSE, showPhylum = FALSE,
              cluster_rows = TRUE, cluster_cols = TRUE, cutree_cols = NA, cutree_rows = NA,
              name = ensure_report_token("taxa heatmap"), width = 8)

  dir_root <- sprintf("%s_%s", project, today)
  da_labels <- run_da_and_report(ps.rna, project, comparison_var, comparison_orders, dir_root)

  expInfo <- Go_expInfo(Project_name = "Test", Samples_info = "n=10", Sequencing_date = today,
                         Sequencing_platform = "test", RNAseq_reference = "test-ref",
                         kit_number = 1, pos_number = 1, prep_number = 5, spikein_number = 2,
                         authorName1 = "Test Author", authorEmail1 = "test@example.edu")
  dirInfo <- Go_dirInfo(Current_working_dir = getwd(), Image_locations = sprintf("%s/pdf/", dir_root), DA_image_location = sprintf("%s/pdf/DA_plot/", dir_root))
  tabInfo <- Go_tabInfo(RNAseq = "1_out/RNAseq_expression_table.csv", Tract_Tab = "1_out/QC_summary.csv")
  imgInfo <- Go_imgInfo(RNAseq.heatmap = c("taxa heatmap"), DAplot = da_labels)
  sumInfo <- Go_sumInfo(Prep_overview = "test")

  out <- Go_renderReport(family = "RNAseq", project = project, expInfo = expInfo, dirInfo = dirInfo,
                          tabInfo = tabInfo, imgInfo = imgInfo, sumInfo = sumInfo, output_dir = file.path(root, "report"))
  expect_true(file.exists(out))
  expect_gt(file.size(out), 1e5)
  expect_length(list.files(dirInfo$image.wd, pattern = "^pheatmap.*\\(taxa heatmap\\).*\\.png$"), 1)
})
