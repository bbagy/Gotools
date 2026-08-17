test_that("beta report matching supports one, two, and three metrics", {
  root <- tempfile("bdiv-files-")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)

  filenames <- c(
    "ordi.PCoA.PM.bray.Project.Group.(TreatmentGroup).260817.pdf",
    "ordi.PCoA.PM.bray+unifrac.Project.Group.(TreatmentGroup).260817.pdf",
    "ordi.PCoA.PM.bray+unifrac+wunifrac.Project.Group.(TreatmentGroup).260817.pdf"
  )
  file.create(file.path(root, filenames))

  expect_length(Gotools:::.gotools_find_bdiv_pdfs(root, "bray", "TreatmentGroup"), 3)
  expect_length(Gotools:::.gotools_find_bdiv_pdfs(root, "unifrac", "TreatmentGroup"), 2)
  expect_length(Gotools:::.gotools_find_bdiv_pdfs(root, "wunifrac", "TreatmentGroup"), 1)
  expect_length(Gotools:::.gotools_find_bdiv_pdfs(root, "jsd", "TreatmentGroup"), 0)
})

test_that("beta report matching respects exact tokens and subgroups", {
  root <- tempfile("bdiv-token-")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)

  file.create(file.path(root, c(
    "ordi.PCoA.PM.bray.Project.Group.(TreatmentGroup).SiteA.260817.pdf",
    "ordi.PCoA.PM.bray.Project.Group.(TreatmentGroup2).SiteB.260817.pdf"
  )))

  matched <- Gotools:::.gotools_find_bdiv_pdfs(root, "bray", "TreatmentGroup", "SiteA")
  expect_length(matched, 1)
  expect_match(basename(matched), "SiteA", fixed = TRUE)
})
