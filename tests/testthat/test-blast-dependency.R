test_that("BLAST version parsing and compatibility checks work", {
  skip_on_os("windows")
  fake_dir <- tempfile("blast-test-")
  dir.create(fake_dir)
  on.exit(unlink(fake_dir, recursive = TRUE), add = TRUE)

  old_blast <- file.path(fake_dir, "old-blastn")
  new_blast <- file.path(fake_dir, "new-blastn")
  writeLines(c("#!/bin/sh", "echo 'blastn: 2.4.0+'"), old_blast)
  writeLines(c("#!/bin/sh", "echo 'blastn: 2.17.0+'"), new_blast)
  Sys.chmod(c(old_blast, new_blast), "0755")

  expect_equal(Gotools:::.gotools_blast_version(old_blast), "2.4.0")
  expect_false(Gotools:::.gotools_blast_is_compatible(old_blast))
  expect_true(Gotools:::.gotools_blast_is_compatible(new_blast))
})

test_that("the official archive matches the current operating system", {
  archive <- Gotools:::.gotools_blast_archive("2.17.0")
  expect_match(archive$filename, "^ncbi-blast-2\\.17\\.0\\+-.+\\.tar\\.gz$")
  expect_match(archive$url, "/2\\.17\\.0/ncbi-blast-")
  expect_equal(archive$md5_url, paste0(archive$url, ".md5"))
})
