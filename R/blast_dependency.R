.gotools_blast_min_version <- "2.8.0"
.gotools_blast_bundle_version <- "2.17.0"

.gotools_blast_root <- function() {
  file.path(tools::R_user_dir("Gotools", which = "data"), "blast")
}

.gotools_blast_version <- function(blastn) {
  if (!nzchar(blastn) || !file.exists(blastn)) return(NA_character_)
  output <- tryCatch(
    system2(blastn, "-version", stdout = TRUE, stderr = TRUE),
    error = function(e) character(0)
  )
  match <- regexpr("[0-9]+\\.[0-9]+\\.[0-9]+", paste(output, collapse = " "))
  if (match[1] < 0) return(NA_character_)
  regmatches(paste(output, collapse = " "), match)[1]
}

.gotools_blast_is_compatible <- function(blastn, minimum = .gotools_blast_min_version) {
  version <- .gotools_blast_version(blastn)
  !is.na(version) && utils::compareVersion(version, minimum) >= 0
}

.gotools_blast_candidates <- function() {
  executable <- if (.Platform$OS.type == "windows") "blastn.exe" else "blastn"
  managed_root <- .gotools_blast_root()
  managed <- if (dir.exists(managed_root)) {
    list.files(managed_root, pattern = paste0("^", executable, "$"), recursive = TRUE, full.names = TRUE)
  } else {
    character(0)
  }
  unique(c(
    unname(Sys.which("blastn")),
    managed,
    file.path("/opt/homebrew/bin", executable),
    file.path("/usr/local/bin", executable)
  ))
}

.gotools_find_blastn <- function(required = TRUE) {
  candidates <- .gotools_blast_candidates()
  candidates <- candidates[nzchar(candidates) & file.exists(candidates)]
  compatible <- candidates[vapply(candidates, .gotools_blast_is_compatible, logical(1))]
  if (length(compatible) > 0) return(normalizePath(compatible[1], mustWork = TRUE))

  if (required) {
    found <- if (length(candidates) > 0) {
      paste0(candidates[1], " (", .gotools_blast_version(candidates[1]), ")")
    } else {
      "none"
    }
    stop(
      "[Gotools] Compatible NCBI BLAST+ was not found. Found: ", found, ".\n",
      "Run Gotool_dependency() to install BLAST+ ", .gotools_blast_min_version, " or newer."
    )
  }
  ""
}

.gotools_blast_archive <- function(version = .gotools_blast_bundle_version) {
  machine <- tolower(Sys.info()[["machine"]])
  platform <- if (.Platform$OS.type == "windows") {
    "x64-win64"
  } else if (Sys.info()[["sysname"]] == "Darwin") {
    if (grepl("arm64|aarch64", machine)) "aarch64-macosx" else "x64-macosx"
  } else if (grepl("arm64|aarch64", machine)) {
    "aarch64-linux"
  } else {
    "x64-linux"
  }
  filename <- sprintf("ncbi-blast-%s+-%s.tar.gz", version, platform)
  list(
    filename = filename,
    url = sprintf("https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/%s/%s", version, filename),
    md5_url = sprintf("https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/%s/%s.md5", version, filename)
  )
}

.gotools_install_blast <- function() {
  archive <- .gotools_blast_archive()
  install_root <- .gotools_blast_root()
  dir.create(install_root, recursive = TRUE, showWarnings = FALSE)
  download_file <- tempfile(fileext = ".tar.gz")
  checksum_file <- tempfile(fileext = ".md5")
  on.exit(unlink(c(download_file, checksum_file)), add = TRUE)

  old_timeout <- getOption("timeout")
  options(timeout = max(3600, old_timeout))
  on.exit(options(timeout = old_timeout), add = TRUE)

  message("[Gotools] Downloading official NCBI BLAST+ ", .gotools_blast_bundle_version, "...")
  utils::download.file(archive$url, download_file, mode = "wb", quiet = FALSE)
  utils::download.file(archive$md5_url, checksum_file, mode = "wb", quiet = TRUE)
  expected_md5 <- tolower(strsplit(readLines(checksum_file, warn = FALSE)[1], "[[:space:]]+")[[1]][1])
  actual_md5 <- unname(tools::md5sum(download_file))
  if (!nzchar(expected_md5) || !identical(tolower(actual_md5), expected_md5)) {
    stop("[Gotools] NCBI BLAST+ download failed its MD5 integrity check.")
  }
  utils::untar(download_file, exdir = install_root)

  blastn <- .gotools_find_blastn(required = FALSE)
  if (!nzchar(blastn)) {
    stop("[Gotools] BLAST+ was downloaded but blastn could not be validated in ", install_root)
  }
  if (.Platform$OS.type != "windows") Sys.chmod(blastn, mode = "0755")
  message("[Gotools] BLAST+ is ready: ", blastn, " (", .gotools_blast_version(blastn), ")")
  blastn
}

.gotools_ensure_blast <- function(ask = interactive()) {
  blastn <- .gotools_find_blastn(required = FALSE)
  if (nzchar(blastn)) {
    message("[Gotools] Compatible BLAST+ found: ", blastn, " (", .gotools_blast_version(blastn), ")")
    return(TRUE)
  }

  if (ask) {
    answer <- readline(
      paste0("[Gotools] Install NCBI BLAST+ ", .gotools_blast_bundle_version,
             " in your Gotools user data folder? [y/N] ")
    )
    if (!tolower(trimws(answer)) %in% c("y", "yes")) {
      message("[Gotools] BLAST+ installation cancelled.")
      return(FALSE)
    }
  }

  .gotools_install_blast()
  TRUE
}
