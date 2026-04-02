#' @keywords internal
.gotools_signature <- function() {
  if (requireNamespace("crayon", quietly = TRUE)) {
    cat(crayon::blue("#--------------------------------------------------------------# \n"))
    cat(crayon::blue("#------       General analysis Of microbiome (Go)        ------# \n"))
    cat(crayon::blue("#------    Quick statistics and visualization tools      ------# \n"))
    cat(crayon::blue("#------                (with R markdown)                 ------# \n"))
    cat(crayon::blue("#--------------------------------------------------------------# \n"))
  } else {
    cat("#--------------------------------------------------------------# \n")
    cat("#------       General analysis Of microbiome (Go)        ------# \n")
    cat("#------    Quick statistics and visualization tools      ------# \n")
    cat("#------                (with R markdown)                 ------# \n")
    cat("#--------------------------------------------------------------# \n")
  }
}

#' @keywords internal
.gotools_setup_needed <- function() {
  core_pkgs <- c("phyloseq", "microbiome", "DESeq2", "ALDEx2", "ANCOMBC")
  any(!vapply(core_pkgs, requireNamespace, logical(1), quietly = TRUE))
}

#' @keywords internal
.onAttach <- function(libname, pkgname) {
  .gotools_signature()
  if (.gotools_setup_needed()) {
    packageStartupMessage(
      "Gotools first-time setup:\n",
      "  1. devtools::install_github(\"bbagy/Gotools\", force = TRUE)\n",
      "  2. library(Gotools)\n",
      "  3. Gotool_dependency()"
    )
  }
}
