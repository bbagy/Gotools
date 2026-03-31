#' @keywords internal
.gotools_setup_needed <- function() {
  core_pkgs <- c("phyloseq", "microbiome", "DESeq2", "ALDEx2", "ANCOMBC")
  any(!vapply(core_pkgs, requireNamespace, logical(1), quietly = TRUE))
}

#' @keywords internal
.onAttach <- function(libname, pkgname) {
  if (.gotools_setup_needed()) {
    packageStartupMessage(
      "Gotools first-time setup:\n",
      "  1. devtools::install_github(\"bbagy/Gotools\", force = TRUE)\n",
      "  2. library(Gotools)\n",
      "  3. Gotool_dependency()"
    )
  }
}
