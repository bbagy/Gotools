
#' Install or repair Gotools dependencies
#'
#' Checks a curated set of CRAN and Bioconductor packages used by `Gotools`.
#' Packages that are already installed and loadable are left untouched.
#' Missing packages are installed, and installed-but-broken packages are
#' reinstalled.
#'
#' @param ask If `TRUE` (default in interactive sessions), prompt before
#'   installing. Set `FALSE` to install without prompting.
#' @param attach If `TRUE` (default), attach a curated set of essential
#'   analysis packages after dependency checks finish successfully.
#'
#' @details
#' This function is designed to support interactive `Gotools` analyses.
#' It checks whether required packages are installed and loadable, installs or
#' repairs missing packages when needed, and can also attach a curated set of
#' essential analysis libraries commonly used together with the `Gotools`
#' family.
#'
#' It uses `requireNamespace()` to test whether each package is loadable.
#' Packages that pass the check are not reinstalled. When `attach = TRUE`,
#' core analysis packages are attached quietly at the end of the check.
#'
#' `BiocManager` is installed automatically when Bioconductor packages need
#' action.
#'
#' @return Invisibly returns `TRUE` when all required packages are available
#'   after the check, or `FALSE` when installation is cancelled.
#'
#' @examples
#' \dontrun{
#' Gotool_dependency()
#' }
#'
#' @export
Gotool_dependency <- function(ask = interactive(), attach = TRUE) {
  .gotools_signature()

  deps <- list(
    cran = c(
      "rmarkdown", "magick", "ape", "Boruta", "car", "cluster", "CLME",
      "compositions", "cowplot", "crayon", "caret", "colorspace", "digest",
      "data.table", "devtools", "doParallel", "ellipse", "emmeans", "e1071",
      "gplots", "ggplot2", "gridExtra", "ggrepel", "doRNG", "ggalluvial",
      "ggExtra",
      "ggforce", "Hmisc", "irlba", "huge", "irr", "lme4", "lmerTest", "nnet",
      "MLmetrics", "Matrix", "magrittr", "MASS", "missForest", "nlme",
      "phangorn", "pheatmap", "pkgconfig", "dplyr", "pscl", "plotly",
      "patchwork", "pdftools", "pROC", "PRROC", "ranger", "rlang",
      "randomForest", "readr", "readxl", "RColorBrewer", "ROCR", "reshape",
      "reshape2", "scales", "stringr", "survminer", "tibble", "tidyr",
      "tidyselect", "vctrs", "yarrr", "stringi", "tidyverse", "vegan",
      "VGAM", "picante", "xgboost", "zoo", "RcppZiggurat", "Rfast",
      "survival", "withr", "knitr", "kableExtra", "DT", "CVXR", "lightgbm",
      "plyr"
    ),
    bioc = c(
      "phyloseq", "microbiome", "Rhtslib", "dada2", "ggpubr", "ggfortify",
      "genefilter", "ggpmisc", "S4Vectors", "ShortRead", "illuminaio",
      "rstatix", "useful", "DECIPHER", "ComplexHeatmap", "DESeq2", "ALDEx2",
      "scater", "ANCOMBC"
    )
  )

  attach_pkgs <- c(
    "ape", "dplyr", "ggplot2", "grid", "patchwork", "phyloseq", "tidyr",
    "vegan", "circlize", "broom", "glue"
  )

  attach_quietly <- function(pkgs) {
    attached_now <- character(0)

    for (pkg in pkgs) {
      if (!requireNamespace(pkg, quietly = TRUE)) next
      pkg_search_name <- paste0("package:", pkg)
      if (pkg_search_name %in% search()) {
        attached_now <- c(attached_now, pkg)
        next
      }
      suppressPackageStartupMessages(
        library(pkg, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)
      )
      attached_now <- c(attached_now, pkg)
    }

    attached_now
  }

  installed_names <- rownames(utils::installed.packages())

  not_loadable_cran <- deps$cran[!vapply(deps$cran, requireNamespace, logical(1), quietly = TRUE)]
  not_loadable_bioc <- deps$bioc[!vapply(deps$bioc, requireNamespace, logical(1), quietly = TRUE)]

  fresh_cran <- not_loadable_cran[!not_loadable_cran %in% installed_names]
  broken_cran <- not_loadable_cran[not_loadable_cran %in% installed_names]
  fresh_bioc <- not_loadable_bioc[!not_loadable_bioc %in% installed_names]
  broken_bioc <- not_loadable_bioc[not_loadable_bioc %in% installed_names]

  loaded_namespaces <- loadedNamespaces()
  loaded_pkgs_needing_repair <- intersect(
    c(broken_cran, broken_bioc),
    loaded_namespaces
  )

  all_needing_action <- c(fresh_cran, broken_cran, fresh_bioc, broken_bioc)

  if (length(all_needing_action) == 0) {
    message("[Gotools] All dependencies are installed and loadable.")
    if (attach) {
      attached_now <- attach_quietly(attach_pkgs)
      if (length(attached_now) > 0) {
        message("[Gotools] Attached essential libraries: ", paste(attached_now, collapse = ", "))
      }
    }
    return(invisible(TRUE))
  }

  if (length(c(fresh_cran, fresh_bioc)) > 0) {
    message("[Gotools] Missing packages: ", paste(c(fresh_cran, fresh_bioc), collapse = ", "))
  }
  if (length(c(broken_cran, broken_bioc)) > 0) {
    message(
      "[Gotools] Installed but not loadable (will reinstall): ",
      paste(c(broken_cran, broken_bioc), collapse = ", ")
    )
  }

  if (length(loaded_pkgs_needing_repair) > 0) {
    stop(
      "[Gotools] Some packages need reinstall but are currently loaded: ",
      paste(loaded_pkgs_needing_repair, collapse = ", "),
      ".\nRestart the R session, then run Gotool_dependency() again."
    )
  }

  if (ask) {
    answer <- readline("[Gotools] Install/repair packages now? [y/N] ")
    if (!tolower(trimws(answer)) %in% c("y", "yes")) {
      message("[Gotools] Installation cancelled.")
      return(invisible(FALSE))
    }
  }

  if (length(c(fresh_bioc, broken_bioc)) > 0 && !requireNamespace("BiocManager", quietly = TRUE)) {
    message("[Gotools] Installing BiocManager first...")
    utils::install.packages("BiocManager", quiet = TRUE)
  }

  if ("ANCOMBC" %in% c(fresh_bioc, broken_bioc) || "CVXR" %in% c(fresh_cran, broken_cran)) {
    cargo_found <- nzchar(Sys.which("cargo"))
    if (!cargo_found) {
      cargo_paths <- c(
        path.expand("~/.cargo/bin/cargo"),
        "/usr/local/bin/cargo",
        "/opt/homebrew/bin/cargo"
      )
      cargo_found <- any(file.exists(cargo_paths))
    }
    if (!cargo_found) {
      stop(
        "[Gotools] ANCOMBC requires the 'clarabel' package via CVXR, which may need Rust source compilation.\n",
        "Rust toolchain (cargo) was not found on this system.\n\n",
        "Install Rust by running this command in your Terminal:\n",
        "  curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh\n\n",
        "After installation, restart R and run Gotool_dependency() again.\n",
        "If Rust is already installed, make sure 'cargo' is on your PATH."
      )
    }
  }

  if (!requireNamespace("terra", quietly = TRUE)) {
    message(
      "[Gotools] Note: 'terra' is missing but is not a direct Gotools dependency.\n",
      "  If installation later fails due to terra, install system libraries first:\n",
      "    brew install libomp\n",
      "  Then set ~/.R/Makevars:\n",
      "    LDFLAGS += -L/opt/homebrew/opt/libomp/lib -lomp\n",
      "    CPPFLAGS += -I/opt/homebrew/opt/libomp/include -Xpreprocessor -fopenmp\n",
      "  Then retry: install.packages(\"terra\", repos = \"https://rspatial.r-universe.dev\")"
    )
  }

  needs_ancombc <- "ANCOMBC" %in% c(fresh_bioc, broken_bioc)
  if (needs_ancombc) {
    cvxr_ok <- tryCatch({
      if (requireNamespace("CVXR", quietly = TRUE)) {
        "solve" %in% getNamespaceExports("CVXR")
      } else {
        FALSE
      }
    }, error = function(e) FALSE)

    if (!cvxr_ok) {
      if (!requireNamespace("remotes", quietly = TRUE)) {
        message("[Gotools] Installing remotes (needed for CVXR version pinning)...")
        utils::install.packages("remotes", quiet = TRUE)
      }
      message("[Gotools] Installing CVXR 1.0-11 (required for ANCOMBC compatibility)...")
      remotes::install_version(
        "CVXR",
        version = "1.0-11",
        quiet = TRUE,
        repos = "https://cloud.r-project.org"
      )
    }
  }

  if (length(fresh_cran) > 0) {
    message("[Gotools] Installing CRAN packages: ", paste(fresh_cran, collapse = ", "))
    utils::install.packages(fresh_cran, quiet = TRUE)
  }
  if (length(broken_cran) > 0) {
    message("[Gotools] Reinstalling broken CRAN packages: ", paste(broken_cran, collapse = ", "))
    utils::install.packages(broken_cran, quiet = TRUE)
  }

  if (length(fresh_bioc) > 0) {
    message("[Gotools] Installing Bioconductor packages: ", paste(fresh_bioc, collapse = ", "))
    BiocManager::install(fresh_bioc, ask = FALSE, update = FALSE, quiet = TRUE)
  }
  if (length(broken_bioc) > 0) {
    message("[Gotools] Reinstalling broken Bioconductor packages: ", paste(broken_bioc, collapse = ", "))
    BiocManager::install(broken_bioc, ask = FALSE, update = FALSE, force = TRUE, quiet = TRUE)
  }

  remaining_cran <- deps$cran[!vapply(deps$cran, requireNamespace, logical(1), quietly = TRUE)]
  remaining_bioc <- deps$bioc[!vapply(deps$bioc, requireNamespace, logical(1), quietly = TRUE)]
  remaining <- c(remaining_cran, remaining_bioc)

  if (length(remaining) > 0) {
    stop("[Gotools] Some dependencies are still not loadable: ", paste(remaining, collapse = ", "))
  }

  if (requireNamespace("crayon", quietly = TRUE)) {
    cat(crayon::yellow("All required Gotools dependencies are available.\n"))
    cat(crayon::blue("#--------------------------------------------------------------# \n"))
  } else {
    cat("All required Gotools dependencies are available.\n")
    cat("#--------------------------------------------------------------# \n")
  }

  if (attach) {
    attached_now <- attach_quietly(attach_pkgs)
    if (length(attached_now) > 0) {
      message("[Gotools] Attached essential libraries: ", paste(attached_now, collapse = ", "))
    }
  }

  invisible(TRUE)
}
