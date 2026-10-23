#!/usr/bin/env Rscript
# ============================================================
# 00 - Install Required R Packages
# Clonal Gene Expression Shift
# Eric Wafula | 2026
# ============================================================
#
# Required for:
#   Standalone clonal-expression-shifts workflow.
#   Adapted from the original parent-project installer; original preserved.
#
# Check without installing:
#   Rscript 00_install_packages.R --check_only TRUE
#
# Core packages:
#   tidyverse
#   optparse
#   Seurat
#   glmnet
#   ape
#   Matrix
#
# Notes:
#   - For Seurat, we install only Depends/Imports/LinkingTo dependencies.
#     This avoids trying to install many optional Suggests packages.
#   - Optional Seurat helper packages such as BPCells, presto, and glmGamPoi
#     are not required for this Elastic Net analysis.
#
# Usage:
#   Rscript 00_install_packages.R
#
# Optional local library:
#   Rscript 00_install_packages.R --lib R_libs
#
# Optional Seurat extension packages:
#   Rscript 00_install_packages.R --install_seurat_extensions TRUE
#
# ============================================================

suppressPackageStartupMessages({
  if (!requireNamespace("utils", quietly = TRUE)) {
    stop("The 'utils' package is required but not available.")
  }
})

# ============================================================
# Simple command-line argument parser
# ============================================================

args <- commandArgs(trailingOnly = TRUE)

lib_path <- NULL
install_seurat_extensions <- FALSE
check_only <- FALSE

if (length(args) > 0) {
  if ("--check_only" %in% args) {
    idx <- match("--check_only", args)
    if (idx == length(args) || !tolower(args[idx + 1L]) %in% c("true", "false"))
      stop("--check_only must be TRUE or FALSE.")
    check_only <- tolower(args[idx + 1L]) == "true"
  }

  if ("--lib" %in% args) {
    lib_idx <- match("--lib", args)

    if (lib_idx == length(args)) {
      stop("--lib was provided but no library path followed it.")
    }

    lib_path <- args[lib_idx + 1]
  }

  if ("--install_seurat_extensions" %in% args) {
    ext_idx <- match("--install_seurat_extensions", args)

    if (ext_idx == length(args)) {
      stop("--install_seurat_extensions was provided but TRUE/FALSE did not follow it.")
    }

    ext_value <- tolower(args[ext_idx + 1])

    if (!ext_value %in% c("true", "false", "t", "f", "1", "0", "yes", "no")) {
      stop("--install_seurat_extensions must be TRUE or FALSE.")
    }

    install_seurat_extensions <- ext_value %in% c("true", "t", "1", "yes")
  }

  allowed_args <- c(
    "--check_only",
    "--lib",
    lib_path,
    "--install_seurat_extensions",
    as.character(install_seurat_extensions),
    "TRUE",
    "FALSE",
    "true",
    "false"
  )

  unknown_args <- setdiff(args, allowed_args)

  if (length(unknown_args) > 0) {
    stop(
      "Unknown argument(s): ",
      paste(unknown_args, collapse = ", "),
      "\nSupported usage:\n",
      "  Rscript 00_install_packages.R\n",
      "  Rscript 00_install_packages.R --lib R_libs\n",
      "  Rscript 00_install_packages.R --install_seurat_extensions TRUE"
    )
  }
}

# ============================================================
# CRAN mirror
# ============================================================

cran_repo <- "https://cloud.r-project.org"
options(repos = c(CRAN = cran_repo))

message("\n========== R PACKAGE INSTALLATION CHECK ==========")
message("[INFO] CRAN mirror: ", cran_repo)

# ============================================================
# Optional local library setup
# ============================================================

if (!is.null(lib_path)) {
  message("\n[INFO] Requested R library path: ", lib_path)

  dir.create(lib_path, recursive = TRUE, showWarnings = FALSE)

  lib_path <- normalizePath(lib_path, mustWork = TRUE)
  .libPaths(c(lib_path, .libPaths()))

  message("[INFO] Updated .libPaths():")
  print(.libPaths())
} else {
  message("\n[INFO] No --lib path provided. Using current .libPaths():")
  print(.libPaths())
}

# ============================================================
# Required packages
# ============================================================

required_cran_packages <- c(
  "tidyverse",
  "optparse",
  "Seurat",
  "SeuratObject",
  "glmnet",
  "ape",
  "Matrix",
  "phangorn",
  "igraph",
  "lars"
)

required_cran_packages <- unique(required_cran_packages)

# ============================================================
# Helper functions
# ============================================================

is_installed <- function(pkg) {
  requireNamespace(pkg, quietly = TRUE)
}

install_package_safe <- function(pkg) {
  if (is_installed(pkg)) {
    message("[OK] ", pkg, " is already installed.")
    return(TRUE)
  }

  if (check_only) {
    message("[MISSING] ", pkg, " (check-only; no installation attempted)")
    return(FALSE)
  }

  message("[INSTALL] Installing missing package: ", pkg)

  # Use a conservative dependency setting.
  # dependencies = TRUE installs Depends, Imports, LinkingTo, Suggests, and Enhances.
  # For Seurat, that tries many optional packages that are not required here.
  deps_to_install <- c("Depends", "Imports", "LinkingTo")

  install_ok <- tryCatch(
    {
      utils::install.packages(
        pkg,
        dependencies = deps_to_install
      )
      TRUE
    },
    error = function(e) {
      message("[ERROR] Failed to install package: ", pkg)
      message("        Error message: ", e$message)
      FALSE
    }
  )

  if (!install_ok) {
    return(FALSE)
  }

  if (is_installed(pkg)) {
    message("[SUCCESS] Installed package: ", pkg)
    return(TRUE)
  }

  message("[FAILED] Package still not available after install attempt: ", pkg)
  return(FALSE)
}

install_seurat_optional_extensions <- function() {
  message("\n========== OPTIONAL SEURAT EXTENSIONS ==========")

  message("[INFO] Adding Seurat/BPCells r-universe repositories.")

  old_repos <- getOption("repos")

  options(
    repos = c(
      CRAN = cran_repo,
      satijalab = "https://satijalab.r-universe.dev",
      bnprks = "https://bnprks.r-universe.dev"
    )
  )

  optional_pkgs <- c("BPCells", "presto", "glmGamPoi")

  status <- vapply(
    optional_pkgs,
    function(pkg) {
      if (is_installed(pkg)) {
        message("[OK] Optional package already installed: ", pkg)
        return(TRUE)
      }

      message("[INSTALL] Installing optional Seurat extension: ", pkg)

      ok <- tryCatch(
        {
          utils::install.packages(
            pkg,
            dependencies = c("Depends", "Imports", "LinkingTo")
          )
          TRUE
        },
        error = function(e) {
          message("[WARNING] Could not install optional package: ", pkg)
          message("          Error message: ", e$message)
          FALSE
        }
      )

      if (ok && is_installed(pkg)) {
        message("[SUCCESS] Installed optional package: ", pkg)
        return(TRUE)
      }

      message("[WARNING] Optional package not available after install attempt: ", pkg)
      FALSE
    },
    logical(1)
  )

  options(repos = old_repos)

  message("\n[INFO] Optional extension package status:")
  print(
    data.frame(
      package = names(status),
      installed = as.logical(status),
      row.names = NULL
    )
  )

  invisible(status)
}

# ============================================================
# Install missing packages
# ============================================================

message("\n========== CHECKING REQUIRED PACKAGES ==========")

install_status <- vapply(
  required_cran_packages,
  install_package_safe,
  logical(1)
)

# Optional Seurat performance/helper packages
if (install_seurat_extensions && !check_only) {
  install_seurat_optional_extensions()
} else {
  message("\n[INFO] Skipping optional Seurat extensions.")
  message("       To install BPCells, presto, and glmGamPoi, run:")
  message("       Rscript 00_install_packages.R --install_seurat_extensions TRUE")
}

# ============================================================
# Final validation
# ============================================================

message("\n========== FINAL PACKAGE CHECK ==========")

final_status <- data.frame(
  package = required_cran_packages,
  installed = vapply(required_cran_packages, is_installed, logical(1)),
  stringsAsFactors = FALSE
)

print(final_status, row.names = FALSE)

missing_after_install <- final_status$package[!final_status$installed]

if (length(missing_after_install) > 0) {
  message("\n[WARNING] Some required packages could not be installed or loaded:")
  message("  ", paste(missing_after_install, collapse = ", "))

  message("\nRun Rscript 00_install_packages.R to install missing dependencies, or use the fixes below.")
  message("\nPossible fixes:")
  message("  1. Check internet access from this compute node.")
  message("  2. Use a writable local R library:")
  message("       mkdir -p R_libs")
  message("       Rscript 00_install_packages.R --lib R_libs")
  message("  3. Install Seurat directly if needed:")
  message("       install.packages('Seurat')")
  message("  4. If Bioconductor packages are desired for optional Seurat features:")
  message("       install.packages('BiocManager')")
  message("       BiocManager::install(c('glmGamPoi', 'SingleCellExperiment'))")
  message("  5. If Matrix fails, check R version compatibility with CRAN Matrix.")

  quit(status = 1)
}

message("\n[SUCCESS] All required packages are installed and available.")

# ============================================================
# Load test
# ============================================================

message("\n========== LOAD TEST ==========")

for (pkg in required_cran_packages) {
  message("[LOAD] Testing library(", pkg, ")")

  ok <- tryCatch(
    {
      suppressPackageStartupMessages(
        library(pkg, character.only = TRUE)
      )
      TRUE
    },
    error = function(e) {
      message("[ERROR] Could not load package: ", pkg)
      message("        Error message: ", e$message)
      FALSE
    }
  )

  if (!ok) {
    message("\n[FAILED] Package installed but could not be loaded: ", pkg)
    quit(status = 1)
  }
}

message("[SUCCESS] All required packages loaded successfully.")

# ============================================================
# Session information
# ============================================================

message("\n========== SESSION INFO ==========")
print(sessionInfo())

message("\n========== INSTALLATION CHECK COMPLETE ==========")
