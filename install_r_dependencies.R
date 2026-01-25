#!/usr/bin/env Rscript
# ChatSpatial R Dependencies Installation Script
# This script installs all R packages required for ChatSpatial's R-based methods
#
# Usage:
#   Rscript install_r_dependencies.R
#   or in R console: source("install_r_dependencies.R")

cat("\n")
cat("========================================\n")
cat("ChatSpatial R Dependencies Installer\n")
cat("========================================\n\n")

# Set CRAN mirror
options(repos = c(CRAN = "https://cran.r-project.org"))

# Helper function to install package with error handling
install_if_missing <- function(pkg, source = "CRAN", install_cmd = NULL) {
  cat(sprintf("Checking %s (%s)...", pkg, source))

  if (require(pkg, character.only = TRUE, quietly = TRUE)) {
    cat(" ✓ Already installed\n")
    return(TRUE)
  }

  cat(" Installing...\n")

  tryCatch({
    if (!is.null(install_cmd)) {
      eval(parse(text = install_cmd))
    } else {
      install.packages(pkg, quiet = FALSE)
    }

    # Verify installation
    if (require(pkg, character.only = TRUE, quietly = TRUE)) {
      cat(sprintf("  ✓ %s installed successfully\n", pkg))
      return(TRUE)
    } else {
      cat(sprintf("  ✗ %s installation verification failed\n", pkg))
      return(FALSE)
    }
  }, error = function(e) {
    cat(sprintf("  ✗ Failed to install %s: %s\n", pkg, e$message))
    return(FALSE)
  })
}

# Track installation results
failed_packages <- c()
success_count <- 0
total_count <- 0

cat("\n")
cat("Step 1: Installing Base Dependencies\n")
cat("-------------------------------------\n")

# Install devtools (needed for GitHub packages)
total_count <- total_count + 1
if (install_if_missing("devtools", "CRAN")) {
  success_count <- success_count + 1
} else {
  failed_packages <- c(failed_packages, "devtools (CRAN)")
}

# Install BiocManager (needed for Bioconductor packages)
total_count <- total_count + 1
if (install_if_missing("BiocManager", "CRAN")) {
  success_count <- success_count + 1
} else {
  failed_packages <- c(failed_packages, "BiocManager (CRAN)")
}

cat("\n")
cat("Step 2: Installing CRAN Packages\n")
cat("---------------------------------\n")

cran_packages <- c(
  "dplyr",        # Data manipulation (required by scType, Numbat, CellChat)
  "openxlsx",     # Excel file reading (required by scType)
  "HGNChelper",   # Gene name validation (required by scType)
  "sctransform",  # SCTransform v2 normalization (variance-stabilizing)
  "Matrix",       # Sparse matrix operations (required by SCTransform)
  "mclust"        # Model-based clustering (required by GraphST)
)

for (pkg in cran_packages) {
  total_count <- total_count + 1
  if (install_if_missing(pkg, "CRAN")) {
    success_count <- success_count + 1
  } else {
    failed_packages <- c(failed_packages, paste0(pkg, " (CRAN)"))
  }
}

cat("\n")
cat("Step 3: Installing Bioconductor Packages\n")
cat("-----------------------------------------\n")

bioc_packages <- list(
  # SPOTlight and its dependencies
  list(name = "SPOTlight", cmd = "BiocManager::install('SPOTlight', update = FALSE, ask = FALSE)"),
  list(name = "SingleCellExperiment", cmd = "BiocManager::install('SingleCellExperiment', update = FALSE, ask = FALSE)"),
  list(name = "SpatialExperiment", cmd = "BiocManager::install('SpatialExperiment', update = FALSE, ask = FALSE)"),
  list(name = "scran", cmd = "BiocManager::install('scran', update = FALSE, ask = FALSE)"),
  list(name = "scuttle", cmd = "BiocManager::install('scuttle', update = FALSE, ask = FALSE)")
)

for (pkg_info in bioc_packages) {
  total_count <- total_count + 1
  if (install_if_missing(pkg_info$name, "Bioconductor", pkg_info$cmd)) {
    success_count <- success_count + 1
  } else {
    failed_packages <- c(failed_packages, paste0(pkg_info$name, " (Bioconductor)"))
  }
}

cat("\n")
cat("Step 4: Installing GitHub Packages\n")
cat("-----------------------------------\n")
cat("Note: GitHub installations may take longer...\n\n")

github_packages <- list(
  list(
    name = "spacexr",
    repo = "dmcable/spacexr",
    method = "RCTD deconvolution",
    cmd = "devtools::install_github('dmcable/spacexr', build_vignettes = FALSE, upgrade = 'never')"
  ),
  list(
    name = "CARD",
    repo = "YMa-lab/CARD",
    method = "CARD deconvolution",
    cmd = "devtools::install_github('YMa-lab/CARD', upgrade = 'never')"
  ),
  list(
    name = "CellChat",
    repo = "jinworks/CellChat",
    method = "Cell-cell communication",
    cmd = "devtools::install_github('jinworks/CellChat', upgrade = 'never')"
  ),
  list(
    name = "numbat",
    repo = "kharchenkolab/numbat",
    method = "CNV analysis",
    cmd = "devtools::install_github('kharchenkolab/numbat', upgrade = 'never')"
  ),
  list(
    name = "SPARK",
    repo = "xzhoulab/SPARK",
    method = "Spatial variable genes (SPARK-X)",
    cmd = "devtools::install_github('xzhoulab/SPARK', upgrade = 'never')"
  )
)

for (pkg_info in github_packages) {
  cat(sprintf("Installing %s (%s) from GitHub...\n", pkg_info$name, pkg_info$method))
  total_count <- total_count + 1
  if (install_if_missing(pkg_info$name, paste0("GitHub: ", pkg_info$repo), pkg_info$cmd)) {
    success_count <- success_count + 1
  } else {
    failed_packages <- c(failed_packages, paste0(pkg_info$name, " (GitHub: ", pkg_info$repo, ")"))
  }
}

# Print summary
cat("\n")
cat("========================================\n")
cat("Installation Summary\n")
cat("========================================\n")
cat(sprintf("Successfully installed: %d/%d packages\n", success_count, total_count))

if (length(failed_packages) > 0) {
  cat("\n⚠️  Failed packages:\n")
  for (pkg in failed_packages) {
    cat(sprintf("  - %s\n", pkg))
  }
  cat("\n")
  cat("Troubleshooting:\n")
  cat("1. Check your internet connection\n")
  cat("2. Update R to the latest version\n")
  cat("3. Try installing failed packages manually\n")
  cat("4. Check the ChatSpatial documentation for platform-specific instructions\n")
  cat("\n")
} else {
  cat("\n✓ All R dependencies installed successfully!\n\n")
  cat("You can now use ChatSpatial's R-based methods:\n")
  cat("  • RCTD deconvolution (spacexr)\n")
  cat("  • SPOTlight deconvolution (SPOTlight + dependencies)\n")
  cat("  • CARD deconvolution (CARD)\n")
  cat("  • CellChat cell communication (CellChat)\n")
  cat("  • SCTransform normalization (sctransform, Matrix)\n")
  cat("  • scType cell type annotation (dplyr, openxlsx, HGNChelper)\n")
  cat("  • Numbat CNV analysis (numbat)\n")
  cat("  • SPARK-X spatial variable genes (SPARK)\n")
  cat("  • GraphST clustering (mclust)\n")
  cat("\n")
}

# Print R session info for debugging
cat("R Session Information:\n")
cat("---------------------\n")
cat(sprintf("R version: %s\n", R.version.string))
cat(sprintf("Platform: %s\n", R.version$platform))

cat("\nInstallation complete!\n")
