# Install additional R packages needed for sc-type
if (!require("HGNChelper", quietly = TRUE)) {
  install.packages("HGNChelper", repos="https://cran.rstudio.com/", dependencies=TRUE)
}

if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos="https://cran.rstudio.com/", dependencies=TRUE)
}

# Test loading packages
library(HGNChelper)
library(dplyr)
library(openxlsx)

# Test checkGeneSymbols function
test_genes <- c("CD3D", "CD19", "CD68")
result <- checkGeneSymbols(test_genes)
print("✅ checkGeneSymbols working:")
print(result)

cat("✅ All required R packages installed and tested successfully!\n")