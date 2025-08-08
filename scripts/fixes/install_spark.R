# Install SPARK R package
if (!require("SPARK", quietly = TRUE)) {
  install.packages("SPARK", repos = "https://cran.r-project.org/")
  cat("SPARK installed successfully\n")
} else {
  cat("SPARK already installed\n")
}

# Test loading the package
library(SPARK)
cat("SPARK package loaded successfully\n")