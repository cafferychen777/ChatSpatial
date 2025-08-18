# Install required R packages for sc-type
if (!require("dplyr", quietly = TRUE)) {
  install.packages("dplyr", repos="https://cran.rstudio.com/", dependencies=TRUE)
}

if (!require("openxlsx", quietly = TRUE)) {
  install.packages("openxlsx", repos="https://cran.rstudio.com/", dependencies=TRUE)
}

# Test if packages loaded successfully
library(dplyr)
library(openxlsx)

cat("✅ R packages installed successfully!\n")