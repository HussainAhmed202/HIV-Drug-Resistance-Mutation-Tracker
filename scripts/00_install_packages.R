# 00_install_packages.R

# CRAN packages
install.packages(c(
  "rentrez",
  "ape",
  "ggplot2"
))

# Bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install(c(
  "Biostrings",
  "DECIPHER",
))

# Verify everything loaded correctly
library(rentrez)
library(ape)
library(ggplot2)
library(Biostrings)
library(DECIPHER)

cat("All packages loaded successfully.\n")