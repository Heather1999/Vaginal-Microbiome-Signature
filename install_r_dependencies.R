# R dependencies for permanova_maaslin2.R
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

BiocManager::install("Maaslin2")
install.packages(c("vegan", "dplyr", "readr"))
