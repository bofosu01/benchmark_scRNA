# Install all packages needed to run the scRNA-seq analysis

# Define the list of CRAN packages
cran_packages <- c(
  "ggplot2",    # plotting
  "dplyr",      # data manipulation
  "tidyr",      # reshaping data
  "readr",      # reading/writing tabular files
  "purrr",      # functional programming helpers
  "stringr",    # string manipulation
  "Matrix",     # sparse matrix support
  "data.table", # fast data frames
  "BiocManager" # manages Bioconductor package installation
)

# Install the CRAN packages
install.packages(cran_packages)

# Switch Bioconductor to the development branch so Splatter is available
BiocManager::install(version = "devel")

# Define the list of Bioconductor packages
bioc_packages <- c(
  "MAST",                 # hurdle-model DE method
  "splatter",             # simulation of scRNA-seq counts
  "SingleCellExperiment", # container class for single-cell data
  "scater"                # single-cell QC and plotting utilities
)

# Install the Bioconductor packages without prompting or updating old ones
BiocManager::install(bioc_packages, update = FALSE, ask = FALSE)

# Install remotes, which fetches packages from GitHub 
install.packages("remotes")

# Install NEBULA and seurat from its GitHub repository (negative-binomial mixed model DE)
remotes::install_github("lhe17/nebula")
 #"Seurat for scRNA-seq analysis framework (also used for Wilcoxon DE)
remotes::install_github("satijalab/seurat", ref = "main")