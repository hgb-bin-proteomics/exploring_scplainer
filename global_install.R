if (!requireNamespace("BiocManager"))
  install.packages("BiocManager")
BiocManager::install(version = "devel")
BiocManager::install("remotes")
stopifnot(BiocManager::valid())
BiocManager::install("UCLouvain-CBIO/scp")
BiocManager::install("UCLouvain-CBIO/scpdata")
BiocManager::install("ensembldb")
BiocManager::install("EnsDb.Hsapiens.v86")
BiocManager::install("bluster")
BiocManager::install("scater")
install.packages("patchwork")
install.packages("tidyverse")