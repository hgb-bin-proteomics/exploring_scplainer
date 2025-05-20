req <- c(
  "UCLouvain-CBIO/scp",
  "UCLouvain-CBIO/scpdata",
  "ensembldb",
  "EnsDb.Hsapiens.v86",
  "bluster",
  "scater",
  "tidyverse",
  "patchwork"
)
renv::install(req)
renv::snapshot()
