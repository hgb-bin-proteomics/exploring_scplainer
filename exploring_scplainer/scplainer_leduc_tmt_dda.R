library(tidyverse)
library(scp)
library(scpdata)
library(patchwork)
library(ensembldb)
library(EnsDb.Hsapiens.v86)
library(bluster)
library(scater)

sce <- readr::read_rds("sce.RDS")

# We follow the workflow as documented here:
# https://uclouvain-cbio.github.io/SCP.replication/articles/scplainer_leduc2022.html