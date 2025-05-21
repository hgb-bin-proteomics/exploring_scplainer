library(tidyverse)
library(scp)
library(scpdata)
library(patchwork)
library(ensembldb)
library(EnsDb.Hsapiens.v86)
library(bluster)
library(scater)

### We follow the workflow as documented here:
### https://uclouvain-cbio.github.io/SCP.replication/articles/scplainer_leduc2022.html
### with some extra details sprinkled in from here:
### https://uclouvain-cbio.github.io/scp/articles/scp_data_modelling.html

# The data here is already preprocessed (see documentation).
# This involves:
# - Data cleaning
#   + Replacing zeros with missing values (NA)
# - Feature quality control
#   + Remove contaminants
#   + Remove decoys
#   + Remove PSMs with low spectral purity
#   + Remove peptides not within the FDR threshold
#   + Remove peptides above 5% single cell to carrier ratio
# - Sample quality control
#   + Remove samples (cells) with low number of peptides
#   + Remove samples with low median intensity
#   + Remove samples with high median CV
#   + (Remove samples we are not interested in, e.g. negative controls, carriers, references, empty wells)
# - Peptide data assembly
#   + Aggregating PSMs to peptides
# - Log2-transformation

sce <- readr::read_rds("sce.RDS")

## Data Overview

# The data set primarily consists of melanoma cells and monocytes.
# The data set also contains carrier samples, negative control samples, reference samples and empty wells (Unused).
# Carrier samples, negative control samples, reference samples and empty wells are already filtered out here, because
# the data is preprocessed.

table(sce$SampleType)

## The data were acquired using TMT-18 labelling (and in DDA mode).

levels(sce$Channel)

# The data were acquired as part of 134 MS acquisition batches.
# 4 were already filtered out in preprocessing.

length(unique(sce$Set))

cat(head(sce$Set), "...", tail(sce$Set))

# Finally, samples were prepared with the nPOP protocol using 2 glass slides. The information is stored under lcbatch.

table(sce$lcbatch)

## Fitting the Regression Model

# First we need to define a formula for the model to fit and choose which variables we want to include.
# Here we use the following:
# - MedianIntensity: The normalization factor used to correct for cell-specific technical differences.
# - Channel: Used to correct for TMT effects.
# - Set: Used to perform batch correction. We consider each acquisition run to be a batch.
# - SampleType: The biological variable of interest. It capture the difference between macrophages and monocytes.

f <- ~ 1 + ## intercept
  ## normalization
  MedianIntensity +
  ## batch effects
  Channel + Set + 
  ## biological variability
  SampleType

# Now we fit the regression model using our formula.

# Un-comment to re-run the model fit
#sce <- scpModelWorkflow(sce, formula = f)
#readr::write_rds(sce, "sce_fitted.RDS", compress = "xz", version = 2, text = F)

# Please download the fitted model from: http://u.pc.cd/rUN
sce <- readr::read_rds("sce_fitted.RDS")

# You can always retrieve the formula that was used to fit model with:

scpModelFormula(sce)

# The data that is modelled by each variable are contained in the so-called effect matrices:

(effects <- scpModelEffects(sce))

# We can investigate the size of these effect matrices with dim()
# We get:
# - 16 670 rows (peptides)
# - 1 656 columns (cells)

dim(effects$MedianIntensity)

# These dimensions are the same as our model input data:

dim((model_input <- scpModelInput(sce)))

# However, they differ from our original data due to peptide filtering that is done
# before the model is fitted (?filtering peptides that have too many missing values to be modelled)

dim(sce)

# Data that could not be modelled by our variables is captured in the residual matrix.

dim((residual <- scpModelResiduals(sce)))

