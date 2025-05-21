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

# Where does this data come from: look at "preprocess_leduc_tmt_dda.R"
# The steps are explained here: https://uclouvain-cbio.github.io/SCP.replication/articles/scplainer_leduc2022.html#minimal-data-processing
sce <- readr::read_rds("sce.RDS")

## The SingleCellExperiment Class

# The SingleCellExperiment class stores all the necessary information we need to analyse data with
# SCPlainer.
# You can read more about the SingleCellExperiment class here:
# - https://bioconductor.org/packages/release/bioc/html/SingleCellExperiment.html
# -

# Essentially the the SingleCellExperiment consists of three main data objects, that we need to define:
# - Assay data (one or more)
#   + In the proteomics context this is a matrix of quantification values
#     The matrix consists of:
#     - Rows: PSMS / peptides / precursors / proteins
#     - Columns: samples / cells
#     - Values: quantification, can be normalized/transformed/etc.
# - colData
#   + Meta-data about the columns, so sample / cell information. For example:
#     - Acquisition date
#     - Acquisition batch
#     - Sample / cell type
#     - Channel
#     - Technical descriptors
#   + The colnames of Assay data need to match with the rownames of the colData.
#     The colData matrix basically consists of:
#     - Rows: samples / cells
#     - Columns: Descriptors of the information, e.g. Acquisition batch, Channel, etc.
# - rowData
#   + Meta-data about the rows, so PSM / peptide / precursor / protein information. For example:
#     - CScore
#     - qvalue
#     - Protein group
#     - Gene
#   + The rownames of Assay data need to match with the rownames of the rowData.
#     The rowData matrix basically consists of:
#     - Rows: PSMs / peptides / precursors / proteins
#     - Columns: Descriptors of the information, e.g. Protein group, qvalue

# See also "SingleCellExperiment.[pptx|pdf]"

# Reading of Spectronaut data into a SingleCellExperiment is demonstrated in:
# "read_SCE.R"

## Data Overview

# The data set primarily consists of melanoma cells and monocytes.
# The data set also contains carrier samples, negative control samples, reference samples and empty wells (Unused).
# Carrier samples, negative control samples, reference samples and empty wells are already filtered out here, because
# the data is preprocessed.

table(sce$SampleType)

# The data were acquired using TMT-18 labelling (and in DDA mode).

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

## Peptide filtering

# Our data contains 23 637 548 missing values and 3 967 972 non-missing values
# -> 85% of our data is missing

table(missing = is.na(assay(sce)))

# We can compute the n over p ration for each peptide as the following:

head((np_peptides <- scpModelFilterNPRatio(sce)))

# The n over p ratio is:
# - n: number of cells or samples with measured intensity
# - p: number of model coefficients to estimate
# In our case p = 4, so any peptide that is observed in less than 4 cells would not satisfy
# the n/p > 1 threshold:

length(np_peptides)

# We can visualize np ratio like this:

scpModelFilterPlot(sce)

# We can also limit our data to peptides with np ratio > 3 by:

scpModelFilterThreshold(sce) <- 3
scpModelFilterPlot(sce)

# As you can see, our data now shrunk down to 6854 peptides:

length((np_peptides <- scpModelFilterNPRatio(sce)))

## Model Analysis - Variance

# Variance analysis explores the proportion of data captured by each variable in the model.

(vaRes <- scpVarianceAnalysis(sce))

# We can look at the variance captured by each variable for each peptide like this, here
# demonstrated for variable MedianIntensity. The columns are as follows:
# - feature: the peptide
# - SS: captured variance
# - df: residual degrees of freedom for estimating the variance
# - percentExplainedVar: the percentage of total variance explained

vaRes$MedianIntensity

# We can annotate this data with the information of our original data (sce)
# - by: refers to the column name of our vaRes table
# - by2: refers to the matching column name in the rowData of our original data (sce)

vaRes <- scpAnnotateResults(
  vaRes, rowData(sce), by = "feature", by2 = "Sequence"
)

# We can then explore the variance for our whole dataset:

scpVariancePlot(vaRes)

# We can also explore this for specific peptides, e.g. the top 20 peptides in terms of explained
# variance by sample type:

scpVariancePlot(
  vaRes, top = 20, by = "percentExplainedVar", effect = "SampleType",
  decreasing = TRUE, combined = FALSE
)

# We can also group these peptides by their corresponding proteins

scpVariancePlot(
  vaRes, top = 20, by = "percentExplainedVar", effect = "SampleType",
  decreasing = TRUE, combined = FALSE, fcol = "gene"
)

# Alternatively, we can generate protein level results by aggregating peptide level results:

vaProtein <- scpVarianceAggregate(vaRes, fcol = "gene")
scpVariancePlot(
  vaProtein, effect = "SampleType", top = 10, combined = FALSE
)

## Model Analysis - Differential Abundance

# Next, we explore the model output to understand the differences between melanoma cells and monocytes.
# The difference of interest is specified using the contrast argument:
# - The first element points to the variable to test and 
# - The two following element are the groups of interest to compare.
# - You can provide multiple contrasts in a list.

(daRes <- scpDifferentialAnalysis(
  sce, contrast = list(c("SampleType", "Melanoma", "Monocyte"))
))

# Similarly to variance analysis, the results are a list of tables, one table for each contrast.
# Each table reports for each feature (peptide):
# - Estimate: the estimated difference between the two groups.
# - SE: the standard error associated to the estimation.
# - Df: the degrees of freedom
# - tstatistic: the t-statistics
# - pvalue: the associated p-value
# - padj: the p-value FDR-adjusted for multiple testing across all peptides

daRes$SampleType_Melanoma_vs_Monocyte

# Again, to better explore the results, we add the annotations available in the rowData:

daRes <- scpAnnotateResults(
  daRes, rowData(sce), 
  by = "feature", by2 = "Sequence"
)

# We then visualize the results using a volcano plot. The function returns a volcano plot for each contrast.

scpVolcanoPlot(daRes)

# To help interpretation of the results, we will label the peptides with their protein name.
# -> textBy = "gene"
# Also we increase the number of labels shown on the plot.
# -> top = 30
# Finally, we can add colors to the plot. For instance, let's explore the impact of the number of observations using
# the np ratio. We create a new annotation table, add it to the results and redraw the plot. The np ratio is retrieved
# using scpModelFilterNPRatio that we used previously.

daRes <- scpAnnotateResults(
  daRes, data.frame(feature = names(np_peptides), npRatio = np_peptides), 
  by = "feature"
)
scpVolcanoPlot(
  daRes, top = 30, textBy = "gene", 
  pointParams = list(aes(colour = npRatio))
)

# As expected, higher number of observations (higher n/p) lead to increased statistical power and hence to more significant results.
# Proteins such as VIM, LGALS3, CALU, LMNA, CTTN, are more abundant in melanoma cells compared to monocytes.
# On the other hand, proteins such as LCP1, CORO1A, ARHGDIB, TMSB4X are more abundant in monocytes compared to melanoma cells.

# Finally, we can report results at the protein level:

scpDifferentialAggregate(daRes, fcol = "gene") |> 
  scpVolcanoPlot(top = 30, textBy = "gene")

## Model Analysis - Component Analysis

# We can perform component analysis to link the modelled effects to the cellular heterogeneity. 
# In the following we run an APCA+ (extended ANOVA-simultaneous principal component analysis) for the sample type effect.
# In other words, we perform a PCA on the data that is captured by the sample type variable along with the residuals (unmodelled data).

(caRes <- scpComponentAnalysis(
  sce, ncomp = 2, method = "APCA", effect = "SampleType"
))

# The results are contained in a list with 2 elements:
# - bySample: contains the PC scores, that is the component results in sample space
# - byFeature: contains the eigenvectors, that is the component results in feature space

caRes$bySample

# Each of the two elements contains components results for:
# - the data before modelling (unmodelled)
# - the residuals
# - the APCA on the sample type variable (APCA_SampleType)

caRes$bySample$APCA_SampleType

# We can explore the component analysis results in cell space.
# Similarly to the previous explorations, we annotate the results beforehand:

caResCells <- caRes$bySample # rownames of this are the names of our cells
sce$cell <- colnames(sce) # we create a new meta-data column in our sce colData called cell, which just holds all cell names
caResCells <- scpAnnotateResults(caResCells, colData(sce), by = "cell") # matching rownames with our cell column in colData

# We then can generate the component plot, colouring by SampleType. To assess the impact of batch effects, we shape the points
# according to the plate batch / lcbatch as well:

scpComponentPlot(
  caResCells,
  pointParams = list(aes(colour = SampleType, shape = lcbatch))
) |>
  wrap_plots() +
  plot_layout(guides = "collect")

# While the data before modelling is mainly driven by batch effects, the APCA clearly separates the two cell populations.
# Interestingly, the PCA on the residuals suggests that there is a small subpopulation that we did not model.
# This is further explored in: https://uclouvain-cbio.github.io/SCP.replication/articles/scplainer_leduc2022.html#downstream-analysis

# We can use the same approach to explore the component results in feature / peptide space.

caResPeps <- caRes$byFeature
caResPeps <- scpAnnotateResults(
  caResPeps, rowData(sce), by = "feature", by2 = "Sequence"
)
plCApeps <- scpComponentPlot(
  caResPeps, pointParams = list(size = 0.8, alpha = 0.4)
) |>
  wrap_plots()

# We can combine the exploration of the components in cell and peptide space using biplots.
# On the peptide level:

biplots <- scpComponentBiplot(
  caResCells, caResPeps, 
  pointParams = list(aes(colour = SampleType, shape = lcbatch)),
  labelParams = list(size = 1.5, max.overlaps = 20),
  textBy = "gene", top = 20
) |>
  wrap_plots(guides = "collect")

# On the protein level:
caResProts <- scpComponentAggregate(caResPeps, fcol = "gene")
biplots <- scpComponentBiplot(
  caResCells, caResProts, 
  pointParams = list(aes(colour = SampleType, shape = lcbatch)),
  labelParams = list(size = 1.5, max.overlaps = 20),
  textBy = "gene", top = 20
) |>
  wrap_plots(guides = "collect")

## Batch Correction

# Based on the fitted model we can generate batch-corrected data, that is data with only the effect of cell type and the residual data.
# We can also remove the intercept (?).

(batch_corrected <- scpRemoveBatchEffect(
  sce, effects = c("Set", "Channel", "MedianIntensity"),
  intercept = TRUE
))

## End Note

# Other down-stream analyses can be found here: 
# https://uclouvain-cbio.github.io/SCP.replication/articles/scplainer_leduc2022.html#downstream-analysis
#
# Other data set explorations with SCPlainer can be found here:
# https://uclouvain-cbio.github.io/SCP.replication/index.html#scp-data-re-analysis-using-scplainer
#