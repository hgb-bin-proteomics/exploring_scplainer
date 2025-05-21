library(SingleCellExperiment)

### Examplary reading of Spectronaut TSV result file into a SingleCellExperiment

# You can retrieve the example file from here:
# https://www.ebi.ac.uk/pride/archive/projects/PXD049412
#
# or directly using this link:
# https://ftp.pride.ebi.ac.uk/pride/data/archive/2025/01/PXD049412/JB_HeLa_10ng_Celegans_noMBR_Report_Peptides_JB_Pivot_.tsv

spectronaut <- read.table("JB_HeLa_10ng_Celegans_noMBR_Report_Peptides_JB_Pivot_.tsv",
                          header = T,
                          sep = "\t")
colnames(spectronaut)

# For the sake of simplicity, lets only take a small sample of this data:

spectronaut <- head(spectronaut, 25)

## ASSAY DATA

# We now need to create our assay data:

# First we select only the quantification columns that we are interested in,
# we can use a simple regex for that:

quant_cols <- grep(pattern = "*PEP.MS1Quantity$", colnames(spectronaut))
quant_data <- spectronaut[,quant_cols]

# We might now want to rename our columns to the sample / cell names, because
# we currently have these names:

colnames(quant_data)

# Let's name them SPD1, SPD2, and SPD3:

colnames(quant_data) <- c("SPD1", "SPD2", "SPD3")

# We also currently only have numeric values as rownames:

rownames(quant_data)[1:10]

# We need to rename this, in our case to precursor names:

rownames(quant_data) <- spectronaut$EG.PrecursorId
rownames(quant_data)[1:10]

# Note that rownames have to be unique, data needs to be correctly aggregated on PSM / peptide / precursor / protein level!

# Finally we have to coerce our data into a matrix:

assay <- as.matrix(quant_data)

## COLDATA

# Now we need to define the colData:

# Let's pretend for the sake of this example that we have different cell types, acquistion batches, and
# we will calculate the median intensity for each cell as additional information:

col_data <- data.frame(
  CellType=c("Erythrocyte", "Lymphocyte", "Endothelial"),
  Batch=c(1, 2, 3),
  MedianIntensity=apply(assay, MARGIN = 2, FUN = median, na.rm = T)
)

# Remember that our rownames need to match to the colnames of our assay, where:
# The data in row 1 corresponds to col 1, row 2 to col 2, and so on:

rownames(col_data) <- colnames(assay)

# Therefor our SPD1 cell has the following data associated:

col_data["SPD1",]

## ROWDATA

# In our rowData we can store additional information about our precursors, for example Protein group information
# and if it's gene specific. We can directly select this from our Spectronaut result:

row_data <- spectronaut[,c("PG.ProteinGroups", "PEP.IsGeneSpecific")]

# Remember that our rownames need to match the rownames of our assay:

rownames(row_data) <- rownames(assay)

## SINGLECELLEXPERIMENT

# Finally we can create our SingleCellExperiment:

sce <- SingleCellExperiment(
  list(quant=assay), 
  colData = col_data, 
  rowData = row_data
)

# Inspecting our create SingleCellExperiment:

sce
sce$CellType
sce$Batch
sce$MedianIntensity
assay(sce, "quant")
colData(sce)
rowData(sce)
