library(tidyverse)
library(scp)
library(scpdata)
library(patchwork)
library(ensembldb)
library(EnsDb.Hsapiens.v86)
library(bluster)
library(scater)

leduc <- leduc2022_pSCoPE()
assaysToRemove <- c("peptides", "peptides_log", "proteins_norm2", "proteins_processed")
leduc <- removeAssay(leduc, assaysToRemove)
requiredRowData <- c(
  "Sequence", "Leading.razor.protein.symbol", 
  "Leading.razor.protein.id", "Reverse", "Potential.contaminant",
  "Leading.razor.protein", "PIF", "dart_qval"
)
leduc <- selectRowData(leduc, requiredRowData)
leduc <- zeroIsNA(leduc, i = names(leduc))
leduc <- computeSCR(
  leduc, names(leduc), colvar = "SampleType", 
  samplePattern = "Mel|Macro", carrierPattern = "Carrier",
  sampleFUN = "mean", rowDataName = "MeanSCR"
)
leduc <- filterFeatures(
  leduc, ~ Reverse != "+" &
    Potential.contaminant != "+" &
    !grepl("REV|CON", Leading.razor.protein) &
    !is.na(PIF) & PIF > 0.6 &
    dart_qval < 0.01 &
    !is.na(MeanSCR) & MeanSCR < 0.05
)
leduc <- countUniqueFeatures(
  leduc, i = names(leduc), groupBy = "Sequence",
  colDataName = "NumberPeptides"
)
MedianIntensity <- lapply(experiments(leduc), function(x) {
  out <- colMedians(log(assay(x)), na.rm = TRUE)
  names(out) <- colnames(x)
  out
})
names(MedianIntensity) <- NULL
MedianIntensity <- unlist(MedianIntensity)
colData(leduc)[names(MedianIntensity), "MedianIntensity"] <- MedianIntensity
leduc <- medianCVperCell(
  leduc, i = names(leduc), groupBy = "Leading.razor.protein.symbol",
  nobs = 3, na.rm = TRUE, colDataName = "MedianCV", norm = "SCoPE2"
)
ggplot(data.frame(colData(leduc))) +
  aes(
    y = MedianIntensity,
    x = NumberPeptides,
    color = MedianCV,
    shape = SampleType
  ) +
  geom_point(size = 2) +
  scale_color_continuous(type = "viridis")
passQC <- !is.na(leduc$MedianCV) & leduc$MedianCV < 0.6 &
  leduc$MedianIntensity > 6 & leduc$MedianIntensity < 8 &
  leduc$NumberPeptides > 750 &
  grepl("Mono|Mel", leduc$SampleType)
leduc <- subsetByColData(leduc, passQC)
peptideAssays <- paste0("peptides_", names(leduc))
leduc <- aggregateFeatures(leduc,
                           i = names(leduc),
                           fcol = "Sequence",
                           name = peptideAssays,
                           fun = colMedians,
                           na.rm = TRUE)
ppMap <- rbindRowData(leduc, i = grep("^pep", names(leduc))) %>%
  data.frame %>%
  group_by(Sequence) %>%
  ## The majority vote happens here
  mutate(Leading.razor.protein.symbol =
           names(sort(table(Leading.razor.protein.symbol),
                      decreasing = TRUE))[1],
         Leading.razor.protein.id =
           names(sort(table(Leading.razor.protein.id),
                      decreasing = TRUE))[1]) %>%
  dplyr::select(Sequence, Leading.razor.protein.symbol, Leading.razor.protein.id) %>%
  dplyr::filter(!duplicated(Sequence, Leading.razor.protein.symbol))
consensus <- lapply(peptideAssays, function(i) {
  ind <- match(rowData(leduc)[[i]]$Sequence, ppMap$Sequence)
  DataFrame(Leading.razor.protein.symbol =
              ppMap$Leading.razor.protein.symbol[ind],
            Leading.razor.protein.id = 
              ppMap$Leading.razor.protein.id[ind])
})
names(consensus) <- peptideAssays
rowData(leduc) <- consensus
leduc <- joinAssays(leduc, i = peptideAssays, 
                    name = "peptides")
proteinIds <- rowData(leduc)[["peptides"]]$Leading.razor.protein.id
proteinConversionDf <- transcripts(
  EnsDb.Hsapiens.v86, 
  columns = "gene_name",
  return.type = "data.frame",
  filter = UniprotFilter(proteinIds)
)
matchedIndex <- match(proteinIds, proteinConversionDf$uniprot_id)
geneName <- proteinConversionDf$gene_name[matchedIndex]
rowData(leduc)[["peptides"]]$gene <- geneName
leduc <- logTransform(leduc, i = "peptides", name = "peptides_log")
plot(leduc)
sce <- getWithColData(leduc, "peptides_log")

readr::write_rds(sce, "sce.RDS", compress = "xz", version = 2, text = F)
