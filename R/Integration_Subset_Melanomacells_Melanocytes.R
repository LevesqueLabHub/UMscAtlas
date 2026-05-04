
#  This script is for integrating subsets of melanoma cells and melanocytes. Harmony, dims = 50, theta = 2, clustering res = 0.5

# Libraries

library(Seurat)
library(harmony)
library(ggplot2)
library(qs)
library(patchwork)
library(edgeR)
library(randomcoloR)
library(Polychrome)
library(pals)
library(dplyr)


# Load data
sce <- qread("../13_Melanomasubset_finalData.qs")


# Adjustments
sce$samplename <- as.character(sce$samplename)

new_names <- c(
  "pUM09"       = "pUM09_1",
  "pUM09_16S"   = "pUM09_2",
  "pUM11"       = "pUM11_1",
  "pUM11_16S"   = "pUM11_2",
  "pUM12"       = "pUM12_1",
  "pUM12_16S"   = "pUM12_2",
  "pUM15"       = "pUM15_1",
  "pUM15_16S"   = "pUM15_2",
  "pUM17"       = "pUM17_1",
  "pUM17_16S"   = "pUM17_2",
  "mUM13_FNA"   = "mUM13_2",
  "mUM13_CNA"   = "mUM13_1"
)


# adjust
sce$samplename <- recode(
  as.character(sce$samplename),
  !!!new_names
)

table(sce$samplename)

# delete 2 cells

sce <- subset(
  sce,
  subset = !(samplename == "mUM15" & majority_celltype == "Melanocytes")
)



# Run Harmony integration
melanoma_integrated <- RunHarmony(
  object = sce,
  group.by.vars = c("Condition", "orig.ident"),
  reduction.use = "pca",
  dims.use = 1:50,
  assay.use = "RNA",
  theta = c(2, 2)
)


# Run UMAP
melanoma_integrated <- RunUMAP(
  object = melanoma_integrated,
  reduction = "harmony",
  dims = 1:50
)


# Clustering
melanoma_integrated <- FindNeighbors(object = melanoma_integrated, reduction = "harmony", dims = 1:50)
melanoma_integrated <- FindClusters(object = melanoma_integrated, resolution  = c(0.2, 0.3, 0.4,0.5), verbose = FALSE)


# Save Melanoma / Melanocyte subset
qsave(melanoma_integrated, "../01_integrated_melanoma_melanocyte_subset.qs")

# Dimplot:
DimPlot(melanoma_integrated, reduction = "umap",
        group.by = "seurat_clusters",
        label = TRUE,
        repel = TRUE)

DimPlot(melanoma_integrated, reduction = "umap",
        group.by = "Condition",
        split.by = "Condition")

DimPlot(melanoma_integrated, reduction = "umap",
        group.by = "Condition")

DimPlot(melanoma_integrated, reduction = "umap",
        group.by = "orig.ident")







