
# # Cell type prediction and stemness scoring for external UM datasets
# using Seurat label transfer and UCell
# 
# Libraries
library(Seurat)
library(qs)
library(ggplot2)
library(patchwork)
library(clustree)
library(UCell)
library(dplyr)


# Load data
sce_Pandiani <- qread("../03_dataset_merged_seurat_Clustering.qs")
sce_Murali   <- qread("../03_dataset_merged_seurat_Clustering.qs")
sce_scAtlas  <- qread("../13_final_Dataset_reclustered_no_MuM15.qs")


# UMAP
sce_Murali   <- RunUMAP(sce_Murali,   dims = 1:26, seed.use = 42, return.model = TRUE)
sce_Pandiani <- RunUMAP(sce_Pandiani, dims = 1:26, seed.use = 42, return.model = TRUE)
sce_scAtlas  <- RunUMAP(sce_scAtlas,  dims = 1:50, seed.use = 42, return.model = TRUE)


# Cell type prediction
# Finds anchors between reference (scAtlas) and query dataset,
# then projects query onto reference UMAP space

predicting_celltypes <- function(dataset, name) {
  
  anchors <- FindTransferAnchors(
    reference           = sce_scAtlas,
    query               = dataset,
    dims                = 1:26,
    reference.reduction = "pca"
  )
  
  query <- MapQuery(
    anchorset          = anchors,
    reference          = sce_scAtlas,
    query              = dataset,
    refdata            = list(celltype = "majority_celltype"),
    reference.reduction = "pca",
    reduction.model    = "umap"
  )
  
  # Plot 1: reference + projected query side by side
  p1 <- DimPlot(sce_scAtlas, reduction = "umap", group.by = "majority_celltype",
                label = TRUE, label.size = 3, repel = TRUE) +
    ggtitle("Reference celltypes") + coord_fixed()
  
  p2 <- DimPlot(query, reduction = "ref.umap",
                label = TRUE, label.size = 3, repel = TRUE) +
    ggtitle(paste0(name, " and scAtlas")) + coord_fixed()
  
  # Plot 2: query projected onto reference (reference in grey)
  ref_coords   <- as.data.frame(Embeddings(sce_scAtlas, "umap"))
  ref_coords$UMAP_1 <- ref_coords[, 1]
  ref_coords$UMAP_2 <- ref_coords[, 2]
  
  query_coords          <- as.data.frame(Embeddings(query, "ref.umap"))
  query_coords$UMAP_1   <- query_coords[, 1]
  query_coords$UMAP_2   <- query_coords[, 2]
  query_coords$celltype <- query$predicted.celltype
  
  p3 <- ggplot() +
    geom_point(data = ref_coords,
               aes(x = UMAP_1, y = UMAP_2),
               color = "grey90", size = 0.1, alpha = 0.3) +
    geom_point(data = query_coords,
               aes(x = UMAP_1, y = UMAP_2, color = celltype),
               size = 0.5, alpha = 0.7) +
    theme_minimal() +
    coord_fixed() +
    labs(title  = paste0(name, " Projected onto scAtlas (grey)"),
         color  = "Predicted Cell Type")
  
  print(p1 + p2)
  print(p3)
  
  return(query)
}


# Run cell type prediction
Murali   <- predicting_celltypes(sce_Murali,   "Murali")
Pandiani <- predicting_celltypes(sce_Pandiani, "Pandiani")

# Save
qsave(Murali,   "../04_predicted_celltypes.qs")
qsave(Pandiani, "../04_predicted_celltypes_Pandiani.qs")


# Mean prediction score per cell type
aggregate(predicted.celltype.score ~ predicted.celltype, data = Murali@meta.data,   FUN = mean)
aggregate(predicted.celltype.score ~ predicted.celltype, data = Pandiani@meta.data, FUN = mean)

# Stemness score function
stemness_genes <- c(
  "HSPD1", "PABPC1", "RPLP0",  "HNRNPA1", "RPS19",  "RPL22",  "RPL36A", "NPM1",
  "HSP90AB1", "MYC",  "NAP1L1", "RPLP2",   "DNAJA1", "RPS16",  "RPL11",  "RPL4",
  "RPL39",  "RPS21", "YBX1",   "RPL27A",  "ID2",    "TAF1D",  "RPL6",   "HSPE1",
  "RPL23",  "EIF5",  "HSPH1",  "RPL18",   "RPL22L1","JUN",    "EEF1G",  "UQCRH",
  "H3-3A",  "UBC",   "SRSF3",  "RPL31",   "RPS2",   "JUND",   "RPL7",   "CCT2",
  "RPL18A", "EIF3E", "EEF1D",  "FOS",     "RPL38",  "HSPA8",  "RPL5",   "YBX3",
  "RPL37",  "RAN"
)

# Check gene availability
cat("Stemness genes in Murali:\n");   print(stemness_genes %in% rownames(sce_Murali))
cat("Stemness genes in Pandiani:\n"); print(stemness_genes %in% rownames(sce_Pandiani))

majority_celltype_col <- c(
  "Glial cells"         = "#5A5156",
  "B/Plasma cells"      = "lightslateblue",
  "Myeloid cells"       = "lightblue",
  "T/NK cells"          = "#325A9B",
  "Melanoma cells"      = "#B00068",
  "Melanocytes"         = "#C075A6",
  "Endothelial cells"   = "chartreuse4",
  "Epithelial cells"    = "#1CBE4F",
  "Fibroblasts"         = "yellowgreen",
  "Mesenchymal cells"   = "darkgreen",
  "Pericytes"           = "olivedrab4",
  "Smooth muscle cells" = "#1CFFCE",
  "Photoreceptor cells" = "peachpuff3"
)

stemness_score <- function(sce, name) {
  
  seurat_obj <- AddModuleScore_UCell(
    sce,
    features = list(Stemness_UM_UCell = stemness_genes),
    ncores   = 8
  )
  
  # Order cell types by median stemness score (descending)
  order <- seurat_obj@meta.data %>%
    group_by(predicted.celltype) %>%
    summarise(median = median(Stemness_UM_UCell_UCell, na.rm = TRUE)) %>%
    arrange(desc(median)) %>%
    pull(predicted.celltype)
  
  seurat_obj$predicted.celltype <- factor(seurat_obj$predicted.celltype, levels = order)
  
  p <- VlnPlot(seurat_obj,
               features  = "Stemness_UM_UCell_UCell",
               group.by  = "predicted.celltype",
               pt.size   = 0,
               alpha     = 0.1,
               cols      = majority_celltype_col) +
    ggtitle(paste0(name, " - Stemness Score")) +
    ylim(0, 1)
  print(p)
  
  return(seurat_obj)
}

# Run stemness scoring 
stemness_Pandiani <- stemness_score(Pandiani, "Pandiani")
stemness_Murali   <- stemness_score(Murali,   "Murali")




