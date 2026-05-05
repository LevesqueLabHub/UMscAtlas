
# This script is used to identify genes, in the melanoma cells that differ significantly inside the locations in the eye (iris, ciliary body, choroid) using DeSeq2, DEGpatterns, GO, KEGG


# Libraries
library(Seurat)
library(ggplot2)
library(DESeq2)
library(qs)
library(patchwork)
library(edgeR)
library(dplyr)
library(tibble)
library(randomcoloR)
library(Polychrome)
library(pals)
library(DEGreport)
library(org.Hs.eg.db)
library(clusterProfiler)
library(pheatmap)
library(enrichplot)


# Load and subset data
sce     <- qread("../../../../data/01_Subsets/02_melanoma_cells/integrated/01_integrated_melanoma_subset.qs")
sce_eye <- subset(sce, subset = location %in% c("choroid", "ciliary body", "iris"))


# Pseudobulk
pb <- AggregateExpression(sce_eye,
                          group.by = c("orig.ident", "location"),
                          assays   = "RNA",
                          slot     = "counts")$RNA

meta_data <- sce_eye@meta.data %>%
  distinct(orig.ident, location, Condition)
rownames(meta_data) <- paste0(meta_data$orig.ident, "_", meta_data$location)

stopifnot(all(colnames(pb) == rownames(meta_data)))


# Cell counts per patient
cell_counts <- sce_eye@meta.data %>%
  dplyr::count(orig.ident, location)
print(cell_counts)


# DeSeq2 object
dds <- DESeqDataSetFromMatrix(
  countData = pb,
  colData   = meta_data,
  design    = ~ Condition + location
)

smallestGroupSize <- 3
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds  <- dds[keep, ]


# Run DeSeq2

dds <- DESeq(dds)


# Pairwise comparison
res_chor_cil  <- results(dds, contrast = c("location", "choroid",      "ciliary body"))
res_cil_iris  <- results(dds, contrast = c("location", "ciliary body", "iris"))
res_chor_iris <- results(dds, contrast = c("location", "choroid",      "iris"))


# Significant DE genes (padj < 0.05, |LFC| > 1)
sig_genes <- function(res) {
  rownames(res[which(res$padj < 0.05 & abs(res$log2FoldChange) > 1), ])
}

de_genes_union <- unique(c(
  sig_genes(res_chor_cil),
  sig_genes(res_cil_iris),
  sig_genes(res_chor_iris)
))
cat("Number of DE genes:", length(de_genes_union), "\n")


# DEGPatterns
rld <- rlog(dds, blind = FALSE)
mat <- assay(rld)[de_genes_union, ]

location_colors <- c("iris" = "#00778BFF", "choroid" = "#211747FF", "ciliary body" = "#B6E2E0FF")

patterns <- degPatterns(mat,
                        metadata = colData(dds),
                        time     = "location",
                        col      = NULL,
                        minc     = 10)


# Custom DEGpatterns plot
metadata      <- as.data.frame(colData(dds))
metadata$loc  <- metadata$location

newtry <- degPatterns(mat,
                      metadata = metadata,
                      time     = "location",
                      col      = "loc",
                      minc     = 10)

newtry$normalized$location_num <- as.numeric(factor(newtry$normalized$location,
                                                    levels = c("choroid", "ciliary body", "iris")))
newtry$summarise$location_num  <- as.numeric(factor(newtry$summarise$location,
                                                    levels = c("choroid", "ciliary body", "iris")))

cluster_labels <- c("1" = "Group 1: 82 genes",
                    "2" = "Group 2: 36 genes",
                    "3" = "Group 3: 1511 genes",
                    "4" = "Group 4: 1048 genes")

ggplot(newtry$normalized, aes(x = location_num, y = value, color = loc, fill = loc)) +
  geom_boxplot(aes(group = interaction(location_num, loc)), alpha = 0.3, outlier.shape = NA) +
  geom_line(aes(group = genes), color = "black", linewidth = 0.1, alpha = 0.1) +
  geom_point(position = position_jitterdodge(), alpha = 0.3, size = 0.5) +
  scale_x_continuous(breaks = 1:3, labels = c("choroid", "ciliary body", "iris")) +
  scale_color_manual(values = location_colors) +
  scale_fill_manual(values = location_colors) +
  facet_wrap(~cluster, labeller = labeller(cluster = cluster_labels)) +
  theme_bw() +
  labs(x = "location", y = "value")



# Save DeSeq2 object and matrix
saveRDS(dds, "../DeSeq2_location_melanoma_cells.RDS")
write.csv(mat, "../Matrix_DeSeq2_location_melanoma_cells.csv")

# Heatmap
cluster_info <- patterns$df
print(table(cluster_info$cluster))

colnames(mat) <- sub("_.*", "", colnames(mat))

annotation_col <- data.frame(
  location  = colData(dds)$location,
  row.names = colnames(mat)
)
annotation_colors <- list(location = location_colors)

cat("Number of genes in heatmap:", nrow(mat), "\n")

pheatmap(mat, cluster_rows = TRUE,  cluster_cols = TRUE,  show_rownames = FALSE,
         annotation_col = annotation_col, scale = "row", main = "Clusters across Locations",
         annotation_colors = annotation_colors, angle_col = 90)


# Sorted Heatmap by location
sample_order          <- order(colData(dds)$location)
mat_sorted            <- mat[, sample_order]
annotation_col_sorted <- annotation_col[sample_order, , drop = FALSE]

pheatmap(mat_sorted, cluster_rows = TRUE, cluster_cols = FALSE, show_rownames = FALSE,
         annotation_col = annotation_col_sorted, scale = "row", main = "Sorted by Location",
         annotation_colors = annotation_colors, angle_col = 90)


# Extract cluster gene list
cluster_genes <- patterns$df
print(table(cluster_genes$cluster))

n_clusters <- sort(unique(cluster_genes$cluster))

genes_per_cluster <- lapply(n_clusters, function(cl) {
  cluster_genes$genes[cluster_genes$cluster == cl]
})
names(genes_per_cluster) <- paste0("cluster", n_clusters)

# cluster 3 = iris-high, cluster 4 = iris-low
cluster3_genes <- genes_per_cluster$cluster3
cluster4_genes <- genes_per_cluster$cluster4


# GO
GO_function <- function(cluster_genes, cluster_name) {
  ego <- enrichGO(gene          = cluster_genes,
                  OrgDb         = org.Hs.eg.db,
                  keyType       = "SYMBOL",
                  ont           = "BP",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.05)
  
  print(dotplot(ego, showCategory = 15, font.size = 10,
                title = paste("GO Enrichment - Group", cluster_name)))
  print(barplot(ego, showCategory = 15, font.size = 10,
                title = paste("Top GO Terms - Group", cluster_name)))
  
  top20 <- ego@result[order(ego@result$pvalue), ][1:20, ]
  print(
    ggplot(top20, aes(x = Count, y = reorder(Description, -pvalue))) +
      geom_point(aes(size = Count, color = pvalue)) +
      scale_color_gradient(low = "red", high = "blue") +
      labs(x = "Gene Count", y = NULL, color = "p-value",
           title = paste("Top GO Terms - Group", cluster_name)) +
      theme_minimal() +
      theme(axis.text.y = element_text(size = 9))
  )
  
  print(
    cnetplot(ego, showCategory = 10, colorEdge = TRUE, circular = FALSE,
             node_label = "all", cex_label_category = 0.8, cex_label_gene = 0.6) +
      ggtitle(paste("Gene-Concept Network - Group", cluster_name))
  )
}


#KEGG

KEGG_function <- function(cluster_genes, cluster_name) {
  gene_entrez <- bitr(cluster_genes,
                      fromType = "SYMBOL",
                      toType   = "ENTREZID",
                      OrgDb    = org.Hs.eg.db)
  
  kegg <- enrichKEGG(gene          = gene_entrez$ENTREZID,
                     organism      = "hsa",
                     pvalueCutoff  = 0.05,
                     pAdjustMethod = "BH")
  
  print(dotplot(kegg, showCategory = 15, font.size = 10,
                title = paste("KEGG Pathways - Group", cluster_name)))
  print(barplot(kegg, showCategory = 15, font.size = 10,
                title = paste("KEGG Pathways - Group", cluster_name)))
  
  top20_kegg <- kegg@result[order(kegg@result$pvalue), ][1:min(20, nrow(kegg@result)), ]
  print(
    ggplot(top20_kegg, aes(x = Count, y = reorder(Description, -pvalue))) +
      geom_point(aes(size = Count, color = pvalue)) +
      scale_color_gradient(low = "red", high = "blue") +
      labs(x = "Gene Count", y = NULL, color = "p-value",
           title = paste("Top KEGG Pathways - Group", cluster_name)) +
      theme_minimal() +
      theme(axis.text.y = element_text(size = 9))
  )
}


# Run GO and KEGG
GO_function(cluster3_genes,   cluster_name = "iris-high 3")
GO_function(cluster4_genes,   cluster_name = "iris-low 4")

KEGG_function(cluster3_genes, cluster_name = "iris-high 3")
KEGG_function(cluster4_genes, cluster_name = "iris-low 4")








