#!/usr/bin/env Rscript
# CytoTRACE2 + stemness-correlated genes on the October 2026 atlas (349,696 cells).
# 2026-10-05, meeting follow-up. Old scores (UM_atlas_cytotrace2.qs2, January object) cannot be
# transferred: counts changed in every library. Writes only small tables, read by 01b_CytoTRACE2_2026b.Rmd.
#   cytotrace2_scores.qs2        one row per cell (cell, Sample, CytoTRACE2_* columns)
#   gene_cor_<set>.csv           Pearson r of CytoTRACE2_Score with every gene, per cell set
# CytoTRACE2 is run per Sample, so its KNN smoothing only ever sees cells of one library.
suppressPackageStartupMessages({ library(Seurat); library(CytoTRACE2); library(Matrix); library(data.table) })

OBJ <- "../new_data_2026b/Objects/full Dataset/final_Subtyping_iCNV_reclustered.qs2"
OUT <- "../Outputs/CytoTRACE2_2026b"
NCORES <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "8"))
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

obj <- qs2::qs_read(OBJ, nthreads = 1)
stopifnot(ncol(obj) == 349696)
obj[["RNA"]] <- JoinLayers(obj[["RNA"]])
counts <- LayerData(obj, assay = "RNA", layer = "counts")
meta   <- obj@meta.data[, c("Sample", "origin", "Malignant", "majority_celltype")]
rm(obj); invisible(gc())

score_file <- file.path(OUT, "cytotrace2_scores.qs2")
if (!file.exists(score_file)) {
  # CytoTRACE2's KNN smoothing splits cells over ncores chunks; on a small library some chunks come out
  # empty and it dies with "cut(): 'x' must be numeric" (reproduced: 290 cells fail at 32 cores, pass at 1).
  # So use about one core per 500 cells. Per-library caches let a rerun resume.
  lib_dir <- file.path(OUT, "per_library"); dir.create(lib_dir, showWarnings = FALSE)
  res <- rbindlist(lapply(split(colnames(counts), meta$Sample), function(cells) {
    f <- file.path(lib_dir, paste0(meta[cells[1], "Sample"], ".qs2"))
    if (file.exists(f)) return(qs2::qs_read(f, nthreads = 1))
    r <- cytotrace2(as.matrix(counts[, cells, drop = FALSE]), is_seurat = FALSE, species = "human",
                    ncores = min(NCORES, max(1L, length(cells) %/% 500L)), batch_size = 10000, seed = 14)
    r <- data.table(cell = rownames(r), r)
    qs2::qs_save(r, f, nthreads = 1)
    r
  }))
  res[, Sample := meta[cell, "Sample"]]
  stopifnot(nrow(res) == ncol(counts), !anyNA(res$CytoTRACE2_Score), setequal(res$cell, colnames(counts)))
  qs2::qs_save(res, score_file, nthreads = 1)
}
res <- qs2::qs_read(score_file, nthreads = 1)
score <- setNames(res$CytoTRACE2_Score, res$cell)[colnames(counts)]

# Pearson r of score vs log1p(CP10k) per gene, from sparse sums (no dense genes x cells matrix)
norm <- counts %*% Diagonal(x = 1e4 / colSums(counts)); norm@x <- log1p(norm@x)
dimnames(norm) <- dimnames(counts)  # %*% Diagonal() drops the cell names
gene_cor <- function(cells) {
  X <- norm[, cells, drop = FALSE]; y <- score[cells]; n <- length(y)
  sx <- rowSums(X); sxx <- rowSums(X^2); sxy <- as.numeric(X %*% y)
  r <- (n * sxy - sx * sum(y)) / sqrt((n * sxx - sx^2) * (n * sum(y^2) - sum(y)^2))
  data.table(gene = rownames(X), correlation = r, pct_expressed = rowSums(X > 0) / n)[
    pct_expressed >= 0.01][order(-correlation)]
}
sets <- list(
  malignant          = rownames(meta)[meta$Malignant == "Malignant"],
  healthy_melanocyte = rownames(meta)[meta$origin == "non-malignant" & meta$majority_celltype == "Melanocytes"],
  all_cells          = rownames(meta))
stopifnot(lengths(sets) > 1000)
for (s in names(sets)) {
  g <- gene_cor(sets[[s]])
  # guard: the sparse formula must equal cor() on one gene
  top <- g$gene[1]
  stopifnot(abs(g$correlation[1] - cor(as.numeric(norm[top, sets[[s]]]), score[sets[[s]]])) < 1e-8)
  fwrite(g, file.path(OUT, sprintf("gene_cor_%s.csv", s)))
}
cat("done:", nrow(res), "cells scored\n")
