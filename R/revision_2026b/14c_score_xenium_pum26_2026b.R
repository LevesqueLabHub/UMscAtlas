#!/usr/bin/env Rscript
# 14c scoring step: PLAID scores of the October 2026 stemness signatures on pUM26 (Xenium), primary eye
# tumour (Region5_0038498) and liver metastasis (Region1_0038471), from the SPLIT-purified objects with the
# refined RCTD labels (Xenium_pUM26/SPLIT_Reannotation). iCNV per cell from 14c_extract_pum26_icnv_obs.py.
# Same gene sets and scorer as 14b. Writes the per-cell table 14c_Xenium_Stemness_pUM26_2026b.Rmd reads.
# Run: Rscript 14c_score_xenium_pum26_2026b.R   (~2 min, ~25 GB)
suppressPackageStartupMessages({ library(Seurat); library(qs2); library(data.table); library(Matrix) })
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))))
source(file.path(here, "plaid_score.R"))
source(file.path(here, "14b_stemness_common.R"))

work    <- "../Xenium_pUM26/SPLIT_Reannotation"
out_dir <- "../Outputs/Xenium_Stemness_pUM26_2026b"
pairs <- fread(file.path(here, "../../Xenium_pUM26/SPLIT_Reannotation/pairs.tsv"))[patient == "pUM26"]
stopifnot(nrow(pairs) == 2)

objs <- lapply(pairs$xenium_region, function(r) qs_read(file.path(work, "regions", r, "scData_SPLIT_immune.qs2"), nthreads = 1))
counts <- lapply(objs, function(o) SeuratObject::LayerData(o, assay = "Xenium", layer = "counts"))
stopifnot(identical(rownames(counts[[1]]), rownames(counts[[2]])))
counts <- do.call(cbind, counts)   # pooled: z-scores relative to all eye + liver cells of pUM26

G <- stemness_sets(rownames(counts))
print(lengths(G$sets))
fwrite(G$coverage, file.path(out_dir, "panel_coverage.tsv"), sep = "\t")
S <- plaid_score(counts, G$sets)

cells <- rbindlist(lapply(seq_along(objs), function(i) {
  md <- objs[[i]]@meta.data
  data.table(cell = colnames(objs[[i]]), sample = pairs$pair[i], region = pairs$xenium_region[i], tissue = pairs$tissue[i],
             refined_immune_type = md$refined_immune_type, rctd_after = as.character(md$rctd_after),
             nFeature = md$nFeature_Xenium, nCount = md$nCount_Xenium)
}))
stopifnot(all(cells$cell == rownames(S)))
cells <- cbind(cells, as.data.table(S))

icnv <- rbindlist(lapply(pairs$pair, function(p) fread(file.path(out_dir, paste0(p, "_icnv_obs.csv")))))
stopifnot(!anyDuplicated(icnv$barcode))
cells[icnv, on = c(cell = "barcode"), `:=`(chr3_cnv = i.chr3_cnv, chr8q_cnv = i.chr8q_cnv, cnv_rms = i.cnv_rms_score,
                                          x = i.x_centroid, y = i.y_centroid, icnv_cell_type = i.cell_type)]
cells[, cnv_directed := chr8q_cnv - chr3_cnv]   # UM hallmark: chr8q gain minus chr3 loss
cat("cells:", nrow(cells), " with iCNV + coordinates:", sum(!is.na(cells$cnv_directed) & !is.na(cells$x)), "\n")
stopifnot(sum(!is.na(cells$cnv_directed)) == nrow(cells), all(cells$icnv_cell_type == cells$rctd_after))

qs_save(c(list(cells = cells, pairs = pairs, date = Sys.time()), G), file.path(out_dir, "pum26_scores_2026b.qs2"), nthreads = 1)
cat("Saved", file.path(out_dir, "pum26_scores_2026b.qs2"), "\n")
