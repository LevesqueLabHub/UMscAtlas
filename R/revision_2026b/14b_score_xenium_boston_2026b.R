#!/usr/bin/env Rscript
# 14b scoring step: PLAID scores of the October 2026 stemness signatures (01b) on the Boston/Sullivan
# Xenium cohort, plus the January signature (report 14), cell cycle and heat-shock / immediate-early sets.
# Writes one per-cell table (coordinates, cell type, CNV, scores) that 14b_Xenium_Stemness_Boston_2026b.Rmd reads.
# Run: Rscript 14b_score_xenium_boston_2026b.R   (~2 min, ~40 GB)
suppressPackageStartupMessages({ library(Seurat); library(qs2); library(data.table); library(Matrix) })
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))))
source(file.path(here, "plaid_score.R"))
source(file.path(here, "14b_stemness_common.R"))

base    <- ".."
out_dir <- file.path(base, "Outputs/Xenium_Stemness_Boston_2026b")
dir.create(file.path(out_dir, "inputs"), recursive = TRUE, showWarnings = FALSE)

# Same cells and counts as Boston cohort merged_xenium_all_samples.qs2 (layers joined), plus CellType_Display and the
# January UCell scores of report 14.
obj_path <- "../Outputs/Xenium_Stemness/xenium_Boston_Sullivan_with_scores.qs2"
cnv_path <- "../Xenium_Boston_Sullivan/06_InSituCNV_Outputs/cnv_scores_09method_all_samples.csv"

# January signature (report 14 inputs; untracked files in the main checkout) - keep a copy next to the cache.
jan_src <- "../Outputs/Xenium_Stemness"
for (f in c("cytotrace2_top_genes_xenium.csv", "cytotrace2_bottom_genes_xenium.csv"))
  if (!file.exists(file.path(out_dir, "inputs", f))) file.copy(file.path(jan_src, f), file.path(out_dir, "inputs", f))

cat("Loading", obj_path, "\n")
obj <- qs_read(obj_path, nthreads = 1)
counts <- SeuratObject::LayerData(obj, assay = "Xenium", layer = "counts")
stopifnot(inherits(counts, "dgCMatrix"), ncol(counts) == ncol(obj),
          all(abs(Matrix::colSums(counts) - obj$nCount_Xenium) < 1e-6))   # counts are the raw Xenium counts

G <- stemness_sets(rownames(counts))
if (!file.exists(panel_sig_file)) fwrite(G$panel_sig, panel_sig_file)   # gene,correlation,direction
print(lengths(G$sets))
fwrite(G$coverage, file.path(out_dir, "panel_coverage.tsv"), sep = "\t")

S <- plaid_score(counts, G$sets)   # all sets in ONE call: z-scores relative to all 248k Boston cells

# Spatial centroids, one FOV per sample
coords <- rbindlist(lapply(names(obj@images), function(f) {
  d <- as.data.table(GetTissueCoordinates(obj, image = f, which = "centroids")); d[, fov := f]; d
}))
stopifnot(!anyDuplicated(coords$cell))

md <- obj@meta.data
cells <- data.table(cell = colnames(obj), sample = md$Sample, tissue = md$Tissue_Type,
                    timepoint = md$Timepoint, response = md$Best_Response,
                    celltype = as.character(md$CellType_Display), nFeature = md$nFeature_Xenium,
                    nCount = md$nCount_Xenium, jan_UCell_Top = md$CytoTRACE2_Top,
                    jan_UCell_Bottom = md$CytoTRACE2_Bottom, Cycling_UCell = md$Cycling)
cells <- cbind(cells, as.data.table(S))
cells[coords, on = "cell", `:=`(x = i.x, y = i.y, fov = i.fov)]

cnv <- fread(cnv_path)
cnv[, cell := paste0(sample, "_", cell_barcode)]
stopifnot(!anyDuplicated(cnv$cell))
cells[cnv, on = "cell", cnv_score := i.cnv_score]
cat("cells:", nrow(cells), " with coordinates:", sum(!is.na(cells$x)), " with CNV:", sum(!is.na(cells$cnv_score)), "\n")
stopifnot(sum(!is.na(cells$cnv_score)) > 0.9 * nrow(cells))   # join worked

qs_save(c(list(cells = cells, obj_path = obj_path, cnv_path = cnv_path, date = Sys.time()), G),
        file.path(out_dir, "boston_scores_2026b.qs2"), nthreads = 1)
cat("Saved", file.path(out_dir, "boston_scores_2026b.qs2"), "\n")
