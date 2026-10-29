#!/usr/bin/env Rscript
# 02b: PLAID scores of the 67 Weizmann 3CA malignant metaprograms on the October 2026 atlas (Paul, 2026-10-13),
# for 02b_Cancer_Metaprograms_2026b.Rmd (re-run of report 02, Victoria: "repeat Weizman metaprograms analysis with
# new samples, on melanoma cells and melanocytes").
# Cells: majority_celltype "Melanoma cells" (malignant) + "Melanocytes" (healthy). Both groups are scored in ONE
# plaid_score() call on raw counts, so every score is a z-score relative to the same melanocytic pool.
# Outputs (Outputs/Metaprograms_3CA_2026b/): per-cell scores + metadata, per-set gene coverage. Stats are in the report.
# Run: Rscript 02b_score_3CA_2026b.R   (~30 GB object, read at nthreads = 1)
suppressPackageStartupMessages({ library(Seurat); library(data.table) })
source("plaid_score.R")
VA  <- ".."
OBJ <- file.path(VA, "new_data_2026b/Objects/full Dataset/final_Subtyping_iCNV_reclustered.qs2")
OUT <- file.path(VA, "Outputs/Metaprograms_3CA_2026b"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
CELLTYPES <- c("Melanoma cells", "Melanocytes")

mp_df <- readxl::read_excel("../reference/Weizman_3CA_meta_programs_2025-01-29.xlsx", sheet = "Malignant")
sets <- lapply(mp_df, function(g) unique(g[!is.na(g)]))
names(sets) <- stringr::str_squish(names(sets))          # "MP38  NPC (glioma)" has a double space
stopifnot(length(sets) == 67, all(grepl("^MP[0-9]+ ", names(sets))))

t0 <- Sys.time()
obj <- qs2::qs_read(OBJ, nthreads = 1)
cat("object read in", format(Sys.time() - t0), ":", ncol(obj), "cells\n")
stopifnot(ncol(obj) == 349696, all(CELLTYPES %in% obj$majority_celltype))
keep <- colnames(obj)[obj$majority_celltype %in% CELLTYPES]
meta <- obj@meta.data[keep, c("orig.ident", "Sample", "Patient_nr", "origin", "majority_celltype", "Malignant",
                              "primary_mutation", "secondary_mutation", "organ", "primary_location",
                              "healthy_location", "metastatic_location", "TreatmentStage_ofProcessedSample")]
obj[["RNA"]] <- JoinLayers(obj[["RNA"]])
cnt <- LayerData(obj, assay = "RNA", layer = "counts")[, keep]   # dgCMatrix subset, not subset() on the Assay5
rm(obj); invisible(gc())
stopifnot(inherits(cnt, "dgCMatrix"), ncol(cnt) == length(keep), all(colnames(cnt) == keep))
cat("counts:", nrow(cnt), "genes x", ncol(cnt), "cells; integer-valued:", all(cnt@x == round(cnt@x)),
    "; max", max(cnt@x), "\n")

cover <- data.table(MP = names(sets), n_genes = lengths(sets),
                    n_found = vapply(sets, function(g) sum(g %in% rownames(cnt)), integer(1)))
cover[, missing := vapply(sets, function(g) paste(setdiff(g, rownames(cnt)), collapse = ","), character(1))]
fwrite(cover, file.path(OUT, "mp_gene_coverage.csv"))
print(cover[n_found < n_genes, .(MP, n_genes, n_found)])

t1 <- Sys.time()
sc <- plaid_score(cnt, sets)
cat("PLAID scored in", format(Sys.time() - t1), "\n")
stopifnot(nrow(sc) == length(keep), all(rownames(sc) == keep), !anyNA(sc), all(colnames(sc) == names(sets)))

cells <- data.table(cell = keep, meta, sc, check.names = FALSE)
qs2::qs_save(cells, file.path(OUT, "plaid_3CA_scores_per_cell.qs2"), nthreads = 1)
cat("done:", nrow(cells), "cells,", ncol(sc), "metaprograms,", uniqueN(cells$Patient_nr), "patients\n")
print(table(cells$majority_celltype, cells$origin))
