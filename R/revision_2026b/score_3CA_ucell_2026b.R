#!/usr/bin/env Rscript
# UCell scores of the 67 Weizmann 3CA malignant metaprograms on the October 2026 atlas (Paul, 2026-10-15), for the
# manuscript panels (manuscript_panels_2026b_cytotrace_mp.R). The submitted Fig. 3 panels were UCell, so their October
# versions are UCell too; the 02b report uses PLAID (mean z), which must never be labelled UCell.
# Cells and gene-set parsing as 02b_score_3CA_2026b.R: majority_celltype "Melanoma cells" + "Melanocytes", sheet
# "Malignant", names squished. UCell settings as old 02 (AddModuleScore_UCell defaults: maxRank 1500, w_neg 1,
# ties "average", chunk.size 100). UCell ranks genes within each cell, so raw counts and log-normalised data give the
# same score (old 02 ran on the "data" layer).
# Outputs (Outputs/Metaprograms_3CA_2026b/): ucell_3CA_scores_per_cell.qs2 (cell + metadata + 67 scores) and
# ucell_vs_plaid_spearman.csv (per-program cell-level Spearman rho against the PLAID cache; flag rho < 0.5).
# Run: Rscript score_3CA_ucell_2026b.R   (~30 GB object, read at nthreads = 1)
suppressPackageStartupMessages({ library(Seurat); library(data.table) })
VA  <- ".."
OBJ <- file.path(VA, "new_data_2026b/Objects/full Dataset/final_Subtyping_iCNV_reclustered.qs2")
OUT <- file.path(VA, "Outputs/Metaprograms_3CA_2026b")
CELLTYPES <- c("Melanoma cells", "Melanocytes")
NCORES <- 32
RHO_FLAG <- 0.5

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
stopifnot(inherits(cnt, "dgCMatrix"), ncol(cnt) == length(keep), all(colnames(cnt) == keep), length(keep) == 200634)

t1 <- Sys.time()
sc <- UCell::ScoreSignatures_UCell(cnt, features = sets, name = "", ncores = NCORES)
cat("UCell scored in", format(Sys.time() - t1), "\n")
sc <- sc[keep, names(sets), drop = FALSE]
stopifnot(nrow(sc) == length(keep), !anyNA(sc), all(sc >= 0 & sc <= 1))

cells <- data.table(cell = keep, meta, sc, check.names = FALSE)
qs2::qs_save(cells, file.path(OUT, "ucell_3CA_scores_per_cell.qs2"), nthreads = 1)
cat("done:", nrow(cells), "cells,", ncol(sc), "metaprograms,", uniqueN(cells$Patient_nr), "patients\n")

# Agreement with PLAID, cell level, same cells.
plaid <- qs2::qs_read(file.path(OUT, "plaid_3CA_scores_per_cell.qs2"), nthreads = 1)
stopifnot("PLAID cache has other cells" = nrow(plaid) == nrow(cells) && all(plaid$cell == cells$cell))
agree <- rbindlist(lapply(names(sets), function(m) data.table(
  MP = m, spearman_rho = cor(cells[[m]], plaid[[m]], method = "spearman"),
  ucell_mean = mean(cells[[m]]), ucell_frac_zero = mean(cells[[m]] == 0))))
agree[, flag_rho_below_0.5 := spearman_rho < RHO_FLAG]
setorder(agree, spearman_rho)
fwrite(agree, file.path(OUT, "ucell_vs_plaid_spearman.csv"))
cat("Spearman UCell vs PLAID: median", round(median(agree$spearman_rho), 3), "range",
    paste(round(range(agree$spearman_rho), 3), collapse = " - "), ";", sum(agree$flag_rho_below_0.5), "programs below", RHO_FLAG, "\n")
print(agree[flag_rho_below_0.5 == TRUE])
