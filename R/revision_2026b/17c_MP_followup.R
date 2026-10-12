#!/usr/bin/env Rscript
# 17c: scores of the nMP = 3 and nMP = 6 meta-programs for 17_GeneNMF_Cancer_Gene_Modules.Rmd.
# Scorer: PLAID (mean z-score of the MP genes, spatial-scribe's fast scorer; plaid_score.R),
# z-scored over the pooled malignant cells + non-malignant melanocytes so the two are on one scale.
# Agreement with the earlier UCell scores is written to mp_scores_plaid_vs_ucell.csv.
# Scores the malignant cells of the cached scored subset (samples with > 500 malignant cells) and the
# non-malignant melanocytes of the October 2026 atlas, and writes the nMP = 6 cut if it is missing.
# Output: Outputs/GeneNMF_2026b/followup/mp_scores_malignant_and_melanocytes.qs2 (one row per cell).
# All figures and statistics are made in the report; this script only does the scoring (slow, 400 GB).
suppressPackageStartupMessages({ library(Seurat); library(GeneNMF); library(dplyr) })
source("plaid_score.R")
BASE <- "../Outputs/GeneNMF_2026b"
# GN_DIR: another program pool (e.g. the blocklist run of 17e); the malignant cells always come from BASE.
GN  <- Sys.getenv("GN_DIR", BASE)
OUT <- file.path(GN, "followup"); dir.create(OUT, showWarnings = FALSE)
OBJ <- "../new_data_2026b/Objects/full Dataset/final_Subtyping_iCNV_reclustered.qs2"
NC  <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "8"))

programs <- qs2::qs_read(file.path(GN, "geneNMF_metaprograms.qs2"), nthreads = 1)$programs
sets <- list()
for (n in c(3L, 6L)) {
  f <- file.path(GN, sprintf("geneNMF_metaprograms_nMP%d.qs2", n))
  if (!file.exists(f)) qs2::qs_save(GeneNMF::getMetaPrograms(programs, nMP = n), f, nthreads = 1)
  g <- qs2::qs_read(f, nthreads = 1)$metaprograms.genes
  stopifnot(length(g) == n)
  sets <- c(sets, setNames(g, paste0("nMP", n, "_", names(g))))
}
# guard: the nMP = 6 cut must reproduce the 17b sweep
if (file.exists(file.path(GN, "nMP_sweep/nMP_sweep_per_MP.csv"))) {
  sweep <- read.csv(file.path(GN, "nMP_sweep/nMP_sweep_per_MP.csv"))
  m6 <- qs2::qs_read(file.path(GN, "geneNMF_metaprograms_nMP6.qs2"), nthreads = 1)$metaprograms.metrics
  s6 <- sweep[sweep$requested_nMP == 6, ]
  stopifnot(all(abs(as.numeric(m6[s6$MP, "sampleCoverage"]) - s6$sampleCoverage) < 1e-9))
}

keep <- c("orig.ident", "Patient_nr", "origin", "secondary_mutation", "primary_location", "healthy_location")
mal <- qs2::qs_read(file.path(BASE, "malignant_sub_with_mp_scores_nMP3.qs2"), nthreads = 1)
mal[["RNA"]] <- JoinLayers(mal[["RNA"]])
cnt_mal <- LayerData(mal, assay = "RNA", layer = "counts")
mal_meta <- mal@meta.data[, intersect(keep, names(mal@meta.data))]
for (k in setdiff(keep, names(mal_meta))) mal_meta[[k]] <- NA_character_  # healthy_location may be absent from the malignant subset
rm(mal); invisible(gc())
atlas <- qs2::qs_read(OBJ, nthreads = 1)
mel <- subset(atlas, subset = origin == "non-malignant" & majority_celltype == "Melanocytes")
rm(atlas); invisible(gc())
mel[["RNA"]] <- JoinLayers(mel[["RNA"]])
cnt_mel <- LayerData(mel, assay = "RNA", layer = "counts")
stopifnot("gene order differs between the two count matrices" = identical(rownames(cnt_mal), rownames(cnt_mel)))
sc <- plaid_score(cbind(cnt_mal, cnt_mel), sets)
stopifnot(nrow(sc) == ncol(cnt_mal) + ncol(cnt_mel), !anyNA(sc))
scores <- bind_rows(cbind(mal_meta, sc[rownames(mal_meta), names(sets)]),
                    cbind(mel@meta.data[, keep], sc[colnames(mel), names(sets)])) |>
  tibble::rownames_to_column("cell") |> tibble::as_tibble()
qs2::qs_save(scores, file.path(OUT, "mp_scores_malignant_and_melanocytes.qs2"), nthreads = 1)

# Healthy melanocytes scored on their own (z over healthy cells only), for the comparison between
# intra-ocular locations within non-malignant tissue (asked 2026-10-08).
sc_mel <- plaid_score(cnt_mel, sets)
stopifnot(nrow(sc_mel) == ncol(cnt_mel), !anyNA(sc_mel))
qs2::qs_save(cbind(mel@meta.data[, keep], sc_mel[colnames(mel), names(sets)]) |>
               tibble::rownames_to_column("cell") |> tibble::as_tibble(),
             file.path(OUT, "mp_scores_melanocytes_only.qs2"), nthreads = 1)

# Agreement with the earlier UCell scores (same cells, same gene sets): per cell and per patient.
# Only for the original pool, which has UCell scores.
if (file.exists(file.path(OUT, "mp_scores_ucell_malignant_and_melanocytes.qs2"))) {
  uc <- qs2::qs_read(file.path(OUT, "mp_scores_ucell_malignant_and_melanocytes.qs2"), nthreads = 1)
  uc <- uc[match(scores$cell, uc$cell), ]
  stopifnot(all(uc$cell == scores$cell))
  agree <- bind_rows(lapply(names(sets), function(k) {
    pp <- tibble(p = scores$Patient_nr, o = scores$origin, a = scores[[k]], b = uc[[k]]) |>
      group_by(p, o) |> summarise(a = mean(a), b = mean(b), .groups = "drop")
    tibble(MP = k, spearman_cells = cor(scores[[k]], uc[[k]], method = "spearman"),
           spearman_patient_means = cor(pp$a, pp$b, method = "spearman"), n_patient_groups = nrow(pp))
  }))
  readr::write_csv(agree, file.path(OUT, "mp_scores_plaid_vs_ucell.csv"))
  print(agree)
}
cat("done:", nrow(scores), "cells\n")
