#!/usr/bin/env Rscript
# 17f: PLAID scores for 17f_GeneNMF_Blocklist_2026b.Rmd (Paul, 2026-10-08).
# Scores, on the malignant cells of the cached scored subset (samples with > 500 malignant cells, the same
# cells as report 17): the blocklist MPs of both 17e runs (ribosome + mitochondria removed; plus immediate-early /
# heat-shock removed), the original MPs, and reference gene sets from the uveal melanoma literature.
# All sets are scored in ONE plaid_score call so they share one z-scale.
# Outputs (Outputs/GeneNMF_2026b/blocklist_report/): per-cell scores, per-sample means with the iCNV chr3-loss /
# chr8-gain fractions, MP gene annotations (chromosome, span). Statistics and figures are made in the report.
suppressPackageStartupMessages({ library(Seurat); library(dplyr) })
source("plaid_score.R")
GN  <- "../Outputs/GeneNMF_2026b"
OUT <- file.path(GN, "blocklist_report"); dir.create(OUT, showWarnings = FALSE)

runs <- c(ribo_mito = "blocklist_test_ribo_mito", ribo_mito_stress = "blocklist_test")
mp_sets <- list(); mp_tab <- list()
for (r in names(runs)) for (n in c(3L, 6L)) {
  mp <- qs2::qs_read(file.path(GN, runs[[r]], sprintf("geneNMF_metaprograms_blocklist_nMP%d.qs2", n)), nthreads = 1)
  stopifnot(length(mp$metaprograms.genes) == n)
  ids <- sprintf("%s_nMP%d_%s", r, n, names(mp$metaprograms.genes))
  mp_sets <- c(mp_sets, setNames(mp$metaprograms.genes, ids))
  m <- mp$metaprograms.metrics; m[] <- lapply(m, function(x) as.numeric(as.character(x)))
  mp_tab[[paste(r, n)]] <- tibble(id = ids, run = r, nMP = n, MP = names(mp$metaprograms.genes),
                                  sampleCoverage = m[names(mp$metaprograms.genes), "sampleCoverage"],
                                  silhouette = m[names(mp$metaprograms.genes), "silhouette"],
                                  n_genes = lengths(mp$metaprograms.genes))
}
for (n in c(3L, 6L)) {
  g <- qs2::qs_read(file.path(GN, sprintf("geneNMF_metaprograms_nMP%d.qs2", n)), nthreads = 1)$metaprograms.genes
  mp_sets <- c(mp_sets, setNames(g, sprintf("original_nMP%d_%s", n, names(g))))
}
# Reference sets, each with its source (cited in the report). Single genes are scored as one-gene sets.
ref_sets <- list(
  ref_HLA_class_I  = c("HLA-A", "HLA-B", "HLA-C", "HLA-E", "B2M"),          # inflammatory phenotype, Maat 2008
  ref_HLA_class_II = c("HLA-DRA", "HLA-DRB1", "HLA-DPA1", "HLA-DPB1", "HLA-DQA1", "CD74"),  # Maat 2008
  ref_PRAME        = "PRAME",                                               # Field 2016 (SF3B1, class 1 PRAME+)
  ref_MIF          = "MIF",                                                 # Repp 2000; Ambrosini 2022
  ref_melanocytic  = c("MITF", "PMEL", "TYR", "DCT", "MLANA", "TYRP1"))     # melanocytic differentiation
sets <- c(mp_sets, ref_sets)

mal <- qs2::qs_read(file.path(GN, "malignant_sub_with_mp_scores_nMP3.qs2"), nthreads = 1)
mal[["RNA"]] <- JoinLayers(mal[["RNA"]])
cnt <- LayerData(mal, assay = "RNA", layer = "counts")
miss <- setdiff(unlist(ref_sets), rownames(cnt))
if (length(miss)) cat("reference genes absent from the object:", miss, "\n")
keep <- c("orig.ident", "Patient_nr", "origin", "secondary_mutation", "primary_location", "has_loss_chr3", "has_dupli_chr8")
meta <- mal@meta.data[, keep]
rm(mal); invisible(gc())
sc <- plaid_score(cnt, sets)
stopifnot(nrow(sc) == ncol(cnt), !anyNA(sc), all(names(sets) %in% colnames(sc)))
cells <- cbind(meta, sc[rownames(meta), names(sets)]) |> tibble::rownames_to_column("cell") |> as_tibble()
qs2::qs_save(cells, file.path(OUT, "scores_per_cell.qs2"), nthreads = 1)

# per-sample means (the report averages samples to patients)
samples <- cells |> group_by(orig.ident, Patient_nr, origin, secondary_mutation, primary_location) |>
  summarise(n_cells = n(), chr3_loss = mean(has_loss_chr3, na.rm = TRUE), chr8_gain = mean(has_dupli_chr8, na.rm = TRUE),
            across(all_of(names(sets)), mean), .groups = "drop")
stopifnot(!anyDuplicated(samples$orig.ident))
readr::write_csv(samples, file.path(OUT, "scores_per_sample.csv"))

# MP gene annotation: rank, chromosome, genomic span
gw <- GenomicFeatures::genes(EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86)
gw <- gw[order(-BiocGenerics::width(gw))]; gw <- gw[!duplicated(gw$gene_name)]
ann <- bind_rows(lapply(names(mp_sets), function(k) tibble(id = k, rank = seq_along(mp_sets[[k]]), gene = mp_sets[[k]]))) |>
  mutate(chr = as.character(GenomicRanges::seqnames(gw))[match(gene, gw$gene_name)],
         span_kb = BiocGenerics::width(gw)[match(gene, gw$gene_name)] / 1000)
readr::write_csv(ann, file.path(OUT, "mp_genes_annotated.csv"))
readr::write_csv(bind_rows(mp_tab), file.path(OUT, "mp_metrics.csv"))
cat("done:", nrow(cells), "cells,", nrow(samples), "samples,", length(sets), "sets\n")
