#!/usr/bin/env Rscript
# 17d: does ANY number of GeneNMF meta-programs give a mutation-associated program?
# 2026-10-05, follow-up to the meeting question "is 6 MPs enough to find a BAP1 program".
# Re-cuts the cached program pool at nMP = 3..10 (as 17b), scores every MP with PLAID (mean z-score of its genes,
# plaid_score.R; was UCell) on the malignant
# cells (same cells as 17), averages per patient, and tests secondary mutation (BAP1/SF3B1/EIF1AX) by
# Kruskal-Wallis in choroidal primaries (main question, no location confound) and in all primaries.
# BH correction across all MPs of all cuts within each cell set.
suppressPackageStartupMessages({
  library(GeneNMF); library(Seurat); library(dplyr); library(tidyr); library(ggplot2)
})
BASE <- "../Outputs/GeneNMF_2026b"
# GN_DIR: another program pool (e.g. the blocklist run of 17e); the malignant cells always come from BASE.
GN  <- Sys.getenv("GN_DIR", BASE)
OUT <- file.path(GN, "nMP_mutation_sweep"); dir.create(OUT, showWarnings = FALSE)
NC  <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "8"))

programs <- qs2::qs_read(file.path(GN, "geneNMF_metaprograms.qs2"), nthreads = 1)$programs
cuts <- setNames(lapply(3:10, function(n) GeneNMF::getMetaPrograms(programs, nMP = n)), paste0("nMP", 3:10))
sets <- unlist(lapply(names(cuts), function(k) setNames(cuts[[k]]$metaprograms.genes,
                                                        paste0(k, "_", names(cuts[[k]]$metaprograms.genes)))),
               recursive = FALSE)
metrics <- bind_rows(lapply(names(cuts), function(k) {
  m <- cuts[[k]]$metaprograms.metrics
  data.frame(set = paste0(k, "_", rownames(m)), sampleCoverage = as.numeric(m$sampleCoverage),
             silhouette = as.numeric(m$silhouette))
}))
# guard: must reproduce the 17b sweep (same pool, deterministic hclust)
sw <- read.csv(file.path(GN, "nMP_sweep/nMP_sweep_per_MP.csv"))
sw$set <- paste0("nMP", sw$requested_nMP, "_", sw$MP)
stopifnot(all(abs(metrics$sampleCoverage - sw$sampleCoverage[match(metrics$set, sw$set)]) < 1e-9))

mal <- qs2::qs_read(file.path(BASE, "malignant_sub_with_mp_scores_nMP3.qs2"), nthreads = 1)
source("plaid_score.R")
mal[["RNA"]] <- JoinLayers(mal[["RNA"]])
sc <- plaid_score(LayerData(mal, assay = "RNA", layer = "counts"), sets)
stopifnot(all(names(sets) %in% colnames(sc)))
meta <- mal@meta.data[rownames(sc), c("Patient_nr", "origin", "secondary_mutation", "primary_location")]
rm(mal); invisible(gc())

pat <- cbind(meta, sc[, names(sets)]) |>
  dplyr::filter(origin == "primary", secondary_mutation %in% c("BAP1", "SF3B1", "EIF1AX")) |>
  group_by(Patient_nr, secondary_mutation, primary_location) |>
  summarise(across(all_of(names(sets)), mean), .groups = "drop")
write.csv(pat, file.path(OUT, "patient_mean_scores.csv"), row.names = FALSE)

test <- function(d, label) {
  long <- pivot_longer(d, all_of(names(sets)), names_to = "set", values_to = "score")
  long |> group_by(set) |>
    summarise(kruskal_p = kruskal.test(score ~ factor(secondary_mutation))$p.value,
              n_patients = n(),
              groups = paste(sprintf("%s %d", names(table(secondary_mutation)), table(secondary_mutation)), collapse = ", "),
              .groups = "drop") |>
    mutate(cells = label, padj_BH = p.adjust(kruskal_p, "BH"))
}
res <- bind_rows(test(dplyr::filter(pat, primary_location == "choroid"), "choroidal primaries"),
                 test(pat, "all primaries")) |>
  left_join(metrics, by = "set") |>
  mutate(nMP = as.integer(sub("nMP(\\d+)_.*", "\\1", set)))
write.csv(res, file.path(OUT, "mutation_tests_all_cuts.csv"), row.names = FALSE)

p <- ggplot(res, aes(factor(nMP), -log10(kruskal_p))) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "grey50") +
  geom_jitter(aes(colour = silhouette < 0, size = sampleCoverage), width = 0.15, height = 0, alpha = 0.8) +
  facet_wrap(~ cells) +
  scale_colour_manual(values = c(`FALSE` = "black", `TRUE` = "firebrick"),
                      labels = c(`FALSE` = "separated (silhouette > 0)", `TRUE` = "not separated (silhouette < 0)"), name = NULL) +
  scale_size_continuous(range = c(1, 4), name = "share of samples\nwith the program") +
  labs(title = "Secondary mutation (BAP1 / SF3B1 / EIF1AX) vs every meta-program, nMP = 3 to 10",
       subtitle = sprintf("One point per MP; per-patient mean MP score (PLAID), Kruskal-Wallis. Smallest BH-adjusted p: choroid %.2g, all primaries %.2g",
                          min(res$padj_BH[res$cells == "choroidal primaries"]), min(res$padj_BH[res$cells == "all primaries"])),
       x = "number of meta-programs (nMP)", y = "-log10 p (unadjusted)") +
  theme_minimal(base_size = 12) + theme(legend.position = "bottom")
for (ext in c("png", "pdf")) ggsave(file.path(OUT, paste0("mutation_vs_nMP.", ext)), p, width = 12, height = 6, dpi = 300, bg = "white")
print(res |> group_by(cells) |> slice_min(kruskal_p, n = 5) |> as.data.frame())
cat("done\n")
