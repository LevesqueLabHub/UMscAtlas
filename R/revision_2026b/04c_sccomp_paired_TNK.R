# UM Atlas October 2026: a closer look at the paired sccomp analysis (Paul, 2026-10-07).
#  (1) Does the patient term change anything? The 8 tumour + adjacent-normal pairs refitted on
#      majority_celltype with full MCMC (hmc), with and without (1 | Patient_nr), next to the report's
#      pathfinder fit (sccomp's default, an approximation).
#  (2) T/NK subtypes (SubTyping) inside the T/NK compartment, non-malignant vs primary of the same
#      patient, ~ origin + (1 | Patient_nr), hmc. Cross-check: paired Wilcoxon on subtype proportions.
# Counts come from the TME metadata dump (tme_meta_full.csv.gz), so no Seurat object is loaded.
# Output: Outputs/sccomp_paired_TNK_2026b/
suppressPackageStartupMessages({ library(sccomp); library(data.table); library(dplyr); library(ggplot2) })
stopifnot("CMDSTAN is unset - was R started with --vanilla?" = nzchar(Sys.getenv("CMDSTAN")))
base <- ".."
OUT <- file.path(base, "Outputs/sccomp_paired_TNK_2026b"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
draws <- file.path(OUT, "draws"); dir.create(draws, showWarnings = FALSE)
seed <- 20261005L
min_tnk <- 30  # a sample with fewer T/NK cells cannot inform 14 subtype proportions

m <- fread(file.path(base, "new_data_2026b/tme_meta_full.csv.gz"),
           select = c("Sample", "orig.ident", "Patient_nr", "origin", "majority_celltype", "SubTyping"))
stopifnot(uniqueN(m$Sample) == 63, nrow(m) == 149062)
pairs <- m[origin %in% c("non-malignant", "primary"), uniqueN(origin), by = Patient_nr][V1 == 2, Patient_nr]
stopifnot(length(pairs) == 8)
p <- m[Patient_nr %in% pairs & origin %in% c("non-malignant", "primary")]
p[, origin := factor(origin, levels = c("non-malignant", "primary"))]

fit <- function(d, formula, cell_group, method) {
  set.seed(seed)
  sccomp_estimate(d, formula_composition = formula, sample = "Sample", cell_group = cell_group, abundance = "count",
                  bimodal_mean_variability_association = TRUE, inference_method = method,
                  output_directory = draws, mcmc_seed = seed, cores = 8) |>
    sccomp_remove_outliers(inference_method = method, output_directory = draws, mcmc_seed = seed, cores = 8) |>
    sccomp_test()
}
tidy <- function(r, cell_group, model) r |> dplyr::filter(parameter == "originprimary") |>
  transmute(cell_group = .data[[cell_group]], model, c_lower, c_effect, c_upper, c_pH0, c_FDR)
counts <- function(d, col) d[, .(count = .N), by = c("Sample", "Patient_nr", "origin", col)] |>
  as.data.frame() |> tidyr::complete(tidyr::nesting(Sample, Patient_nr, origin), !!rlang::sym(col), fill = list(count = 0L))

# ---- (1) majority_celltype: pathfinder (report) vs hmc, with and without the patient term ----
ct <- counts(p, "majority_celltype")
res1 <- bind_rows(
  tidy(readRDS(file.path(base, "sccomp_Outputs_2026b/sccomp_paired_results.rds")), "majority_celltype", "pathfinder, paired (report)"),
  tidy(fit(ct, ~ origin + (1 | Patient_nr), "majority_celltype", "hmc"), "majority_celltype", "hmc, paired"),
  tidy(fit(ct, ~ origin, "majority_celltype", "hmc"), "majority_celltype", "hmc, unpaired"))
fwrite(res1, file.path(OUT, "celltype_paired_models.csv"))

# ---- (2) T/NK subtypes within the T/NK compartment ----
tnk_types <- c("CD4_N", "CD4_EM", "CD4_REG", "CD4_EX", "CD8_N", "CD8_EM", "CD8_EMRA", "CD8_RM", "CD8_EX",
               "MAIT", "gdT", "NK_REST", "NK_CYTO", "Proliferating T cells")
t <- p[SubTyping %in% tnk_types]
n_tnk <- t[, .N, by = .(Patient_nr, Sample, orig.ident, origin)]
keep_pt <- n_tnk[, all(N >= min_tnk) && .N == 2, by = Patient_nr][V1 == TRUE, Patient_nr]
fwrite(n_tnk[, kept := Patient_nr %in% keep_pt][order(Patient_nr, origin)], file.path(OUT, "tnk_cells_per_sample.csv"))
t <- t[Patient_nr %in% keep_pt]
# subtypes seen in at least 3 of the kept samples
tnk_keep <- t[, uniqueN(Sample), by = SubTyping][V1 >= 3, SubTyping]
tc <- counts(t[SubTyping %in% tnk_keep], "SubTyping")
r_tnk <- fit(tc, ~ origin + (1 | Patient_nr), "SubTyping", "hmc")
saveRDS(r_tnk, file.path(OUT, "sccomp_tnk_paired_hmc.rds"))
res2 <- bind_rows(tidy(r_tnk, "SubTyping", "hmc, paired"),
                  tidy(fit(tc, ~ origin, "SubTyping", "hmc"), "SubTyping", "hmc, unpaired"))

# paired Wilcoxon on the share of each subtype among T/NK cells (exact, floor p = 2 / 2^n_pairs)
prop <- as.data.table(tc)[, share := count / sum(count), by = Sample]
wx <- prop[, {
  w <- dcast(.SD, Patient_nr ~ origin, value.var = "share")
  d <- w$primary - w$`non-malignant`
  .(n_pairs = nrow(w), n_up = sum(d > 0), median_diff = median(d),
    wilcox_p = suppressWarnings(wilcox.test(w$primary, w$`non-malignant`, paired = TRUE, exact = TRUE)$p.value))
}, by = SubTyping]
res2 <- merge(res2, wx, by.x = "cell_group", by.y = "SubTyping", all.x = TRUE)
fwrite(res2, file.path(OUT, "tnk_subtypes_paired.csv"))
fwrite(prop, file.path(OUT, "tnk_subtype_shares.csv"))

# figure: per-patient share in normal vs tumour, lines join the two samples of a patient
lab <- res2 |> dplyr::filter(model == "hmc, paired") |>
  transmute(SubTyping = cell_group, lab = sprintf("%s\nsccomp p = %.3g, Wilcoxon p = %.2g", cell_group, c_pH0, wilcox_p))
pd <- merge(prop, lab, by = "SubTyping")
g <- ggplot(pd, aes(origin, share)) + geom_line(aes(group = Patient_nr), colour = "grey60") +
  geom_point(aes(colour = origin), size = 2) + facet_wrap(~ lab, scales = "free_y", ncol = 4) +
  scale_colour_manual(values = c(`non-malignant` = "#4DAF4A", primary = "#E41A1C"), guide = "none") +
  scale_y_continuous(labels = scales::percent) +
  labs(title = sprintf("T/NK subtypes, tumour vs adjacent normal of the same patient (%d pairs)", length(keep_pt)),
       x = NULL, y = "Share of T/NK cells") + theme_bw(base_size = 10)
ggsave(file.path(OUT, "tnk_paired_shares.png"), g, width = 13, height = 3 * ceiling(length(tnk_keep) / 4), dpi = 300, bg = "white")
ggsave(file.path(OUT, "tnk_paired_shares.pdf"), g, width = 13, height = 3 * ceiling(length(tnk_keep) / 4))
print(as.data.frame(res1)); print(as.data.frame(res2))
