# UM Atlas October 2026: T/NK STATE groups across the WHOLE cohort (Paul, 2026-10-09).
# 04d refitted the 7 tumour-normal pairs only. Here every sample with >= 30 T/NK cells is used
# (non-malignant, primary, metastasis), ~ origin + (1 | Patient_nr), hmc. The patient intercept carries the
# pairing where there is one (7 patients with a non-malignant + a primary) and the repeat libraries of one
# tumour; unpaired samples inform the origin effect. Contrasts: primary vs non-malignant, metastasis vs
# non-malignant, metastasis vs primary (no patient has both a primary and a metastasis, so the last is unpaired).
# State grouping and seed as in 04d; counts from the TME metadata dump, so no Seurat object is loaded.
# Output: Outputs/sccomp_paired_TNK_2026b/cohort_state_*.csv
suppressPackageStartupMessages({ library(sccomp); library(data.table); library(dplyr) })
stopifnot("CMDSTAN is unset - was R started with --vanilla?" = nzchar(Sys.getenv("CMDSTAN")))
base <- ".."
OUT <- file.path(base, "Outputs/sccomp_paired_TNK_2026b")
draws <- file.path(OUT, "draws_cohort"); dir.create(draws, showWarnings = FALSE)
seed <- 20261005L
min_tnk <- 30

state <- c(CD4_N = "Naive", CD8_N = "Naive",
           CD4_EM = "Memory/effector", CD8_EM = "Memory/effector", CD8_EMRA = "Memory/effector", CD8_RM = "Memory/effector",
           CD4_EX = "Exhausted", CD8_EX = "Exhausted", CD4_REG = "Treg",
           MAIT = "MAIT + gdT", gdT = "MAIT + gdT", NK_REST = "NK", NK_CYTO = "NK",
           `Proliferating T cells` = "Proliferating T")

m <- fread(file.path(base, "new_data_2026b/tme_meta_full.csv.gz"), select = c("Sample", "Patient_nr", "origin", "SubTyping"))
stopifnot(uniqueN(m$Sample) == 63, nrow(m) == 149062, setequal(unique(m$origin), c("non-malignant", "primary", "metastasis")))
t <- m[SubTyping %in% names(state)][, grp := state[SubTyping]]
n_tnk <- t[, .N, by = .(Sample, Patient_nr, origin)]
keep <- n_tnk[N >= min_tnk, Sample]
fwrite(merge(unique(m[, .(Sample, Patient_nr, origin)]), n_tnk, by = c("Sample", "Patient_nr", "origin"), all.x = TRUE)[
  , `:=`(N = fifelse(is.na(N), 0L, N), kept = Sample %in% keep)][order(origin, Patient_nr)], file.path(OUT, "cohort_state_tnk_per_sample.csv"))
t <- t[Sample %in% keep]
stopifnot(uniqueN(t$Sample) == 51)  # 12 + 22 + 17, from the count check on 2026-10-09

d <- t[, .(count = .N), by = .(Sample, Patient_nr, origin, grp)] |> as.data.frame() |>
  tidyr::complete(tidyr::nesting(Sample, Patient_nr, origin), grp, fill = list(count = 0L)) |> as.data.table()
stopifnot(sum(d$count) == nrow(t), uniqueN(d$grp) == 7)
d[, origin := factor(origin, levels = c("non-malignant", "primary", "metastasis"))]
d[, share := count / sum(count), by = Sample]

set.seed(seed)
r <- sccomp_estimate(as.data.frame(d), formula_composition = ~ origin + (1 | Patient_nr), sample = "Sample", cell_group = "grp",
                     abundance = "count", bimodal_mean_variability_association = TRUE, inference_method = "hmc",
                     output_directory = draws, mcmc_seed = seed, cores = 8) |>
  sccomp_remove_outliers(inference_method = "hmc", output_directory = draws, mcmc_seed = seed, cores = 8) |>
  sccomp_test(contrasts = c("originprimary", "originmetastasis", "originmetastasis - originprimary"))
res <- r |> transmute(contrast = recode(as.character(parameter), originprimary = "primary vs non-malignant",
                                        originmetastasis = "metastasis vs non-malignant",
                                        `originmetastasis - originprimary` = "metastasis vs primary"),
                      cell_group = grp, c_lower, c_effect, c_upper, c_pH0, c_FDR)
stopifnot(!anyNA(res$contrast), nrow(res) == 21)

# cross-check: patient-level mean share per origin, unpaired rank-sum (ignores the pairing)
pm <- d[, .(share = mean(share)), by = .(Patient_nr, origin, grp)]
cmp <- list(c("primary", "non-malignant"), c("metastasis", "non-malignant"), c("metastasis", "primary"))
wx <- rbindlist(lapply(cmp, function(cc) pm[origin %in% cc, {
  a <- share[origin == cc[1]]; b <- share[origin == cc[2]]
  .(contrast = sprintf("%s vs %s", cc[1], cc[2]), n_a = length(a), n_b = length(b), median_a = median(a), median_b = median(b),
    rank_sum_p = suppressWarnings(wilcox.test(a, b, exact = FALSE)$p.value))
}, by = .(cell_group = grp)]))
res <- as.data.table(merge(res, wx, by = c("contrast", "cell_group")))  # sccomp returns a tibble; the final print needs data.table

# how much each group rests on: cells, samples and patients with >= 1 cell
pc <- d[, .(count = sum(count), samples = sum(count > 0)), by = .(origin, cell_group = grp, Patient_nr)]
nn <- pc[, .(cells = sum(count), samples = sum(samples), patients = sum(count > 0),
             top_patient = Patient_nr[which.max(count)], top_patient_share = max(count) / sum(count)), by = .(origin, cell_group)]
fwrite(res, file.path(OUT, "cohort_state_tests.csv"))
fwrite(d, file.path(OUT, "cohort_state_shares.csv"))
fwrite(nn, file.path(OUT, "cohort_state_n.csv"))
print(as.data.frame(res[order(contrast, c_pH0)]), digits = 3); print(as.data.frame(nn), digits = 3)
