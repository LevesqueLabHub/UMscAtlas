# UM Atlas October 2026: discover in the matched pairs, validate in samples from OTHER patients (Paul, 2026-10-09).
# Discovery = the 7 tumour-normal pairs, T/NK state groups (04d, tnk_groups_paired.csv, grouping == "state").
# Validation = every kept sample from patients NOT in those 7 pairs: 5 non-malignant (4 patients), 15 primary
# (11 patients), 17 metastasis (13 patients). Same groups, same model (~ origin + (1 | Patient_nr), hmc), so the
# validation fit shares no patient with the discovery fit. 04e used all samples, so its pairs were in both.
# Declared before the fit: a discovery call "holds" when the validation effect has the same sign and sccomp
# c_pH0 < 0.05; "consistent" = same sign only.
# Counts: 04e's cohort_state_shares.csv (>= 30 T/NK cells per sample). Output: Outputs/sccomp_paired_TNK_2026b/validation_*.csv
suppressPackageStartupMessages({ library(sccomp); library(data.table); library(dplyr) })
stopifnot("CMDSTAN is unset - was R started with --vanilla?" = nzchar(Sys.getenv("CMDSTAN")))
base <- ".."
OUT <- file.path(base, "Outputs/sccomp_paired_TNK_2026b")
draws <- file.path(OUT, "draws_validation"); dir.create(draws, showWarnings = FALSE)
seed <- 20261005L

d <- fread(file.path(OUT, "cohort_state_shares.csv"))
disc <- fread(file.path(OUT, "tnk_groups_paired.csv"))[grouping == "state"]
pairs <- intersect(d[origin == "non-malignant", Patient_nr], d[origin == "primary", Patient_nr])
stopifnot(length(pairs) == 7, nrow(disc) == 7)
v <- d[!Patient_nr %in% pairs]
stopifnot(uniqueN(v$Sample) == 37, uniqueN(v[origin == "non-malignant", Sample]) == 5, !any(v$Patient_nr %in% pairs))
v[, origin := factor(origin, levels = c("non-malignant", "primary", "metastasis"))]

set.seed(seed)
r <- sccomp_estimate(as.data.frame(v[, .(Sample, Patient_nr, origin, grp, count)]), formula_composition = ~ origin + (1 | Patient_nr),
                     sample = "Sample", cell_group = "grp", abundance = "count", bimodal_mean_variability_association = TRUE,
                     inference_method = "hmc", output_directory = draws, mcmc_seed = seed, cores = 8) |>
  sccomp_remove_outliers(inference_method = "hmc", output_directory = draws, mcmc_seed = seed, cores = 8) |>
  sccomp_test(contrasts = c("originprimary", "originmetastasis", "originmetastasis - originprimary"))
val <- as.data.table(r |> transmute(contrast = recode(as.character(parameter), originprimary = "primary vs non-malignant",
                                                      originmetastasis = "metastasis vs non-malignant",
                                                      `originmetastasis - originprimary` = "metastasis vs primary"),
                                    cell_group = grp, v_lower = c_lower, v_effect = c_effect, v_upper = c_upper, v_pH0 = c_pH0, v_FDR = c_FDR))
stopifnot(!anyNA(val$contrast), nrow(val) == 21)
fwrite(val, file.path(OUT, "validation_tests.csv"))

# discovery vs validation, primary vs non-malignant
cmp <- merge(disc[, .(cell_group, d_lower = c_lower, d_effect = c_effect, d_upper = c_upper, d_pH0 = c_pH0)],
             val[contrast == "primary vs non-malignant"], by = "cell_group")
cmp[, `:=`(discovery_call = d_pH0 < 0.05, same_sign = sign(d_effect) == sign(v_effect),
           holds = d_pH0 < 0.05 & sign(d_effect) == sign(v_effect) & v_pH0 < 0.05,
           disc_in_val_ci = d_effect >= v_lower & d_effect <= v_upper)]
fwrite(cmp, file.path(OUT, "validation_vs_discovery.csv"))
fwrite(v[, .(cells = sum(count)), by = .(Sample, Patient_nr, origin)][order(origin, Patient_nr)], file.path(OUT, "validation_samples.csv"))
print(as.data.frame(cmp[order(d_pH0)]), digits = 3); print(as.data.frame(val[contrast != "primary vs non-malignant"][order(contrast, v_pH0)]), digits = 3)
