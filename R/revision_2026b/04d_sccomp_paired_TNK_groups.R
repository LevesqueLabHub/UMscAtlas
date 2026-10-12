# UM Atlas October 2026: paired sccomp on COMBINED T/NK subtypes (Paul, 2026-10-09).
# Follow-up to 04c: the 14 SubTyping levels are summed into (a) 5 lineage groups and (b) 7 state groups,
# then refitted ~ origin + (1 | Patient_nr), full MCMC, same settings as 04c. The groupings are defined here,
# on biology, before any p-value of the combined groups was seen; they were chosen AFTER the 14-subtype result.
# Counts come from 04c's tnk_subtype_shares.csv (same 7 pairs, same >= 30 T/NK cell rule), so nothing is re-derived.
# Output: Outputs/sccomp_paired_TNK_2026b/tnk_groups_*.csv
suppressPackageStartupMessages({ library(sccomp); library(data.table); library(dplyr) })
stopifnot("CMDSTAN is unset - was R started with --vanilla?" = nzchar(Sys.getenv("CMDSTAN")))
base <- ".."
OUT <- file.path(base, "Outputs/sccomp_paired_TNK_2026b")
draws <- file.path(OUT, "draws_groups"); dir.create(draws, showWarnings = FALSE)
seed <- 20261005L  # same seed as 04c

s <- fread(file.path(OUT, "tnk_subtype_shares.csv"))
stopifnot(uniqueN(s$SubTyping) == 14, uniqueN(s$Patient_nr) == 7, uniqueN(s$Sample) == 14)

groupings <- list(
  lineage = c(CD4_N = "CD4", CD4_EM = "CD4", CD4_REG = "CD4", CD4_EX = "CD4",
              CD8_N = "CD8", CD8_EM = "CD8", CD8_EMRA = "CD8", CD8_RM = "CD8", CD8_EX = "CD8",
              MAIT = "MAIT + gdT", gdT = "MAIT + gdT", NK_REST = "NK", NK_CYTO = "NK",
              `Proliferating T cells` = "Proliferating T"),
  state = c(CD4_N = "Naive", CD8_N = "Naive",
            CD4_EM = "Memory/effector", CD8_EM = "Memory/effector", CD8_EMRA = "Memory/effector", CD8_RM = "Memory/effector",
            CD4_EX = "Exhausted", CD8_EX = "Exhausted", CD4_REG = "Treg",
            MAIT = "MAIT + gdT", gdT = "MAIT + gdT", NK_REST = "NK", NK_CYTO = "NK",
            `Proliferating T cells` = "Proliferating T"))
stopifnot(setequal(unique(s$SubTyping), names(groupings$lineage)), setequal(names(groupings$lineage), names(groupings$state)))

fit <- function(d, method = "hmc") {
  set.seed(seed)
  sccomp_estimate(d, formula_composition = ~ origin + (1 | Patient_nr), sample = "Sample", cell_group = "grp", abundance = "count",
                  bimodal_mean_variability_association = TRUE, inference_method = method,
                  output_directory = draws, mcmc_seed = seed, cores = 8) |>
    sccomp_remove_outliers(inference_method = method, output_directory = draws, mcmc_seed = seed, cores = 8) |>
    sccomp_test()
}

res <- list(); shares <- list()
for (g in names(groupings)) {
  d <- copy(s)[, grp := groupings[[g]][SubTyping]][, .(count = sum(count)), by = .(Sample, Patient_nr, origin, grp)]
  d[, origin := factor(origin, levels = c("non-malignant", "primary"))]
  stopifnot(sum(d$count) == sum(s$count), !anyNA(d$grp))  # combining must not lose or invent cells
  r <- fit(as.data.frame(d)) |> dplyr::filter(parameter == "originprimary") |>
    transmute(grouping = g, cell_group = grp, c_lower, c_effect, c_upper, c_pH0, c_FDR)
  d[, share := count / sum(count), by = Sample]
  wx <- d[, {
    w <- dcast(.SD, Patient_nr ~ origin, value.var = "share"); dd <- w$primary - w$`non-malignant`
    .(n_pairs = nrow(w), n_up = sum(dd > 0), median_diff = median(dd),
      wilcox_p = suppressWarnings(wilcox.test(w$primary, w$`non-malignant`, paired = TRUE, exact = TRUE)$p.value))
  }, by = grp]
  res[[g]] <- merge(r, wx, by.x = "cell_group", by.y = "grp")
  shares[[g]] <- d[, grouping := g]
}
fwrite(rbindlist(res), file.path(OUT, "tnk_groups_paired.csv"))
fwrite(rbindlist(shares), file.path(OUT, "tnk_groups_shares.csv"))
print(as.data.frame(rbindlist(res)[order(grouping, c_pH0)]), digits = 3)
