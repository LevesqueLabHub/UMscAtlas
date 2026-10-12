#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# 25 - ECODA (scECODA CLR composition) vs GloScope on the UM atlas
#
# Question: on this cohort, does the compositional representation recover the
# known sample groupings as well as GloScope's KL divergence between samples?
# Halter, Andreatta & Carmona (bioRxiv 10.64898/2026.03.27.714811) report that
# CLR-transformed cell-type proportions beat GloScope across 11 cohorts. This
# script tests that claim on the UM atlas with the paper's own metric set.
#
# Reuses the GloScope KL divergence matrices already cached by
# 03b_GloScope_Full_Cohort.Rmd - GloScope is NOT re-run.
#
# Scope note: only groupings whose cached GloScope matrix was computed on ALL
# cells of each sample are included. The malignant-only / melanocytic-only
# matrices are cell subsets, where "cell-type composition" is close to
# degenerate, so comparing the two representations there is not meaningful.
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(scECODA)
  library(data.table)
  library(dplyr)
  library(ggplot2)
})

set.seed(123)

work_dir   <- ".."
gs_dir     <- file.path(work_dir, "GloScope_Full_Outputs_2026")
out_dir    <- file.path(work_dir, "ECODA_Benchmark_2026")
meta_gz    <- file.path(work_dir, "new_data_2026/new_atlas_meta_full.csv.gz")
dir.create(out_dir, showWarnings = FALSE)

# --- 1. inputs -------------------------------------------------------------

sample_meta <- fread(file.path(gs_dir, "sample_metadata.csv"), data.table = FALSE)

anno_cols <- c(
  ECODA_LR            = "majority_celltype",   # 13 broad types (report default)
  ECODA_HR            = "SubTyping",           # 41 types incl. T/macrophage subsets
  ECODA_unsup_res0.8  = "RNA_snn_res.0.8",     # 51 Leiden clusters, no annotation
  ECODA_unsup_res0.4  = "RNA_snn_res.0.4",     # 36 clusters - resolution control
  ECODA_unsup_res1    = "seurat_clusters"      # 54 clusters - resolution control
)

cell_meta <- fread(
  cmd = paste("zcat", shQuote(meta_gz)),
  select = c("Sample", unname(anno_cols)),
  data.table = FALSE
)
cat("cells:", nrow(cell_meta), " samples:", length(unique(cell_meta$Sample)), "\n")

# --- 2. representations ----------------------------------------------------

# CLR composition matrix (samples x cell types) via the scECODA pipeline.
clr_matrix <- function(df, celltype_col) {
  df <- df[!is.na(df[[celltype_col]]) & df[[celltype_col]] != "NA", ]
  counts <- get_celltype_counts(df, sample_col = "Sample",
                                celltype_col = celltype_col)
  clr <- as.matrix(calc_clr(replace_zeros(counts)))
  # get_celltype_counts returns celltypes x samples; orient samples as rows
  if (!all(rownames(clr) %in% sample_meta$Sample)) clr <- t(clr)
  # positive control: a CLR vector sums to zero across its components
  stopifnot(
    "CLR rows are not samples" = all(rownames(clr) %in% sample_meta$Sample),
    "CLR not centred per sample" = max(abs(rowSums(clr))) < 1e-8
  )
  clr
}

# raw proportions, Euclidean - isolates what the CLR transform itself buys
prop_matrix <- function(df, celltype_col) {
  df <- df[!is.na(df[[celltype_col]]) & df[[celltype_col]] != "NA", ]
  counts <- get_celltype_counts(df, sample_col = "Sample",
                                celltype_col = celltype_col)
  m <- as.matrix(counts)
  if (!all(rownames(m) %in% sample_meta$Sample)) m <- t(m)
  m / rowSums(m)
}

clr_lr    <- clr_matrix(cell_meta, anno_cols[["ECODA_LR"]])
clr_hr    <- clr_matrix(cell_meta, anno_cols[["ECODA_HR"]])
clr_u08   <- clr_matrix(cell_meta, anno_cols[["ECODA_unsup_res0.8"]])
clr_u04   <- clr_matrix(cell_meta, anno_cols[["ECODA_unsup_res0.4"]])
clr_u10   <- clr_matrix(cell_meta, anno_cols[["ECODA_unsup_res1"]])
prop_hr   <- prop_matrix(cell_meta, anno_cols[["ECODA_HR"]])

# HVC selection: keep the cell types carrying 50% of the compositional variance
hr_var  <- sort(apply(clr_hr, 2, var), decreasing = TRUE)
hr_keep <- names(hr_var)[cumsum(hr_var) / sum(hr_var) <= 0.5]
if (length(hr_keep) < 2) hr_keep <- names(hr_var)[1:2]
clr_hvc <- clr_hr[, hr_keep, drop = FALSE]
cat("HVC kept", length(hr_keep), "of", ncol(clr_hr), "cell types:",
    paste(hr_keep, collapse = ", "), "\n")

feature_spaces <- list(
  ECODA_HR            = clr_hr,
  ECODA_HR_HVC        = clr_hvc,
  ECODA_LR            = clr_lr,
  ECODA_unsup_res0.4  = clr_u04,
  ECODA_unsup_res0.8  = clr_u08,
  ECODA_unsup_res1    = clr_u10,
  Proportions_noCLR   = prop_hr
)

# --- 3. groupings (whole-sample GloScope matrices only) --------------------

jobs <- tibble::tribble(
  ~matrix,                                      ~group_col,             ~label,
  "divergence_origin_full_cohort",              "origin",               "Disease origin",
  "divergence_healthy_vs_primary",              "origin",               "Non-malignant vs primary",
  "divergence_organ_full_cohort",               "organ",                "Organ",
  "divergence_gender_full_cohort",              "Gender",               "Gender",
  "divergence_age_full_cohort",                 "Age_group",            "Age group",
  "divergence_primary_mutation_full",           "primary_mutation",     "Primary mutation",
  "divergence_secondary_mutation_full",         "secondary_mutation",   "Secondary mutation (all)",
  "divergence_secondary_mutation_pUM_full",     "secondary_mutation",   "Secondary mutation (pUM)",
  "divergence_secondary_mutation_choroid_only", "secondary_mutation",   "Secondary mutation (choroid)",
  "divergence_location_pUM_full",               "LocationPrimaryTumor", "Tumour location (pUM)",
  "divergence_healthy_uveal_sites",             "healthy_location",     "Non-malignant uveal site"
)

# --- 4. scoring ------------------------------------------------------------
# Paper's metric set: ANOSIM R, group-adjusted modularity (k = sqrt(n)),
# mean silhouette width, ARI of hierarchical + PAM clustering.
score_all <- function(dist_obj, labels) {
  n <- length(labels)
  a <- vegan::anosim(dist_obj, labels, permutations = 999)
  data.frame(
    anosim     = a$statistic,
    anosim_p   = a$signif,
    modularity = calc_modularity(dist_obj, labels,
                                 knn_k = max(3, round(sqrt(n))), digits = 6),
    silhouette = calc_sil(dist_obj, labels, digits = 6),
    ari        = calc_ari(dist_obj, labels, digits = 6)
  )
}

rows <- list()
for (i in seq_len(nrow(jobs))) {
  mf <- file.path(gs_dir, paste0(jobs$matrix[i], ".qs2"))
  if (!file.exists(mf)) { warning("missing: ", mf); next }
  dm <- qs2::qs_read(mf, nthreads = 1)
  s  <- rownames(dm)
  g  <- sample_meta[[jobs$group_col[i]]][match(s, sample_meta$Sample)]
  keep <- !is.na(g) & s %in% rownames(clr_hr)
  s <- s[keep]; g <- as.character(g[keep])
  if (length(unique(g)) < 2 || length(s) < 4) next

  reps <- c(
    list(GloScope = as.dist(dm[s, s])),
    lapply(feature_spaces, function(m) dist(m[s, , drop = FALSE]))
  )
  # shuffled labels on the best-resolution ECODA space = null floor
  shuffled <- sample(g)

  for (rep_name in names(reps)) {
    rows[[length(rows) + 1]] <- cbind(
      label = jobs$label[i], n_samples = length(s),
      n_groups = length(unique(g)), method = rep_name,
      score_all(reps[[rep_name]], g)
    )
  }
  rows[[length(rows) + 1]] <- cbind(
    label = jobs$label[i], n_samples = length(s),
    n_groups = length(unique(g)), method = "Shuffled_baseline",
    score_all(reps$ECODA_HR, shuffled)
  )
  cat("done:", jobs$label[i], "(n =", length(s), ")\n")
}

res <- bind_rows(rows)
fwrite(res, file.path(out_dir, "ecoda_vs_gloscope_scores.csv"))

# --- 5. aggregate exactly as the paper does: mean rank per method ----------
metrics <- c("anosim", "modularity", "silhouette", "ari")

ranked <- res |>
  tidyr::pivot_longer(all_of(metrics), names_to = "metric", values_to = "value") |>
  group_by(label, metric) |>
  mutate(rank = rank(-value, ties.method = "average")) |>
  ungroup()

summary_tbl <- ranked |>
  group_by(method) |>
  summarise(avg_rank = mean(rank), .groups = "drop") |>
  left_join(
    ranked |> group_by(method, metric) |>
      summarise(mean_value = mean(value), .groups = "drop") |>
      tidyr::pivot_wider(names_from = metric, values_from = mean_value),
    by = "method"
  ) |>
  arrange(avg_rank)

fwrite(summary_tbl, file.path(out_dir, "ecoda_vs_gloscope_summary.csv"))
print(as.data.frame(summary_tbl), digits = 3)

per_metric_rank <- ranked |>
  group_by(method, metric) |>
  summarise(avg_rank = mean(rank), .groups = "drop") |>
  tidyr::pivot_wider(names_from = metric, values_from = avg_rank)
fwrite(per_metric_rank, file.path(out_dir, "ecoda_vs_gloscope_rank_by_metric.csv"))

# Groupings where ANY method finds real signal - ranking over groupings that are
# pure noise for every method just adds noise to the average rank.
signal_labels <- res |>
  group_by(label) |>
  summarise(min_p = min(anosim_p), .groups = "drop") |>
  filter(min_p < 0.05) |> pull(label)
cat("\nsignal-bearing groupings (", length(signal_labels), "of",
    length(unique(res$label)), "):", paste(signal_labels, collapse = "; "), "\n")

summary_signal <- ranked |>
  filter(label %in% signal_labels) |>
  group_by(method) |>
  summarise(avg_rank = mean(rank), .groups = "drop") |>
  arrange(avg_rank)
fwrite(summary_signal, file.path(out_dir, "ecoda_vs_gloscope_summary_signal_only.csv"))
cat("\n-- avg rank, signal-bearing groupings only --\n")
print(as.data.frame(summary_signal), digits = 3)

# Head-to-head: best ECODA variant vs GloScope, per grouping
h2h <- res |>
  mutate(family = case_when(method == "GloScope" ~ "GloScope",
                            grepl("^ECODA", method) ~ "ECODA",
                            TRUE ~ "other")) |>
  filter(family != "other") |>
  group_by(label, n_samples, family) |>
  summarise(anosim = max(anosim), ari = max(ari), .groups = "drop") |>
  tidyr::pivot_wider(names_from = family, values_from = c(anosim, ari)) |>
  mutate(winner_anosim = ifelse(anosim_ECODA > anosim_GloScope, "ECODA", "GloScope"),
         winner_ari    = ifelse(ari_ECODA > ari_GloScope, "ECODA", "GloScope")) |>
  arrange(desc(anosim_GloScope))
fwrite(h2h, file.path(out_dir, "ecoda_vs_gloscope_head_to_head.csv"))
cat("\n-- head to head (best ECODA variant vs GloScope) --\n")
print(as.data.frame(h2h), digits = 3)

# --- 6. figure -------------------------------------------------------------
plot_df <- res |>
  filter(method %in% c("GloScope", "ECODA_HR", "ECODA_LR",
                       "ECODA_unsup_res0.8", "Shuffled_baseline")) |>
  tidyr::pivot_longer(all_of(metrics), names_to = "metric", values_to = "value") |>
  mutate(method = factor(method, levels = c("ECODA_unsup_res0.8", "ECODA_HR", "ECODA_LR",
                                            "GloScope", "Shuffled_baseline")),
         metric = factor(metric, levels = metrics,
                         labels = c("ANOSIM R", "Modularity (adj.)",
                                    "Silhouette", "ARI")))

p <- ggplot(plot_df, aes(x = label, y = value, fill = method)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.75) +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  facet_wrap(~metric, ncol = 2, scales = "free_x") +
  coord_flip() +
  scale_fill_manual(values = c(ECODA_unsup_res0.8 = "#1F77B4", ECODA_HR = "#6CACE4",
                               ECODA_LR = "#AEC7E8", GloScope = "#E15759",
                               Shuffled_baseline = "grey70")) +
  labs(title = "Cell-type composition (scECODA) vs GloScope on the UM atlas",
       subtitle = sprintf("%d sample groupings, %d samples max; GloScope matrices reused from report 03b",
                          length(unique(res$label)), max(res$n_samples)),
       x = NULL, y = "Separation score", fill = NULL) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"),
        legend.position = "top")

ggsave(file.path(out_dir, "ecoda_vs_gloscope.png"), p,
       width = 13, height = 9, dpi = 300, bg = "white")

# cell types driving the composition signal
var_df <- tibble::tibble(cell_type = names(hr_var), clr_variance = as.numeric(hr_var))
fwrite(var_df, file.path(out_dir, "clr_variance_by_celltype_HR.csv"))

cat("\nwrote:", out_dir, "\n")
