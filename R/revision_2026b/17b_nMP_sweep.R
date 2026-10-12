#!/usr/bin/env Rscript
# 17b: Model-selection sweep for the number of GeneNMF meta-programs (nMP).
# Re-cuts the CACHED program pool from 17 (no multiNMF recompute) over nMP=3:10,
# so we can actually justify the MP count instead of asserting nMP=5.
# getMetaPrograms uses deterministic hclust -> no seed needed.

suppressMessages({
  library(qs2); library(GeneNMF); library(dplyr); library(tidyr); library(ggplot2)
})

# Optional arg: the GeneNMF output folder (default = the January run's relative path).
# October 2026 object: Rscript 17b_nMP_sweep.R ../Outputs/GeneNMF_2026b
args     <- commandArgs(trailingOnly = TRUE)
gnmf_dir <- if (length(args)) args[1] else "Outputs/GeneNMF"
cache    <- file.path(gnmf_dir, "geneNMF_metaprograms.qs2")
out_dir  <- file.path(gnmf_dir, "nMP_sweep")
stopifnot("GeneNMF cache not found" = file.exists(cache))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

programs <- qs2::qs_read(cache, nthreads = 8)$programs   # 192 per-sample NMF programs
nMP_grid <- 3:10
CURRENT  <- if (length(args) >= 2) as.integer(args[2]) else 5L   # the nMP report 17 uses (2nd arg)

# Faithful to 17: only nMP varies; every other getMetaPrograms arg = package default
per_mp <- lapply(nMP_grid, function(n) {
  mp <- GeneNMF::getMetaPrograms(programs, nMP = n)
  m  <- mp$metaprograms.metrics
  m$requested_nMP <- n
  m$MP <- rownames(m)
  m
})
per_mp_df <- dplyr::bind_rows(per_mp)

# One-row-per-nMP summary. min silhouette / min coverage are the honest guards:
# a single incoherent or one-patient MP is what over-splitting produces first.
summary_df <- per_mp_df |>
  dplyr::group_by(requested_nMP) |>
  dplyr::summarise(
    realized_nMP     = dplyr::n(),          # remove.empty can drop some
    mean_silhouette  = mean(silhouette),
    min_silhouette   = min(silhouette),
    n_neg_silhouette = sum(silhouette < 0),
    mean_coverage    = mean(sampleCoverage),
    min_coverage     = min(sampleCoverage),
    n_low_coverage   = sum(sampleCoverage < 0.5),
    mean_similarity  = mean(meanSimilarity),
    .groups = "drop"
  )

readr::write_csv(per_mp_df,  file.path(out_dir, "nMP_sweep_per_MP.csv"))
readr::write_csv(summary_df, file.path(out_dir, "nMP_sweep_summary.csv"))
cat("=== nMP sweep summary ===\n"); print(as.data.frame(summary_df), digits = 3)

# Elbow plot: silhouette (coherence) + coverage (recurrence) vs nMP, mean & worst-case.
plot_df <- summary_df |>
  dplyr::select(requested_nMP, mean_silhouette, min_silhouette,
                mean_coverage, min_coverage) |>
  tidyr::pivot_longer(-requested_nMP, names_to = "metric", values_to = "value") |>
  tidyr::separate(metric, into = c("stat", "quantity"), sep = "_", extra = "merge")

p <- ggplot(plot_df, aes(requested_nMP, value, colour = stat)) +
  geom_vline(xintercept = CURRENT, linetype = "dashed", colour = "grey50") +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey70") +
  geom_line(linewidth = 0.9) + geom_point(size = 2.4) +
  facet_wrap(~ quantity, scales = "free_y") +
  scale_x_continuous(breaks = nMP_grid) +
  scale_colour_manual(values = c(mean = "#2c3e50", min = "#c0392b"),
                      labels = c(mean = "mean across MPs", min = "worst MP")) +
  labs(title = "GeneNMF nMP model selection (cached program pool re-cut)",
       subtitle = sprintf("Dashed line = nMP used in the report (%d). Pick the largest nMP where the worst MP stays coherent (silhouette > 0) and recurrent (coverage high).", CURRENT),
       x = "number of meta-programs (nMP)", y = NULL, colour = NULL) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "top", plot.subtitle = element_text(size = 9))

ggsave(file.path(out_dir, "nMP_sweep_elbow.png"), p,
       width = 11, height = 5, dpi = 300, bg = "white")
cat("\nWrote:", file.path(out_dir, "nMP_sweep_elbow.png"), "\n")
