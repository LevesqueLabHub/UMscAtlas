#!/usr/bin/env Rscript
# Fig. 1h GloScope ANOSIM summary heatmap, rebuilt from its table only, for the Levesque lab to restyle
# (Mauro, 2026-10-20). Same plot as manuscript_panels_2026b_gloscope.R; no GloScope or ANOSIM is run here.
#
# Input: anosim_manuscript_panels.csv in the same folder (one row per test; the heatmap uses the rows whose
#   `grouping` is in row_order and whose `cell_set` is in col_order). Columns used: grouping, cell_set, R, p,
#   n_samples, R_1pt, p_1pt (1pt = one sample per patient per group). `groups` lists the samples per group.
# Run:   Rscript Fig1h_heatmap_from_csv.R      -> Fig1h_gloscope_anosim_summary.{png,pdf} in the same folder
# Needs: ggplot2, dplyr, tidyr, readr, scales

suppressPackageStartupMessages({ library(ggplot2); library(dplyr); library(tidyr); library(readr) })

here <- tryCatch(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)))),
                 error = function(e) getwd())
res <- read_csv(file.path(here, "anosim_manuscript_panels.csv"), show_col_types = FALSE)

# Rows (top to bottom) and columns (left to right): edit, reorder or drop entries here.
row_order <- c("Origin", "Non-malignant vs primary", "Secondary mutation (pUM)", "Secondary mutation (choroidal pUM)",
               "Secondary mutation", "Location", "Location (pUM)", "Location (non-malignant)",
               "Primary mutation", "Organ", "Gender", "Age")
col_order <- c("All Cells", "Malignant Cells", "Melanocytic Cells", "Primary Tumors", "Malignant Cells, pUM", "Choroidal pUM")

stars <- function(p) as.character(cut(p, c(-Inf, 0.001, 0.01, 0.05, Inf), c("***", "**", "*", "ns")))
hm <- res |>
  dplyr::filter(grouping %in% row_order, cell_set %in% col_order) |>
  mutate(label = paste0(sprintf("R=%.2f\n%s", R, stars(p)),
                        ifelse(n_samples_1pt < n_samples, sprintf("\n1/pt: R=%.2f %s", R_1pt, stars(p_1pt)), "")))
stopifnot("two results in one cell" = !anyDuplicated(hm[c("grouping", "cell_set")]))
grid_all <- expand_grid(grouping = row_order, cell_set = col_order) |>
  left_join(hm, by = c("grouping", "cell_set")) |>
  mutate(grouping = factor(grouping, levels = rev(row_order)), cell_set = factor(cell_set, levels = col_order))

caption <- paste(strwrap(c(
  "R = ANOSIM R on GloScope KL divergence between samples. *** p <= 0.001 (the floor with 999 permutations), ** p < 0.01, * p < 0.05, ns otherwise; unadjusted.",
  "1/pt = same test with one sample per patient per group (largest library kept), shown where it changes n.",
  "Melanocytic cells = healthy melanocytes in non-malignant samples; for tumours the melanocytic cells are the malignant cells.",
  "Primary Tumors = all cells of the 30 primaries; Malignant Cells, pUM = malignant cells of the primaries only; Choroidal pUM = all cells of the 18 choroidal primaries.",
  "The pUM rows of Malignant Cells are already primaries-only, so they equal the Malignant Cells, pUM column.",
  "Secondary mutation rows: BAP1 / SF3B1 / EIF1AX only, as submitted. Secondary mutation and Location without a qualifier: every sample of that cell set.",
  "Location = sample location (eye sub-site, non-malignant eye sub-site, or metastatic site); groups of one sample (ciliochoroidal, thyroid) are kept.",
  "Age, Melanocytic Cells: non-malignant samples take their patient's age (in the other columns they have none and are left out).",
  "Grey = not applicable to that cell set, or not run."), 112), collapse = "\n")

p <- ggplot(grid_all, aes(cell_set, grouping)) +
  geom_tile(data = dplyr::filter(grid_all, !is.na(R)), aes(fill = R), colour = "white", linewidth = 0.6) +
  geom_tile(data = dplyr::filter(grid_all, is.na(R)), fill = "grey92", colour = "white", linewidth = 0.6) +
  geom_text(data = dplyr::filter(grid_all, !is.na(R)), aes(label = label), size = 3.2, lineheight = 0.9) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0.2, limits = c(-0.1, 0.6),
                       oob = scales::squish, name = "ANOSIM R") +
  labs(title = "GloScope ANOSIM Summary",
       subtitle = "R statistic and significance across cell subsets (63 samples, 45 patients)",
       caption = caption, x = NULL, y = NULL) +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_text(face = "bold", angle = 30, hjust = 1),
        axis.text.y = element_text(face = "bold"),
        panel.grid = element_blank(),
        plot.title = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5, size = 10),
        plot.caption = element_text(hjust = 0, size = 8.5, colour = "grey25"),
        plot.caption.position = "plot",
        plot.margin = margin(5.5, 5.5, 5.5, 20))

ggsave(file.path(here, "Fig1h_gloscope_anosim_summary.pdf"), p, width = 13, height = 13)
ggsave(file.path(here, "Fig1h_gloscope_anosim_summary.png"), p, width = 13, height = 13, dpi = 300, bg = "white")
