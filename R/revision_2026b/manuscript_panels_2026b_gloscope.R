#!/usr/bin/env Rscript
# UM atlas - the manuscript's GloScope panels rebuilt on the October 2026 object (2026b), in the
# submitted design: generate_revised_plots_gloscope.R (A2-A7: cmdscale MDS, theme_classic, A6 labelled MDS,
# A7 KL heatmap) and fix_a1_heatmap.R (A1 ANOSIM summary). Reads the 03b divergence caches only; no GloScope
# run. Panel names = figure IDs (Victoria's numbering).
#
# Why a sub-matrix of a cached matrix is a valid divergence matrix: GloScope 1.6.0 with dens = "KNN",
# dist_mat = "KL" fills entry (i, j) with FNN::KL.dist(X_i, X_j, k)[k] (gloscope -> .calc_dist -> .calc_kl),
# i.e. from the cells of samples i and j only, in the fixed atlas PCA. Checked below as well: every all-cells
# cache equals the matching sub-matrix of the 63-sample matrix exactly.
#
# Run: Rscript manuscript_panels_2026b_gloscope.R   (login node, ~1 min, < 2 GB)

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(ggrepel)
  library(vegan)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))))
source(file.path(script_dir, "config_colors_2026b.R"))

V       <- ".."
G       <- file.path(V, "GloScope_Full_Outputs_2026b")            # 03b caches (October object)
M       <- file.path(V, "GloScope_Full_Outputs")                  # March (January object) caches
fig_dir <- file.path(V, "Figures_2026b", "Manuscript_panels")
out_dir <- file.path(V, "Outputs", "GloScope_manuscript_2026b")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

N_PERM <- 999
rd <- function(f) qs2::qs_read(file.path(G, f), nthreads = 1)

# ---- sample metadata -----------------------------------------------------------------------------------
# 03b's table (pooled none/unknown secondary mutation, LocationPrimaryTumor, Age_group, n_cells) plus the raw
# per-sample columns of the October object (samplename, location, unpooled secondary_mutation).
sm  <- read_csv(file.path(G, "sample_metadata.csv"), show_col_types = FALSE)
raw <- read_csv(file.path(V, "new_data_2026b", "full_sample_metadata.csv"), show_col_types = FALSE,
                col_select = c(Sample, samplename, location, secondary_mutation, origin, Patient_nr))
stopifnot("03b and object sample tables differ" = setequal(sm$Sample, raw$Sample), nrow(sm) == 63,
          all(sm$origin == raw$origin[match(sm$Sample, raw$Sample)]),
          all(sm$Patient_nr == raw$Patient_nr[match(sm$Sample, raw$Sample)]))
meta <- sm |>
  left_join(raw |> select(Sample, samplename, location, smut_raw = secondary_mutation), by = "Sample") |>
  mutate(
    origin = factor(origin, levels = c("non-malignant", "primary", "metastasis")),
    eye_location = case_when(origin == "metastasis" ~ "metastasis",
                             origin == "non-malignant" ~ paste("non-malignant", healthy_location),
                             TRUE ~ primary_location),          # primary_location keeps pUM06 ciliochoroidal
    age_bin = case_when(is.na(Age_numeric) ~ "unknown", Age_numeric < 40 ~ "<40",
                        Age_numeric < 60 ~ "40-60", TRUE ~ "≥60")
  )
stopifnot(!anyNA(meta$samplename), !anyNA(meta$location), !anyNA(meta$smut_raw), !anyNA(meta$origin),
          "manuscript age bins differ from 03b's Age_group" =
            all(c(Young = "<40", Middle = "40-60", Old = "≥60", Unknown = "unknown")[meta$Age_group] == meta$age_bin),
          "samplename not unique" = !anyDuplicated(meta$samplename))
# Age is a patient attribute: non-malignant samples carry NA and take their patient's age (Victoria, 2026-10-16:
# the sample-level age left every non-malignant sample "unknown"). Used by SFig3b and the Fig1h melanocytic age cell.
age_pt <- meta |> dplyr::filter(!is.na(Age_numeric)) |> distinct(Patient_nr, Age_numeric)
stopifnot("a patient has two ages" = !anyDuplicated(age_pt$Patient_nr))
meta$age_bin_pt <- with(list(x = age_pt$Age_numeric[match(meta$Patient_nr, age_pt$Patient_nr)]),
                        case_when(is.na(x) ~ "unknown", x < 40 ~ "<40", x < 60 ~ "40-60", TRUE ~ "≥60"))
pat <- function(s) meta$Patient_nr[match(s, meta$Sample)]
# Keeps factors: 03b's origin is a factor, and the 1/pt subset is ordered by it, which sets the permutation p.
val <- function(s, col) meta[[col]][match(s, meta$Sample)]

# ---- ANOSIM, exactly as 03b ----------------------------------------------------------------------------
# safe_anosim: vegan ANOSIM, 999 permutations, seed 42 per call. "1/pt": one sample per patient per group,
# the one with most cells (03b's ladder level b).
safe_anosim <- function(d, g) {
  set.seed(42)
  a <- vegan::anosim(as.dist(d), g, permutations = N_PERM)
  c(R = unname(a$statistic), p = a$signif)
}
run_test <- function(key, div, group, row, col) {
  s <- rownames(div)[!is.na(group)]; g <- group[!is.na(group)]
  stopifnot(length(s) >= 3, length(unique(g)) >= 2)
  a <- safe_anosim(div[s, s], g)
  b_keep <- tibble(Sample = s, g = g, n = as.numeric(meta$n_cells[match(s, meta$Sample)]), pt = pat(s)) |>
    group_by(pt, g) |> slice_max(n, n = 1, with_ties = FALSE) |> ungroup()
  b <- safe_anosim(div[b_keep$Sample, b_keep$Sample], b_keep$g)
  tab <- table(as.character(g))
  tibble(key = key, grouping = row, cell_set = col,
         groups = paste(sprintf("%s %d", names(tab), tab), collapse = ", "),
         n_samples = length(s), n_patients = n_distinct(pat(s)), R = a[["R"]], p = a[["p"]],
         n_samples_1pt = nrow(b_keep), R_1pt = b[["R"]], p_1pt = b[["p"]])
}

C_ALL <- "All cells, all samples"; C_PRIM <- "All cells, primaries only"
C_MAL <- "Malignant cells";        C_MEL  <- "Melanocytic cells"
C_CHOR <- "All cells, choroidal primaries"
C_MALP <- "Malignant cells, primaries only"

# 03b's comparisons: cache, grouping column of 03b's sample table, heatmap row and column.
# 03b's smut_pUM (all 30 primaries, none/unknown pooled as a 4th group) is left out on purpose: Fig1h uses
# the submitted definition (BAP1/SF3B1/EIF1AX only) for that cell, smut_pUM_called below, so it is not
# part of the 03b-match check. The other secondary-mutation cells (smut_mal, chor, chor_mal) are already
# called-mutation only in 03b and stay checked.
rows_03b <- tibble::tribble(
  ~key,           ~cache,                                                   ~gcol,                  ~row,                                 ~col,
  "origin_all",   "divergence_origin_full_cohort.qs2",                      "origin",               "Origin",                             C_ALL,
  "origin_mal",   "divergence_origin_malignant_full.qs2",                   "origin",               "Origin",                             C_MAL,
  "origin_mel",   "divergence_origin_melanocytic_full.qs2",                 "origin",               "Origin",                             C_MEL,
  "hs_prim",      "divergence_healthy_vs_primary.qs2",                      "origin",               "Non-malignant vs primary",           C_ALL,
  "smut_mal",     "divergence_secondary_mutation_malignant_full.qs2",       "secondary_mutation",   "Secondary mutation (pUM)",           C_MAL,
  "smut_all",     "divergence_secondary_mutation_full.qs2",                 "secondary_mutation",   "Secondary mutation",                 C_ALL,
  "smut_all_mal", "divergence_secondary_mutation_any_origin_malignant.qs2", "secondary_mutation",   "Secondary mutation",                 C_MAL,
  "chor",         "divergence_secondary_mutation_choroid_only.qs2",         "secondary_mutation",   "Secondary mutation (choroidal pUM)", C_PRIM,
  "chor_mal",     "divergence_secondary_mutation_choroid_malignant.qs2",    "secondary_mutation",   "Secondary mutation (choroidal pUM)", C_MAL,
  "loc_pUM",      "divergence_location_pUM_full.qs2",                       "LocationPrimaryTumor", "Location (pUM)",                     C_PRIM,
  "loc_mal",      "divergence_location_malignant.qs2",                      "LocationPrimaryTumor", "Location (pUM)",                     C_MAL,
  "hsites",       "divergence_healthy_uveal_sites.qs2",                     "healthy_location",     "Location (non-malignant)",           C_ALL,
  "hs_mel",       "divergence_healthy_uveal_sites_melanocytes.qs2",         "healthy_location",     "Location (non-malignant)",           C_MEL,
  "pmut_all",     "divergence_primary_mutation_full.qs2",                   "primary_mutation",     "Primary mutation",                   C_ALL,
  "pmut_pUM",     "divergence_primary_mutation_pUM.qs2",                    "primary_mutation",     "Primary mutation",                   C_PRIM,
  "pmut_mal",     "divergence_primary_mutation_malignant.qs2",              "primary_mutation",     "Primary mutation",                   C_MAL,
  "organ_all",    "divergence_organ_full_cohort.qs2",                       "organ",                "Organ",                              C_ALL,
  "organ_mal",    "divergence_organ_malignant.qs2",                         "organ",                "Organ",                              C_MAL,
  "organ_mel",    "divergence_organ_melanocytic.qs2",                       "organ",                "Organ",                              C_MEL,
  "gender_all",   "divergence_gender_full_cohort.qs2",                      "Gender",               "Gender",                             C_ALL,
  "gender_pUM",   "divergence_gender_pUM.qs2",                              "Gender",               "Gender",                             C_PRIM,
  "gender_mal",   "divergence_gender_malignant.qs2",                        "Gender",               "Gender",                             C_MAL,
  "gender_mel",   "divergence_gender_melanocytic.qs2",                      "Gender",               "Gender",                             C_MEL,
  "age_all",      "divergence_age_full_cohort.qs2",                         "Age_group",            "Age",                                C_ALL,
  "age_pUM",      "divergence_age_pUM.qs2",                                 "Age_group",            "Age",                                C_PRIM,
  "age_mal",      "divergence_age_malignant.qs2",                           "Age_group",            "Age",                                C_MAL
)
res_03b <- bind_rows(lapply(seq_len(nrow(rows_03b)), function(i) {
  r <- rows_03b[i, ]; d <- rd(r$cache)
  run_test(r$key, d, val(rownames(d), r$gcol), r$row, r$col)
}))

# Positive control: the recomputation reproduces 03b's own CSVs exactly (R, p, and 1/pt).
ref <- read_csv(file.path(G, "anosim_summary_full_cohort.csv"), show_col_types = FALSE) |>
  select(key, R_ref = statistic, p_ref = p_value, n_ref = n_samples)
ref_b <- read_csv(file.path(G, "anosim_patient_level_sensitivity_all.csv"), show_col_types = FALSE) |>
  dplyr::filter(startsWith(subset, "(b)")) |> select(key, R_b = R, p_b = p, n_b = n_samples)
chk <- res_03b |> inner_join(ref, by = "key") |> inner_join(ref_b, by = "key")
stopifnot("03b rows not all matched" = nrow(chk) == 26,   # 27 in 03b minus smut_pUM (redefined, see above)
          all(chk$key %in% ref$key), !"smut_pUM" %in% chk$key,
          "R differs from 03b" = all(abs(chk$R - chk$R_ref) < 1e-12),
          "p differs from 03b" = all(chk$p == chk$p_ref),
          "n differs from 03b" = all(chk$n_samples == chk$n_ref),
          "1/pt differs from 03b" = all(abs(chk$R_1pt - chk$R_b) < 1e-12 & chk$p_1pt == chk$p_b &
                                          chk$n_samples_1pt == chk$n_b))

# ---- matrices for the new tests and the panels ------------------------------------------------------------
full <- rd("divergence_origin_full_cohort.qs2")          # all cells, 63 samples
mel  <- rd("divergence_origin_melanocytic_full.qs2")     # healthy melanocytes + malignant cells, >= 100 / sample
mal  <- rd("divergence_origin_malignant_full.qs2")
stopifnot(nrow(full) == 63, setequal(rownames(full), meta$Sample))
# Step-0 check, empirical: all-cells caches computed on other sample sets equal the 63-sample sub-matrix.
for (f in c("divergence_location_pUM_full.qs2", "divergence_secondary_mutation_pUM_full.qs2",
            "divergence_healthy_vs_primary.qs2", "divergence_primary_mutation_full.qs2")) {
  d <- rd(f); stopifnot("cache is not a sub-matrix of the full matrix" = max(abs(d - full[rownames(d), rownames(d)])) == 0)
}
# Caption claim "for tumours the melanocytic cells are the malignant cells": same divergences on shared samples.
s_mm <- intersect(rownames(mal), rownames(mel))
stopifnot("melanocytic and malignant matrices differ on tumours" =
            length(s_mm) == nrow(mal) && max(abs(mal[s_mm, s_mm] - mel[s_mm, s_mm])) == 0)

prim   <- meta$Sample[meta$origin == "primary"]
called <- c("BAP1", "SF3B1", "EIF1AX")
chor_s <- intersect(prim, meta$Sample[meta$primary_location %in% "choroid"])
stopifnot("choroidal pUM count" = length(chor_s) == 18)
sub_g  <- function(d, keep, col) { s <- intersect(rownames(d), keep); list(d = d[s, s], g = val(s, col)) }
new_test <- function(key, d, keep, col, row, cset) { x <- sub_g(d, keep, col); run_test(key, x$d, x$g, row, cset) }

res_new <- bind_rows(
  # requested Fig1h additions
  new_test("loc_all", full, meta$Sample, "location", "Location", C_ALL),
  # Victoria 2026-10-20 (yellow cells): Location in the other cell sets, melanocytic age, choroidal pUM column
  new_test("loc_all_mal", mal, meta$Sample, "location", "Location", C_MAL),
  new_test("loc_all_mel", mel, meta$Sample, "location", "Location", C_MEL),
  new_test("loc_all_prim", full, prim, "location", "Location", C_PRIM),
  new_test("age_mel", mel, meta$Sample[meta$age_bin_pt != "unknown"], "age_bin_pt", "Age", C_MEL),
  new_test("smut_chor", full, intersect(chor_s, meta$Sample[meta$smut_raw %in% called]),
           "smut_raw", "Secondary mutation", C_CHOR),
  new_test("pmut_chor", full, intersect(chor_s, meta$Sample[meta$primary_mutation %in% c("GNAQ", "GNA11", "PLCB4")]),
           "primary_mutation", "Primary mutation", C_CHOR),
  new_test("gender_chor", full, chor_s, "Gender", "Gender", C_CHOR),
  new_test("age_chor", full, intersect(chor_s, meta$Sample[meta$Age_group != "Unknown"]), "Age_group", "Age", C_CHOR),
  new_test("pmut_mel", mel, meta$Sample[meta$primary_mutation %in% c("GNAQ", "GNA11", "PLCB4")],
           "primary_mutation", "Primary mutation", C_MEL),
  new_test("smut_mel", mel, intersect(prim, meta$Sample[meta$smut_raw %in% called]),
           "smut_raw", "Secondary mutation (pUM)", C_MEL),
  new_test("hs_prim_mel", mel, meta$Sample[meta$origin %in% c("non-malignant", "primary")],
           "origin", "Non-malignant vs primary", C_MEL),
  # Mauro 2026-10-28: every variable on malignant cells of primary tumours only (the Malignant Cells column mixes in
  # metastases for primary mutation, gender and age). The pUM rows repeat the Malignant Cells column (checked below).
  new_test("smut_malp", mal, intersect(prim, meta$Sample[meta$smut_raw %in% called]), "smut_raw", "Secondary mutation (pUM)", C_MALP),
  new_test("chor_malp", mal, intersect(chor_s, meta$Sample[meta$smut_raw %in% called]), "smut_raw", "Secondary mutation (choroidal pUM)", C_MALP),
  new_test("loc_malp", mal, intersect(prim, meta$Sample[meta$LocationPrimaryTumor %in% c("choroid", "ciliary body", "iris")]),
           "LocationPrimaryTumor", "Location (pUM)", C_MALP),
  new_test("loc_all_malp", mal, prim, "location", "Location", C_MALP),
  new_test("pmut_malp", mal, intersect(prim, meta$Sample[meta$primary_mutation %in% c("GNAQ", "GNA11", "PLCB4")]),
           "primary_mutation", "Primary mutation", C_MALP),
  new_test("gender_malp", mal, prim, "Gender", "Gender", C_MALP),
  new_test("age_malp", mal, intersect(prim, meta$Sample[meta$Age_group != "Unknown"]), "Age_group", "Age", C_MALP),
  # Fig1h primaries cell and Fig1j subtitle: BAP1/SF3B1/EIF1AX only, the submitted A1/A3 definition
  # (replaces 03b's smut_pUM, which pooled none/unknown as a 4th group)
  new_test("smut_pUM_called", full, intersect(prim, meta$Sample[meta$smut_raw %in% called]),
           "smut_raw", "Secondary mutation (pUM)", C_PRIM)
) |> mutate(new_in_2026b = TRUE)
res <- bind_rows(res_03b |> mutate(new_in_2026b = FALSE), res_new)
R_ <- function(k) res[res$key == k, ]
# Positive control for the pUM-only malignant column: its pUM rows must reproduce the Malignant Cells column.
for (k in c("smut", "chor", "loc")) stopifnot("pUM malignant column does not reproduce the Malignant Cells column" =
  abs(R_(paste0(k, "_malp"))$R - R_(paste0(k, "_mal"))$R) < 1e-12 && R_(paste0(k, "_malp"))$n_samples == R_(paste0(k, "_mal"))$n_samples)
# Positive control for the choroidal subset: same samples and grouping as 03b's choroid-only cache.
stopifnot("choroidal pUM column does not reproduce 03b's chor" =
            abs(R_("smut_chor")$R - R_("chor")$R) < 1e-12 && R_("smut_chor")$n_samples == R_("chor")$n_samples)

# ---- formatting helpers ----------------------------------------------------------------------------------
stars <- function(p) as.character(cut(p, c(-Inf, 0.001, 0.01, 0.05, Inf), c("***", "**", "*", "ns")))
fmt_p <- function(p) sprintf("p = %.3g", p)
anosim_lines <- function(k, what = "") {
  r <- R_(k)
  c(sprintf("ANOSIM%s: R = %.2f, %s", what, r$R, fmt_p(r$p)),
    if (r$n_samples_1pt < r$n_samples)
      sprintf("One sample per patient per group (n = %d): R = %.2f, %s", r$n_samples_1pt, r$R_1pt, fmt_p(r$p_1pt)))
}
count_line <- function(s, col) {
  v <- as.character(val(s, col))
  tab <- table(factor(v, levels = unique(c(names(pal_of[[col]]), v))))
  tab <- tab[tab > 0]
  sprintf("n = %d samples, %d patients (%s)", length(s), n_distinct(pat(s)),
          paste(sprintf("%s %d", names(tab), tab), collapse = ", "))
}
wrap <- function(x, w = 100) paste(unlist(lapply(x, function(l) strwrap(l, w))), collapse = "\n")

origin_shapes <- c("non-malignant" = 16, "primary" = 17, "metastasis" = 18)
pal_of <- list(origin = origin_colors[c("non-malignant", "primary", "metastasis")], smut_raw = second_mutations_col, primary_location = primary_location_colors,
               primary_mutation = primary_mutations_col, organ = organ_colors)

save_gg <- function(p, name, w, h) {
  ggsave(file.path(fig_dir, paste0(name, ".pdf")), p, device = cairo_pdf, width = w, height = h)
  ggsave(file.path(fig_dir, paste0(name, ".png")), p, width = w, height = h, dpi = 300, bg = "white")
}

# Manuscript MDS (build_mds_plot of generate_revised_plots_gloscope.R): cmdscale on the sub-matrix of the
# all-cells matrix for the plotted samples, theme_classic(14), legend right.
mds_plot <- function(s, colour, title, subtitle, legend_title, shape = NULL, label = FALSE) {
  xy <- cmdscale(as.dist(full[s, s]), k = 2)
  pal <- pal_of[[colour]]
  df <- tibble(Sample = s, MDS1 = xy[, 1], MDS2 = xy[, 2]) |> left_join(meta, by = "Sample") |>
    mutate(colour = factor(as.character(.data[[colour]]), levels = intersect(names(pal), as.character(.data[[colour]]))))
  stopifnot("level without a colour" = !anyNA(df$colour), nrow(df) == length(s))
  p <- ggplot(df, aes(MDS1, MDS2, colour = colour)) +
    scale_colour_manual(values = pal, name = legend_title) +
    coord_fixed() +
    labs(title = title, subtitle = subtitle) +
    theme_classic(base_size = 14) +
    theme(legend.position = "right",
          plot.title = element_text(face = "bold", size = 16),
          plot.subtitle = element_text(size = 11.5, lineheight = 1.05),
          axis.title = element_text(size = 14), axis.text = element_text(size = 12),
          legend.title = element_text(size = 12), legend.text = element_text(size = 11))
  if (is.null(shape)) return(p + geom_point(size = 6))
  p <- p + geom_point(aes(shape = .data[[shape]]), size = 6) +
    scale_shape_manual(values = origin_shapes, name = "Origin") +
    guides(colour = guide_legend(override.aes = list(shape = 15, size = 5), order = 2),
           shape = guide_legend(order = 1))
  # point.size = the drawn point size, so labels are pushed off the markers, not just off their centres
  if (label) p <- p + geom_text_repel(aes(label = samplename), size = 3, max.overlaps = Inf, point.size = 6,
                                      box.padding = 0.3, min.segment.length = 0.2, seed = 1, show.legend = FALSE)
  p
}

# ---- Fig1i_mds_origin ----------------------------------------------------------------------------------
s <- meta$Sample
save_gg(mds_plot(s, "origin", "GloScope MDS: Disease Origin (All Cells)",
                 wrap(c(count_line(s, "origin"), anosim_lines("origin_all"))), "Origin"),
        "Fig1i_mds_origin", 10, 8)

# ---- Fig1i_mds_primary_mutation (Victoria labelled this one Fig. 1i too) -----------------------------------
pm <- R_("pmut_all")
save_gg(mds_plot(s, "primary_mutation", "GloScope MDS: Primary Mutations (All Cells)",
                 wrap(c(count_line(s, "primary_mutation"), anosim_lines("pmut_all", sprintf(" GNAQ vs GNA11 vs PLCB4 (n = %d)", pm$n_samples)),
                        "none = the 16 non-malignant samples; none and unknown are plotted, not tested.")),
                 "Primary Mutation"),
        "Fig1i_mds_primary_mutation", 10, 8)

# ---- Fig1j_mds_pUM_secondary_mutation: none and unknown as separate levels (raw object column) ---------------
sc <- R_("smut_pUM_called")
save_gg(mds_plot(prim, "smut_raw", "GloScope MDS: Secondary Mutations in Primary UM (All Cells)",
                 wrap(c(count_line(prim, "smut_raw"),
                        anosim_lines("smut_pUM_called", sprintf(" BAP1 vs SF3B1 vs EIF1AX (n = %d)", sc$n_samples)),
                        "none (no mutation detected) and unknown (not assessed) are plotted, not tested.")),
                 "Secondary Mutation"),
        "Fig1j_mds_pUM_secondary_mutation", 10, 8)

# ---- Fig1k_mds_pUM_location: all 30 primaries incl. ciliochoroidal pUM06 (as the March A2) ------------------
lp <- R_("loc_pUM")
save_gg(mds_plot(prim, "primary_location", "GloScope MDS: Anatomical Location in Primary UM (All Cells)",
                 wrap(c(count_line(prim, "primary_location"),
                        anosim_lines("loc_pUM", sprintf(" choroid vs iris vs ciliary body (n = %d)", lp$n_samples)),
                        "Ciliochoroidal (n = 1) is plotted, not tested.")),
                 "Location"),
        "Fig1k_mds_pUM_location", 10, 8)

# ---- SFig3c_mds_samples_labelled (A6 design: colour organ, shape origin, sample names) -------------------------
oa <- R_("organ_all")
save_gg(mds_plot(s, "organ", "GloScope MDS: Organ/Metastasis Site (All Cells)",
                 wrap(c(count_line(s, "organ"),
                        anosim_lines("organ_all", sprintf(" organ (n = %d, thyroid n = 1 excluded)", oa$n_samples))), 110),
                 "Organ", shape = "origin", label = TRUE),
        "SFig3c_mds_samples_labelled", 12, 9.5)

# ---- Fig1h_gloscope_anosim_summary (A1 design) ---------------------------------------------------------------
row_order <- c("Origin", "Non-malignant vs primary", "Secondary mutation (pUM)", "Secondary mutation (choroidal pUM)",
               "Secondary mutation", "Location", "Location (pUM)",
               "Location (non-malignant)", "Primary mutation", "Organ", "Gender", "Age")
col_lab <- c(setNames("All Cells", C_ALL), setNames("Malignant Cells", C_MAL),
             setNames("Melanocytic Cells", C_MEL), setNames("Primary Tumors", C_PRIM),
             setNames("Malignant Cells, pUM", C_MALP), setNames("Choroidal pUM", C_CHOR))
hm <- res |> dplyr::filter(grouping %in% row_order) |>
  mutate(label = paste0(sprintf("R=%.2f\n%s", R, stars(p)),
                        ifelse(n_samples_1pt < n_samples, sprintf("\n1/pt: R=%.2f %s", R_1pt, stars(p_1pt)), "")))
stopifnot("two results in one Fig1h cell" = !anyDuplicated(hm[c("grouping", "cell_set")]))
grid_all <- tidyr::expand_grid(grouping = row_order, cell_set = names(col_lab)) |>
  left_join(hm, by = c("grouping", "cell_set")) |>
  mutate(grouping = factor(grouping, levels = rev(row_order)),
         cell_set = factor(col_lab[cell_set], levels = col_lab))
cap_h <- wrap(c(
  "R = ANOSIM R on GloScope KL divergence between samples. *** p <= 0.001 (the floor with 999 permutations), ** p < 0.01, * p < 0.05, ns otherwise; unadjusted.",
  "1/pt = same test with one sample per patient per group (largest library kept), shown where it changes n.",
  "Melanocytic cells = healthy melanocytes in non-malignant samples; for tumours the melanocytic cells are the malignant cells.",
  "Primary Tumors = all cells of the 30 primaries; Malignant Cells, pUM = malignant cells of the primaries only; Choroidal pUM = all cells of the 18 choroidal primaries.",
  "The pUM rows of Malignant Cells are already primaries-only, so they equal the Malignant Cells, pUM column.",
  "Secondary mutation rows: BAP1 / SF3B1 / EIF1AX only, as submitted. Secondary mutation and Location without a qualifier: every sample of that cell set.",
  "Location = sample location (eye sub-site, non-malignant eye sub-site, or metastatic site); groups of one sample (ciliochoroidal, thyroid) are kept.",
  "Age, Melanocytic Cells: non-malignant samples take their patient's age (in the other columns they have none and are left out).",
  "Grey = not applicable to that cell set, or not run."), 112)
p_h <- ggplot(grid_all, aes(cell_set, grouping)) +
  geom_tile(data = dplyr::filter(grid_all, !is.na(R)), aes(fill = R), colour = "white", linewidth = 0.6) +
  geom_tile(data = dplyr::filter(grid_all, is.na(R)), fill = "grey92", colour = "white", linewidth = 0.6) +
  geom_text(data = dplyr::filter(grid_all, !is.na(R)), aes(label = label), size = 3.2, lineheight = 0.9) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0.2, limits = c(-0.1, 0.6),
                       oob = scales::squish, name = "ANOSIM R") +
  labs(title = "GloScope ANOSIM Summary",
       subtitle = sprintf("R statistic and significance across cell subsets (%d samples, %d patients)",
                          nrow(meta), n_distinct(meta$Patient_nr)),
       caption = cap_h, x = NULL, y = NULL) +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_text(face = "bold", angle = 30, hjust = 1),
        axis.text.y = element_text(face = "bold"),
        panel.grid = element_blank(),
        plot.title = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5, size = 10),
        plot.caption = element_text(hjust = 0, size = 8.5, colour = "grey25"),
        plot.caption.position = "plot",
        plot.margin = margin(5.5, 5.5, 5.5, 20))
save_gg(p_h, "Fig1h_gloscope_anosim_summary", 13, 13)

# ---- ANOSIM CSV, with the March (January-object) value where the same comparison existed ----------------------
march <- tibble::tribble(
  ~key,              ~march_comparison,                  ~march_cache,
  "origin_all",      "disease_origin_full",              "divergence_origin_full_cohort.rds",
  "origin_mal",      "disease_origin_malignant_full",    "divergence_origin_malignant_full.rds",
  "origin_mel",      "disease_origin_melanocytic_full",  "divergence_origin_melanocytic_full.rds",
  "smut_pUM_called", "secondary_mutation_pUM_full",      "divergence_secondary_mutation_pUM_full.rds",
  "smut_mal",        "secondary_mutation_malignant_full", "divergence_secondary_mutation_malignant_full.rds",
  "smut_all",        "secondary_mutation_full",          "divergence_secondary_mutation_full.rds",
  "loc_all",         "location_full",                    NA,
  "loc_pUM",         "anatomical_location_pUM_full",     "divergence_location_pUM_full.rds",
  "loc_mal",         "location_malignant_full",          "divergence_location_malignant.rds",
  "pmut_all",        "primary_mutation_full",            "divergence_primary_mutation_full.rds",
  "pmut_pUM",        "primary_mutation_pUM_full",        "divergence_primary_mutation_pUM.rds",
  "pmut_mal",        "primary_mutation_malignant_full",  "divergence_primary_mutation_malignant.rds",
  "organ_all",       "organ_full",                       "divergence_organ_full_cohort.rds",
  "organ_mal",       "organ_malignant_full",             "divergence_organ_malignant.rds",
  "organ_mel",       "organ_melanocytic_full",           "divergence_organ_melanocytic.rds",
  "gender_all",      "gender_full",                      "divergence_gender_full_cohort.rds",
  "gender_pUM",      "gender_pUM_full",                  "divergence_gender_pUM.rds",
  "gender_mal",      "gender_malignant",                 "divergence_gender_malignant.rds",
  "gender_mel",      "gender_melanocytic_full",          "divergence_gender_melanocytic.rds",
  "age_all",         "age_full",                         "divergence_age_full_cohort.rds",
  "age_pUM",         "age_pUM_full",                     "divergence_age_pUM.rds",
  "age_mal",         "age_malignant",                    "divergence_age_malignant.rds"
)
march_vals <- read_csv(file.path(V, "Revised_Plots_20260325", "GroupA_GloScope", "A1_anosim_summary_data.csv"),
                       show_col_types = FALSE) |>
  select(march_comparison = comparison, R_march = statistic, p_march = p_value)
march <- march |>
  left_join(march_vals, by = "march_comparison") |>
  mutate(n_samples_march = vapply(march_cache, function(f) if (is.na(f)) NA_integer_ else {
    r <- rownames(readRDS(file.path(M, f))); sum(!grepl("mUM15", r)) }, integer(1)),
    march_note = ifelse(is.na(march_cache), "no cache or code found for this March row; grouping not recorded", NA))
stopifnot("March row missing" = !anyNA(march$R_march))

panel_of <- c(origin_all = "Fig1h; Fig1i_mds_origin; SFig3b_kl_heatmap", pmut_all = "Fig1h; Fig1i_mds_primary_mutation",
              smut_pUM_called = "Fig1h; Fig1j_mds_pUM_secondary_mutation", loc_pUM = "Fig1h; Fig1k_mds_pUM_location",
              organ_all = "Fig1h; SFig3c_mds_samples_labelled")
csv <- res |>
  mutate(cell_set = ifelse(cell_set %in% names(col_lab), col_lab[cell_set], cell_set),
         panel = coalesce(panel_of[key], ifelse(grouping %in% row_order, "Fig1h", NA))) |>
  left_join(march |> select(-march_cache), by = "key") |>
  select(key, grouping, cell_set, groups, n_samples, n_patients, R, p,
         n_samples_1pt, R_1pt, p_1pt, new_in_2026b, panel,
         march_comparison, n_samples_march, R_march, p_march, march_note)
write_csv(csv, file.path(out_dir, "anosim_manuscript_panels.csv"))

# ---- SFig3b_kl_heatmap (A7 design, clustered; Heatmap() built inside each open device) ------------------------
a <- meta[match(rownames(full), meta$Sample), ]
fac <- function(x, lev) factor(x, levels = intersect(lev, x))
cat("SFig3b age: patients without an age =", n_distinct(a$Patient_nr[a$age_bin_pt == "unknown"]), "\n")
age_cols <- c("<40" = unname(ageGroup_col["Young"]), "40-60" = unname(ageGroup_col["Middle"]),
              "≥60" = unname(ageGroup_col["Old"]), unknown = "#E1E1E1FF")
anno <- list(Patient = fac(a$Patient_nr, names(patient_cols)), Origin = fac(as.character(a$origin), names(origin_colors)[c(2, 1, 3)]),
             Organ = fac(a$organ, names(organ_colors)), `Location` = fac(a$location, names(location_colors)),
             Gender = fac(a$Gender, names(gender_col)), `Age group` = fac(a$age_bin_pt, names(age_cols)),
             `Primary mutation` = fac(a$primary_mutation, names(primary_mutations_col)),
             `Secondary mutation` = fac(a$smut_raw, names(second_mutations_col)))
cols <- list(Patient = patient_cols, Origin = origin_colors, Organ = organ_colors, `Location` = location_colors,
             Gender = gender_col, `Age group` = age_cols, `Primary mutation` = primary_mutations_col,
             `Secondary mutation` = second_mutations_col)
cols <- Map(function(x, cl) cl[levels(x)], anno, cols[names(anno)])
stopifnot("annotation value without a colour" = !any(vapply(anno, anyNA, logical(1))),
          !any(vapply(cols, anyNA, logical(1))))
oa <- R_("origin_all")
make_ht <- function() Heatmap(
  full, name = "KL divergence",
  col = colorRamp2(c(0, median(full), max(full)), c("white", "orange", "red")),
  top_annotation = do.call(HeatmapAnnotation, c(anno, list(col = cols, annotation_name_gp = gpar(fontsize = 11),
                                                          show_legend = c(FALSE, rep(TRUE, length(anno) - 1))))),
  show_row_names = FALSE, show_column_names = FALSE,
  column_title = sprintf("GloScope Divergence (All Cells, %d samples, %d patients)\nOrigin ANOSIM: R = %.2f, %s; one sample per patient per group: R = %.2f, %s",
                         nrow(full), n_distinct(a$Patient_nr), oa$R, fmt_p(oa$p), oa$R_1pt, fmt_p(oa$p_1pt)),
  column_title_gp = gpar(fontsize = 13, fontface = "bold"),
  heatmap_legend_param = list(title = "KL divergence", legend_height = unit(4, "cm")))
ht_file <- file.path(fig_dir, "SFig3b_kl_heatmap")
png(paste0(ht_file, ".png"), width = 14, height = 12, units = "in", res = 300)
draw(make_ht(), padding = unit(c(5, 5, 5, 15), "mm")); invisible(dev.off())
cairo_pdf(paste0(ht_file, ".pdf"), width = 14, height = 12)
draw(make_ht(), padding = unit(c(5, 5, 5, 15), "mm")); invisible(dev.off())

cat("Panels:", paste(list.files(fig_dir), collapse = ", "), "\nCSV:", file.path(out_dir, "anosim_manuscript_panels.csv"), "\n")
print(as.data.frame(csv |> select(key, cell_set, n_samples, n_patients, R, p, n_samples_1pt, R_1pt, p_1pt, R_march, p_march)),
      digits = 3)
