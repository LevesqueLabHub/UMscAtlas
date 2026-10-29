#!/usr/bin/env Rscript
# Manuscript panels for CytoTRACE2 (Fig. 2, Suppl. Fig. 4) and the Weizmann 3CA metaprograms (Fig. 3, Suppl. Fig. 6)
# rebuilt on the October 2026 atlas (Paul, 2026-10-15). Port of generate_revised_plots_cytotrace_mp.R (March 2026,
# January object); file names are the submitted figure IDs.
# Inputs: Outputs/CytoTRACE2_2026b (scores + atlas UMAP, run_cytotrace2_2026b.R), Outputs/Metaprograms_3CA_2026b/
#   ucell_3CA_scores_per_cell.qs2 (score_3CA_ucell_2026b.R; UCell, as the submitted panels - not the PLAID of report 02b),
#   new_data_2026b/full_meta_full.csv.gz. No October Seurat object is loaded.
# Design: theme_classic. Boxplots draw one point per tumour/tissue (orig.ident, >= MIN_CELLS cells) but test PATIENT
#   means (Patient_nr = unit of replication; two-sided Wilcoxon, unadjusted); groups with < 3 patients are drawn, not tested.
# Outputs: Figures_2026b/Manuscript_panels/<ID>.pdf + .png (300 dpi); Outputs/Manuscript_panels_2026b/tests_cytotrace_mp.csv
#   (patient-level tests, with the March p from pvalues_all_paper_plots.xlsx and its direction where one exists).
# The March direction needs the January per-cell values: cached once to january_inputs_for_direction.qs2 (loads three
#   ~27 GB January objects when the cache is missing).
# Run: Rscript manuscript_panels_2026b_cytotrace_mp.R
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(ggsignif); library(ggbeeswarm); library(scattermore)
  library(ComplexHeatmap); library(circlize)
})
source("config_colors_2026b.R")
A    <- ".."
FIG  <- file.path(A, "Figures_2026b/Manuscript_panels")
OUTT <- file.path(A, "Outputs/Manuscript_panels_2026b")
dir.create(FIG, recursive = TRUE, showWarnings = FALSE); dir.create(OUTT, recursive = TRUE, showWarnings = FALSE)
ORIGINS <- c("non-malignant", "primary", "metastasis")
MUTS <- c("BAP1", "SF3B1", "EIF1AX")
MIN_CELLS <- 10              # per orig.ident, as the March p_b1
UMAP_Q <- c(0.01, 0.99)      # colour cap of the metaprogram UMAPs
CAP <- "points = tumour/tissue (orig.ident); p = Wilcoxon on patient means"
potency_labels <- c("Differentiated", "Unipotent", "Oligopotent", "Multipotent", "Pluripotent", "Totipotent")
potency_colors <- setNames(c("#9E0142", "#F46D43", "#FEE08B", "#E6F598", "#66C2A5", "#5E4FA2"), potency_labels)
AGES <- c("<40", "40-60", "≥60", "unknown")   # manuscript bins; 40-60 = [40, 60)
# ageGroup_col (Young/Middle/Old) from config_colors_2026b.R; "unknown" is dark grey (added here: the palette has none, and
# gender_col's light grey cannot be told apart from the pale "Middle").
age_colors <- setNames(c(unname(ageGroup_col[c("Young", "Middle", "Old")]), "grey30"), AGES)
nfmt <- function(x) formatC(x, format = "d", big.mark = ",")   # format() would pad a vector to a common width
signif_label <- function(p) sprintf("%s  p %s",
  as.character(cut(p, c(-Inf, 0.001, 0.01, 0.05, Inf), c("***", "**", "*", "ns"))),
  ifelse(p < 2.2e-16, "< 2.2e-16", paste("=", formatC(p, digits = 2, format = "g"))))
wrap_cap <- function(x, width_in) paste(strwrap(x, floor(13 * width_in)), collapse = "\n")
thm <- theme_classic(base_size = 12) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5), plot.caption = element_text(hjust = 0, size = 8.5),
        plot.caption.position = "plot", plot.title.position = "plot")
save_gg <- function(p, id, w, h) {
  ggsave(file.path(FIG, paste0(id, ".pdf")), p, width = w, height = h, device = cairo_pdf)
  ggsave(file.path(FIG, paste0(id, ".png")), p, width = w, height = h, dpi = 300, bg = "white")
}
save_ht <- function(make, id, w, h) {   # build the Heatmap inside the device (a stored Heatmap segfaults a later ggsave)
  png(file.path(FIG, paste0(id, ".png")), width = w, height = h, units = "in", res = 300); draw(make(), padding = unit(c(3, 3, 3, 12), "mm")); invisible(dev.off())
  cairo_pdf(file.path(FIG, paste0(id, ".pdf")), width = w, height = h); draw(make(), padding = unit(c(3, 3, 3, 12), "mm")); invisible(dev.off())
}

# ---- data ---------------------------------------------------------------------------------------------------------
meta <- fread(file.path(A, "new_data_2026b/full_meta_full.csv.gz"),
              select = c("cell", "orig.ident", "Sample", "Patient_nr", "origin", "majority_celltype", "Malignant",
                         "primary_location", "location", "healthy_location", "secondary_mutation", "primary_mutation",
                         "Gender", "Age_atDiagnosisY", "TreatmentStage_ofProcessedSample"))
stopifnot(nrow(meta) == 349696, uniqueN(meta$orig.ident) == 56, uniqueN(meta$Sample) == 63, uniqueN(meta$Patient_nr) == 45)
ct <- qs2::qs_read(file.path(A, "Outputs/CytoTRACE2_2026b/cytotrace2_scores.qs2"), nthreads = 1)
d <- merge(meta, ct[, .(cell, CytoTRACE2_Score, CytoTRACE2_Potency)], by = "cell")
stopifnot(nrow(d) == nrow(meta), !anyNA(d$CytoTRACE2_Score), all(d$CytoTRACE2_Potency %in% potency_labels),
          all(d$origin %in% ORIGINS), all(d$Gender %in% names(gender_col)))
# Age is a patient attribute: non-malignant samples carry NA, their patient's tumour carries the age at diagnosis.
age <- unique(d[!is.na(Age_atDiagnosisY), .(Patient_nr, age = Age_atDiagnosisY)])
stopifnot("a patient has two ages" = !anyDuplicated(age$Patient_nr))
d[age, age := i.age, on = "Patient_nr"]
d[, age_class := factor(fcase(is.na(age), "unknown", age < 40, "<40", age < 60, "40-60", default = "≥60"), levels = AGES)]
n_age_unknown <- uniqueN(d[is.na(age), Patient_nr])
stopifnot("expected 8 patients without an age" = n_age_unknown == 8)
oi <- unique(d[, .(orig.ident, location, Gender, age_class, TreatmentStage_ofProcessedSample)])
stopifnot("orig.ident metadata not unique" = !anyDuplicated(oi$orig.ident), nrow(oi) == 56)

uc <- qs2::qs_read(file.path(A, "Outputs/Metaprograms_3CA_2026b/ucell_3CA_scores_per_cell.qs2"), nthreads = 1)
mps <- grep("^MP[0-9]+ ", names(uc), value = TRUE)
stopifnot(nrow(uc) == 200634, length(mps) == 67, all(uc$cell %in% d$cell),
          all(uc$majority_celltype %in% c("Melanoma cells", "Melanocytes")))
mp_of <- function(s) { m <- grep(paste0("^", s, " "), mps, value = TRUE); stopifnot(length(m) == 1); m }
umap <- fread(file.path(A, "Outputs/CytoTRACE2_2026b/atlas_umap.csv.gz"))
cat(sprintf("atlas: %s cells, %d orig.ident, %d Sample, %d patients; melanocytic (UCell): %s cells, %d orig.ident, %d patients\n",
            nfmt(nrow(d)), uniqueN(d$orig.ident), uniqueN(d$Sample), uniqueN(d$Patient_nr),
            nfmt(nrow(uc)), uniqueN(uc$orig.ident), uniqueN(uc$Patient_nr)))

# ---- units, tests, boxplot ----------------------------------------------------------------------------------------
# x: one row per cell with orig.ident, Patient_nr, the grouping column and the value v.
units <- function(x, grp) {
  pts <- x[, .(score = mean(v), n_cells = .N), by = c("orig.ident", "Patient_nr", grp)]
  stopifnot("an orig.ident sits in two groups" = !anyDuplicated(pts$orig.ident))
  drop <- pts[n_cells < MIN_CELLS, orig.ident]
  pts <- merge(pts[n_cells >= MIN_CELLS], oi, by = "orig.ident")
  pat <- x[orig.ident %in% pts$orig.ident, .(score = mean(v), n_cells = .N), by = c("Patient_nr", grp)]
  list(pts = pts, pat = pat, dropped = drop)
}
pat_tests <- function(u, grp, levels, panel, what) {
  n_by <- table(factor(u$pat[[grp]], levels = levels)); lv <- names(n_by)[n_by >= 3]
  if (length(lv) < 2) return(data.table())
  rbindlist(lapply(combn(lv, 2, simplify = FALSE), function(pr) {
    a <- u$pat[get(grp) == pr[1], score]; b <- u$pat[get(grp) == pr[2], score]
    data.table(panel = panel, score = what, group1 = pr[1], group2 = pr[2], n_patients1 = length(a), n_patients2 = length(b),
               median1 = median(a), median2 = median(b), p = suppressWarnings(wilcox.test(a, b)$p.value))
  }))
}
box_panel <- function(u, grp, levels, tests, title, ylab, fill = NULL, colour = NULL, colour_values = NULL,
                      colour_name = NULL, shape = NULL, size_cells = FALSE, note = NULL, width = 7.5) {
  pts <- u$pts[get(grp) %in% levels]; pat <- u$pat[get(grp) %in% levels]
  levels <- levels[levels %in% pts[[grp]]]
  stopifnot(nrow(pts) > 0)
  pts[, x := factor(get(grp), levels = levels)]
  n_pt <- table(pts$x); n_pa <- table(factor(pat[[grp]], levels = levels)); small <- names(n_pa)[n_pa < 3]
  xlab <- setNames(sprintf("%s\nn = %d (%d pat.)", levels, as.integer(n_pt), as.integer(n_pa)), levels)
  cap <- paste(c(paste0(CAP, "."), sprintf("One point per orig.ident with >= %d cells; two-sided, unadjusted.", MIN_CELLS),
                 "n = orig.ident (patients).",
                 if (length(small)) sprintf("Not tested (< 3 patients): %s.", paste(small, collapse = ", ")), note), collapse = " ")
  tt <- tests[group1 %in% levels & group2 %in% levels]
  m <- list(); if (!is.null(colour)) m$colour <- as.name(colour)
  if (size_cells) m$size <- as.name("n_cells"); if (!is.null(shape)) m$shape <- as.name(shape)
  fixed <- c(if (is.null(colour)) list(colour = "black"), if (!size_cells) list(size = 2.2))
  p <- ggplot(pts, aes(x, score)) +
    (if (is.null(fill)) geom_boxplot(fill = "grey85", alpha = 0.7, outlier.shape = NA, width = 0.6)
     else geom_boxplot(aes(fill = x), alpha = 0.45, outlier.shape = NA, width = 0.6, show.legend = FALSE)) +
    do.call(geom_quasirandom, c(list(mapping = do.call(aes, m), width = 0.22, alpha = 0.85), fixed)) +
    scale_x_discrete(labels = xlab) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.06 + 0.11 * nrow(tt)))) +
    labs(title = title, x = NULL, y = ylab, caption = wrap_cap(cap, width)) + thm
  if (nrow(tt)) p <- p + geom_signif(comparisons = unname(Map(c, tt$group1, tt$group2)), annotations = signif_label(tt$p),
                                     step_increase = 0.11, tip_length = 0.01, textsize = 3.2, vjust = -0.2)
  if (!is.null(fill)) p <- p + scale_fill_manual(values = fill)
  if (!is.null(colour)) p <- p + scale_colour_manual(values = colour_values, name = colour_name,
                                                     breaks = intersect(names(colour_values), unique(as.character(pts[[colour]]))))
  if (size_cells) p <- p + scale_size_continuous(name = "Cells", range = c(1.5, 6), breaks = c(100, 1000, 5000, 10000), labels = nfmt)
  if (!is.null(shape)) p <- p + scale_shape_manual(values = c(primary = 16, metastasis = 17), name = NULL)
  p
}
tests <- list()

# ---- Fig. 2c: potency on the atlas UMAP ----------------------------------------------------------------------------
um <- merge(umap, d[, .(cell, CytoTRACE2_Potency)], by = "cell")
stopifnot("UMAP misses cells" = nrow(um) == nrow(d))
um[, pot := factor(CytoTRACE2_Potency, levels = potency_labels)]
setorder(um, pot)   # higher potency drawn last (on top), as DimPlot(order = TRUE) in the March panel
n_pot <- table(um$pot)
p <- ggplot(um, aes(UMAP_1, UMAP_2, colour = pot)) + geom_scattermore(pointsize = 2.8, pixels = c(2400, 2400)) +
  scale_colour_manual(values = potency_colors, drop = FALSE, name = "CytoTRACE2 potency",
                      labels = sprintf("%s (%s)", names(n_pot), nfmt(as.integer(n_pot)))) +
  guides(colour = guide_legend(override.aes = list(size = 4))) +
  labs(title = "CytoTRACE2 potency", x = "UMAP 1", y = "UMAP 2",
       caption = sprintf("All %s cells of the October 2026 atlas (63 libraries, 45 patients).", nfmt(nrow(um)))) +
  thm + theme(aspect.ratio = 1)
save_gg(p, "Fig2c_cytotrace_potency_umap", 9.5, 7.5)

# ---- Fig. 2d: score by cell type -----------------------------------------------------------------------------------
ord <- d[, .(m = median(CytoTRACE2_Score)), by = majority_celltype][order(-m), majority_celltype]
stopifnot("cell type without a colour" = all(ord %in% names(majority_celltype_col)))
n_ct <- d[, .N, by = majority_celltype][, setNames(N, majority_celltype)]
p <- ggplot(d, aes(factor(majority_celltype, levels = ord), CytoTRACE2_Score, fill = majority_celltype)) +
  geom_violin(scale = "width", linewidth = 0.3) +
  stat_summary(fun = median, fun.min = median, fun.max = median, geom = "crossbar", width = 0.5, linewidth = 0.3) +
  scale_fill_manual(values = majority_celltype_col, guide = "none") +
  scale_x_discrete(labels = sprintf("%s (%s)", ord, nfmt(n_ct[ord]))) +
  labs(title = "CytoTRACE2 score by cell type", x = NULL, y = "CytoTRACE2 score",
       caption = sprintf("All %s cells, cell level; bar = median. Descriptive, no test.", nfmt(nrow(d)))) +
  thm + theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.margin = margin(5.5, 5.5, 5.5, 60))
save_gg(p, "Fig2d_cytotrace_violin_celltype", 11, 7)

# ---- Fig. 2e + Suppl. Fig. 4a/b: melanocytic cells by origin -------------------------------------------------------
mel <- d[Malignant == "Malignant" | majority_celltype == "Melanocytes", .(orig.ident, Patient_nr, origin, v = CytoTRACE2_Score)]
u_or <- units(mel, "origin")
t_or <- pat_tests(u_or, "origin", ORIGINS, "Fig2e_cytotrace_by_origin_location", "CytoTRACE2")
n_shared <- length(intersect(u_or$pat[origin == "non-malignant", Patient_nr], u_or$pat[origin == "primary", Patient_nr]))
note_or <- sprintf("Melanoma cells of tumours, melanocytes of non-malignant tissue. %d patients give both a non-malignant and a primary unit.", n_shared)
ylab_or <- "Mean CytoTRACE2 score"
tests$f2e <- t_or
save_gg(box_panel(u_or, "origin", ORIGINS, t_or, "CytoTRACE2 by disease stage", ylab_or, colour = "location",
                  colour_values = location_colors, colour_name = "Location", size_cells = TRUE, note = note_or),
        "Fig2e_cytotrace_by_origin_location", 7.5, 6.5)
tests$s4a <- copy(t_or)[, panel := "SFig4a_cytotrace_origin_age"]
save_gg(box_panel(u_or, "origin", ORIGINS, t_or, "CytoTRACE2 by disease stage, age class", ylab_or, colour = "age_class",
                  colour_values = age_colors, colour_name = "Age at diagnosis", size_cells = TRUE,
                  note = paste(note_or, sprintf("Age class per patient; %d patients have no age (unknown).", n_age_unknown))),
        "SFig4a_cytotrace_origin_age", 7.5, 6.5)
tests$s4b <- copy(t_or)[, panel := "SFig4b_cytotrace_origin_gender"]
save_gg(box_panel(u_or, "origin", ORIGINS, t_or, "CytoTRACE2 by disease stage, gender", ylab_or, colour = "Gender",
                  colour_values = gender_col, colour_name = "Gender", size_cells = TRUE, note = note_or),
        "SFig4b_cytotrace_origin_gender", 7.5, 6.5)

# Reviewer (revision doc L172, 2026-10-28): the same per-patient test for gender and age as for origin. Within each
# disease stage, since origin moves CytoTRACE2 most. Grouping column renamed (g) so units() can still merge oi.
by_stage <- function(col, levels, test_lv, fill, id, title, note) {
  tt <- list(); nlv <- integer()
  loc_all <- intersect(names(location_colors), unique(oi$location))
  ps <- lapply(ORIGINS, function(o) {
    x <- d[(Malignant == "Malignant" | majority_celltype == "Melanocytes") & origin == o,
           .(orig.ident, Patient_nr, g = as.character(get(col)), v = CytoTRACE2_Score)]
    u <- units(x, "g")
    tt[[o]] <<- pat_tests(u, "g", test_lv, paste0(id, ":", o), "CytoTRACE2")
    nlv[o] <<- sum(levels %in% u$pts$g)
    p <- box_panel(u, "g", levels, tt[[o]], o, if (o == ORIGINS[1]) ylab_or else NULL, fill = fill,
                   colour = "location", colour_values = location_colors, colour_name = "Location")
    p$layers[[2]]$show.legend <- TRUE   # draw every location key in every panel, so patchwork collects one legend
    p +
      labs(caption = NULL) + theme(axis.text.x = element_text(size = 9)) + guides(fill = "none") +
      suppressMessages(scale_colour_manual(values = location_colors, name = "Location", drop = FALSE, limits = loc_all))   # one shared legend
  })
  cap <- sprintf("%s. One point per orig.ident with >= %d cells, coloured by location; Wilcoxon on patient means within each stage, two-sided, unadjusted. n = orig.ident (patients). %s",
                 CAP, MIN_CELLS, note)
  save_gg(patchwork::wrap_plots(ps, nrow = 1, widths = nlv, guides = "collect") +
            patchwork::plot_annotation(title = title, caption = wrap_cap(cap, 13), theme = thm), id, 14, 6.5)
  rbindlist(tt)
}
# Gender and location are entangled in the primaries (iris tumours, lower CytoTRACE2, are mostly from women):
# the same patient means in a linear model with location.
pg <- d[Malignant == "Malignant" & origin == "primary" & primary_location %in% c("choroid", "iris", "ciliary body"),
        .(s = mean(CytoTRACE2_Score), n = .N), by = .(Patient_nr, Gender, primary_location)][n >= MIN_CELLS]
stopifnot("a patient has two primary rows" = !anyDuplicated(pg$Patient_nr))
fit <- coef(summary(lm(s ~ Gender + primary_location, pg)))["GenderM", ]
n_iris <- pg[primary_location == "iris", table(factor(Gender, c("F", "M")))]
tests$s4b_test <- by_stage("Gender", c("F", "M"), c("F", "M"), gender_col, "SFig4b_cytotrace_gender_by_stage",
                           "CytoTRACE2 by gender within disease stage",
                           paste(note_or, sprintf("Primaries: %d of %d women and %d of %d men have an iris tumour; adjusted for location (linear model on patient means, %d patients, choroid / iris / ciliary body) the male - female difference is %.3f, p = %.2g.",
                                                  n_iris[["F"]], pg[Gender == "F", .N], n_iris[["M"]], pg[Gender == "M", .N],
                                                  nrow(pg), fit[["Estimate"]], fit[["Pr(>|t|)"]])))
tests$s4a_test <- by_stage("age_class", AGES, AGES[1:3], age_colors, "SFig4a_cytotrace_age_by_stage",
                           "CytoTRACE2 by age class within disease stage",
                           paste(note_or, sprintf("Unknown age (%d patients) shown, not tested.", n_age_unknown)))

# ---- Reviewer L169: fibroblast CytoTRACE2 by site (January object: the high-potency fibroblasts came from iris) ----------
SITES <- c("iris", "ciliary body", "choroid", "metastasis")
fib <- d[majority_celltype == "Fibroblasts"]
fib[, site := fifelse(origin == "metastasis", "metastasis", sub("^non-malignant ", "", location))]
stopifnot(nrow(fib) > 1e4, all(fib$site %in% c(SITES, "ciliochoroidal")))
u_fib <- units(fib[, .(orig.ident, Patient_nr, site, v = CytoTRACE2_Score)], "site")
tests$l169 <- pat_tests(u_fib, "site", SITES, "Rev_L169_fibroblast_cytotrace_site", "CytoTRACE2")
hi <- fib[CytoTRACE2_Score > 0.5, .N, by = site][, setNames(N, site)]
fib_note <- sprintf("Fibroblasts of every sample, non-malignant tissue and tumours; eye sites pool non-malignant and primary. %s fibroblasts, %s from iris samples. Of the %s fibroblasts with CytoTRACE2 > 0.5, %.0f%% are from iris, %.0f%% ciliary body, %.0f%% choroid, %.0f%% metastases.",
                    nfmt(nrow(fib)), nfmt(fib[site == "iris", .N]), nfmt(sum(hi)),
                    100 * hi[["iris"]] / sum(hi), 100 * hi[["ciliary body"]] / sum(hi), 100 * hi[["choroid"]] / sum(hi), 100 * hi[["metastasis"]] / sum(hi))
save_gg(box_panel(u_fib, "site", SITES, tests$l169, "CytoTRACE2 in fibroblasts by site", "Mean CytoTRACE2 score, fibroblasts",
                  fill = c(primary_location_colors[c("iris", "ciliary body", "choroid")], origin_colors["metastasis"]),
                  colour = "location", colour_values = location_colors, colour_name = "Location", size_cells = TRUE,
                  note = fib_note, width = 8.5),
        "Rev_L169_fibroblast_cytotrace_site", 8.5, 6.5)

# ---- Fig. 2i: pUM malignant cells by intra-ocular location ---------------------------------------------------------
LOC_2I <- c("choroid", "iris", "ciliary body")   # ciliochoroidal (1 patient) dropped, Victoria 2026-10-16
u_loc <- units(d[Malignant == "Malignant" & origin == "primary", .(orig.ident, Patient_nr, primary_location, v = CytoTRACE2_Score)],
               "primary_location")
tests$f2i <- pat_tests(u_loc, "primary_location", LOC_2I, "Fig2i_cytotrace_pUM_location", "CytoTRACE2")
save_gg(box_panel(u_loc, "primary_location", LOC_2I, tests$f2i, "CytoTRACE2 in primary UM by location",
                  "Mean CytoTRACE2 score, malignant cells", fill = primary_location_colors, width = 7),
        "Fig2i_cytotrace_pUM_location", 7, 6.5)

# ---- Suppl. Fig. 4c: malignant cells by secondary mutation ---------------------------------------------------------
u_mut <- units(d[Malignant == "Malignant" & secondary_mutation %in% MUTS,
                 .(orig.ident, Patient_nr, secondary_mutation, v = CytoTRACE2_Score)], "secondary_mutation")
u_mut$pts <- merge(u_mut$pts, unique(d[, .(orig.ident, origin)]), by = "orig.ident")
tests$s4c <- pat_tests(u_mut, "secondary_mutation", MUTS, "SFig4c_cytotrace_secondary_mutation_gender", "CytoTRACE2")
n_bap1_met <- u_mut$pts[secondary_mutation == "BAP1" & origin == "metastasis", .N]
save_gg(box_panel(u_mut, "secondary_mutation", MUTS, tests$s4c, "CytoTRACE2 by secondary mutation",
                  "Mean CytoTRACE2 score, malignant cells", fill = second_mutations_col, colour = "Gender",
                  colour_values = gender_col, colour_name = "Gender", shape = "origin",
                  note = sprintf("Primary and metastatic tumours, as the March panel: %d of the BAP1 points are metastases.", n_bap1_met)),
        "SFig4c_cytotrace_secondary_mutation_gender", 7.5, 6.5)

# ---- Fig. 3a: top 30 programs, cell-level correlation --------------------------------------------------------------
mp_mean <- sort(colMeans(uc[, ..mps]), decreasing = TRUE)
top30 <- names(mp_mean)[1:30]
cm <- cor(as.matrix(uc[, ..top30]))
dimnames(cm) <- list(sub(" .*", "", top30), sub(" .*", "", top30))
save_ht(function() Heatmap(cm, name = "Pearson r", col = colorRamp2(c(-1, 0, 1), c("blue", "white", "red")),
                           row_names_gp = gpar(fontsize = 9), column_names_gp = gpar(fontsize = 9),
                           column_title = sprintf("Top 30 metaprograms by mean UCell, cell-level correlation\n(%s melanocytic cells)", nfmt(nrow(uc))),
                           column_title_gp = gpar(fontsize = 11, fontface = "bold")),
        "Fig3a_mp_top30_correlation", 8.5, 8)

# ---- Fig. 3c/e/g UMAPs + 3d/f/h boxplots by origin -----------------------------------------------------------------
umm <- merge(umap, uc[, c("cell", mps), with = FALSE], by = "cell")
stopifnot("UMAP misses melanocytic cells" = nrow(umm) == nrow(uc))
set.seed(14); umm <- umm[sample(.N)]   # random draw order: no group painted over the others
n_mel <- uc[majority_celltype == "Melanocytes", .N]; n_mal <- uc[majority_celltype == "Melanoma cells", .N]
cells_txt <- sprintf("%s melanocytic cells (%s melanoma cells, %s melanocytes)", nfmt(nrow(uc)), nfmt(n_mal), nfmt(n_mel))
for (k in list(c("MP45", "Fig3c_umap_MP45", "Fig3d_box_MP45"), c("MP17", "Fig3e_umap_MP17", "Fig3f_box_MP17"),
               c("MP9", "Fig3g_umap_MP9", "Fig3h_box_MP9"))) {
  m <- mp_of(k[1])
  lim <- quantile(umm[[m]], UMAP_Q); if (diff(lim) <= 0) lim <- range(umm[[m]])
  p <- ggplot(umm, aes(UMAP_1, UMAP_2, colour = .data[[m]])) + geom_scattermore(pointsize = 2.8, pixels = c(2400, 2400)) +
    scale_colour_viridis_c(option = "magma", limits = lim, oob = scales::squish, name = "UCell") +
    labs(title = m, x = "UMAP 1", y = "UMAP 2",
         caption = wrap_cap(sprintf("Atlas UMAP, %s. Colour capped at percentiles %g and %g (%.3f-%.3f); cells outside are drawn at the end colour.",
                                    cells_txt, 100 * UMAP_Q[1], 100 * UMAP_Q[2], lim[1], lim[2]), 7.5)) +
    thm + theme(aspect.ratio = 1)
  save_gg(p, k[2], 7.5, 7)
  u <- units(uc[, .(orig.ident, Patient_nr, origin, v = get(m))], "origin")
  tests[[k[3]]] <- pat_tests(u, "origin", ORIGINS, k[3], m)
  save_gg(box_panel(u, "origin", ORIGINS, tests[[k[3]]], m, "Mean UCell score", fill = origin_colors,
                    note = "Melanoma cells of tumours, melanocytes of non-malignant tissue.", width = 6.5), k[3], 6.5, 6.5)
}

# ---- Fig. 3i/j: interferon/MHC-II programs, one boxplot per program --------------------------------------------------
# Paul, 2026-10-16: boxplots per tumour/tissue replace the submitted Seurat DotPlot, whose dot size (% cells with
# UCell > 0) was 92-100 % in every group and whose colour was a z-score across only 3-5 group means.
ifn <- vapply(c("MP22", "MP23", "MP24"), mp_of, character(1))
no_tests <- data.table(group1 = character(), group2 = character(), p = numeric())
ifn_boxes <- function(x, grp, levels, test, fill, id, title, cap, w, h, flip = FALSE) {
  tt <- list()
  ps <- lapply(seq_along(ifn), function(i) {
    u <- units(x[, .(orig.ident, Patient_nr, g = get(grp), v = get(ifn[i]))], "g")
    tt[[i]] <<- if (test) pat_tests(u, "g", levels, id, ifn[i]) else no_tests
    p <- box_panel(u, "g", if (flip) rev(levels) else levels, tt[[i]], ifn[i], "Mean UCell score", fill = fill) +
      labs(caption = NULL)
    if (flip) p <- p + coord_flip() + if (i > 1) theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
    p
  })
  save_gg(patchwork::wrap_plots(ps, nrow = 1) + patchwork::plot_annotation(title = title, caption = wrap_cap(cap, w), theme = thm),
          id, w, h)
  if (test) rbindlist(tt)
}
tests$f3i <- ifn_boxes(uc, "origin", ORIGINS, TRUE, origin_colors, "Fig3i_box_MP22_24_origin", "Interferon/MHC-II programs by disease stage",
                       sprintf("%s. %s. One point per orig.ident with >= %d cells; two-sided, unadjusted. n = orig.ident (patients). Melanoma cells of tumours, melanocytes of non-malignant tissue.",
                               CAP, cells_txt, MIN_CELLS), 14.5, 6.5)
met <- uc[origin == "metastasis"]
n_tr <- met[, .(n_pat = uniqueN(Patient_nr), n_cells = .N), by = TreatmentStage_ofProcessedSample][order(-n_pat, -n_cells)]
ifn_boxes(met, "TreatmentStage_ofProcessedSample", n_tr$TreatmentStage_ofProcessedSample, FALSE, NULL,
          "Fig3j_box_MP22_24_treatment", "Interferon/MHC-II programs by treatment, metastases",
          sprintf("points = tumour/tissue (orig.ident). %s melanoma cells of %d metastases (%d patients; %s each in two arms). n = orig.ident (patients). %d-%d patients per arm: descriptive, no test.",
                  nfmt(nrow(met)), uniqueN(met$orig.ident), uniqueN(met$Patient_nr),
                  paste(met[, uniqueN(TreatmentStage_ofProcessedSample), by = Patient_nr][V1 > 1, Patient_nr], collapse = " and "),
                  min(n_tr$n_pat), max(n_tr$n_pat)),
          14, 5.5, flip = TRUE)

# Victoria (2026-10-16) wants the submitted DotPlot design back: Seurat DotPlot semantics on the UCell columns
# (colour = per-program z of log1p(mean(expm1(UCell))) across groups, clipped +-2.5; size = % cells with UCell > 0).
# The boxplots above stay in the panel set as the patient-level view.
dot_panel <- function(x, grp, levels, labels, title, id, w, h) {
  long <- melt(x[, c(grp, ifn), with = FALSE], id.vars = grp, variable.name = "MP", value.name = "s")
  dd <- long[, .(avg = mean(expm1(s)), pct = 100 * mean(s > 0)), by = c(grp, "MP")]
  dd[, z := pmin(pmax(as.numeric(scale(log1p(avg))), -2.5), 2.5), by = MP]
  dd[, g := factor(get(grp), levels = levels)][, MP := factor(MP, levels = ifn)]
  stopifnot(!anyNA(dd$g))
  p <- ggplot(dd, aes(g, MP, colour = z, size = pct)) + geom_point() +
    scale_colour_gradient(low = "lightgrey", high = "blue", name = "Average Expression") +
    scale_radius(range = c(1, 6), name = "Percent Expressed") +
    scale_x_discrete(labels = labels) + labs(title = title, x = NULL, y = NULL) + thm +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  save_gg(p, id, w, h)
}
dot_panel(uc, "origin", ORIGINS, setNames(ORIGINS, ORIGINS), NULL, "Fig3i_dotplot_MP22_24_origin", 6.5, 4.5)
tr_lv <- sort(unique(met$TreatmentStage_ofProcessedSample))
dot_panel(met, "TreatmentStage_ofProcessedSample", tr_lv, setNames(sub(" \\(", "\n(", tr_lv), tr_lv), NULL,
          "Fig3j_dotplot_MP22_24_treatment", 8, 5)

# ---- Suppl. Fig. 6a: top 20 programs by mean UCell ------------------------------------------------------------------
b <- data.table(MP = names(mp_mean)[1:20], mean = unname(mp_mean[1:20]))
p <- ggplot(b, aes(mean, reorder(MP, mean))) + geom_col(fill = "steelblue", width = 0.75) +
  geom_text(aes(label = sprintf("%.3f", mean)), hjust = -0.15, size = 3.2) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(title = "Top 20 metaprograms by mean UCell score", x = "Mean UCell score", y = NULL,
       caption = wrap_cap(sprintf("All %s; cell-level mean.", cells_txt), 8)) + thm
save_gg(p, "SFig6a_mp_top20_mean_ucell", 8, 7)

# ---- Suppl. Fig. 6b: MP45 in primary UM by location -----------------------------------------------------------------
LOC_6B <- c("choroid", "ciliary body", "iris")
m45 <- mp_of("MP45")
u6 <- units(uc[Malignant == "Malignant" & origin == "primary", .(orig.ident, Patient_nr, primary_location, v = get(m45))],
            "primary_location")
tests$s6b <- pat_tests(u6, "primary_location", LOC_6B, "SFig6b_box_MP45_pUM_location", m45)
save_gg(box_panel(u6, "primary_location", LOC_6B, tests$s6b, paste(m45, "in primary UM"), "Mean UCell score, malignant cells",
                  fill = primary_location_colors, note = "Ciliochoroidal (1 patient) not shown, as in the March panel.", width = 7),
        "SFig6b_box_MP45_pUM_location", 7, 6.5)

# ---- tests table, with the March p and its direction ----------------------------------------------------------------
JAN <- file.path(OUTT, "january_inputs_for_direction.qs2")
if (!file.exists(JAN)) {   # one-off, ~10 min: per-cell values of the three January objects the March tests used
  suppressPackageStartupMessages(library(SeuratObject))
  o <- qs2::qs_read(file.path(A, "UM_atlas_cytotrace2.qs2"), nthreads = 1)
  ct2 <- setNames(o$CytoTRACE2_Score, colnames(o)); rm(o); invisible(gc())
  o <- qs2::qs_read(file.path(A, "13_final_Dataset_reclustered_28012025.qs2"), nthreads = 1)
  jc <- data.table(cell = colnames(o), o@meta.data[, c("Sample", "origin", "Malignant", "majority_celltype",
                                                       "secondary_mutation", "primary_location")])
  rm(o); invisible(gc()); jc[, CytoTRACE2_Score := unname(ct2[cell])]
  o <- qs2::qs_read(file.path(A, "melanocytic_metaprograms.qs2"), nthreads = 1)
  jm <- data.table(cell = colnames(o), o@meta.data[, c("origin", "Malignant", "LocationPrimaryTumor",
                                                       grep("^MP(45|17|9) ", colnames(o@meta.data), value = TRUE))])
  rm(o); invisible(gc())
  qs2::qs_save(list(ct = jc, mp = jm), JAN, nthreads = 1)
}
jan <- qs2::qs_read(JAN, nthreads = 1); jc <- jan$ct; jm <- jan$mp
jmp <- function(s) { m <- grep(paste0("^", s, " "), names(jm), value = TRUE); stopifnot(length(m) == 1); m }
# What each March test compared (export_pvalues.R): B1 = sample means over ALL cells (Sample >= 10 cells), the rest cell level.
old_vals <- list(
  B1 = jc[, .(v = mean(CytoTRACE2_Score, na.rm = TRUE), n = .N), by = .(Sample, g = origin)][n >= 10],
  B2 = jc[Malignant == "Malignant" & secondary_mutation %in% MUTS, .(v = CytoTRACE2_Score, g = secondary_mutation)],
  D3 = jc[(Malignant == "Malignant" | majority_celltype == "Melanocytes") & origin == "primary" &
            sub("chiliochoroidal", "ciliochoroidal", primary_location) %in% c("iris", "choroid", "ciliary body"),
          .(v = CytoTRACE2_Score, g = primary_location)],
  B4 = jm[origin == "primary" & Malignant == "Malignant" & LocationPrimaryTumor %in% c("iris", "choroid", "ciliary body"),
          .(v = get(jmp("MP45")), g = LocationPrimaryTumor)],
  C1 = jm[Malignant == "Malignant", .(v = get(jmp("MP45")), g = origin)],
  C2 = jm[Malignant == "Malignant", .(v = get(jmp("MP17")), g = origin)],
  C3 = jm[Malignant == "Malignant", .(v = get(jmp("MP9")), g = origin)])
march_panel <- c(B1 = "Fig2e_cytotrace_by_origin_location", B2 = "SFig4c_cytotrace_secondary_mutation_gender",
                 D3 = "Fig2i_cytotrace_pUM_location", B4 = "SFig6b_box_MP45_pUM_location",
                 C1 = "Fig3d_box_MP45", C2 = "Fig3f_box_MP17", C3 = "Fig3h_box_MP9")
march_unit <- c(B1 = "sample means of all cells (Sample >= 10 cells)", B2 = "cells", D3 = "cells", B4 = "cells",
                C1 = "cells", C2 = "cells", C3 = "cells")
dir_of <- function(g1, g2, m1, m2) ifelse(m2 > m1, paste(g2, ">", g1), paste(g1, ">", g2))
pair_key <- function(a, b) vapply(seq_along(a), function(i) paste(sort(c(a[i], b[i])), collapse = " | "), character(1))
xl <- as.data.table(readxl::read_excel(file.path(A, "Revised_Plots_20260325/pvalues_all_paper_plots.xlsx"), sheet = "All_comparisons"))
xl <- xl[Figure %in% names(old_vals)]
old <- rbindlist(lapply(seq_len(nrow(xl)), function(i) {
  r <- xl[i]; v <- old_vals[[r$Figure]]
  a <- na.omit(v[g == r$Group1, v]); b <- na.omit(v[g == r$Group2, v])
  g1 <- sub("^healthy$", "non-malignant", r$Group1); g2 <- sub("^healthy$", "non-malignant", r$Group2)
  data.table(panel = march_panel[[r$Figure]], pair = pair_key(g1, g2), march_figure = r$Figure,
             march_test = paste0(r$Description, "; unit = ", march_unit[[r$Figure]]),
             march_n1 = r$N_Group1, march_n2 = r$N_Group2, march_p = r$p_value,
             march_p_recomputed = suppressWarnings(wilcox.test(a, b, exact = FALSE)$p.value),
             n_reproduced = length(a) == r$N_Group1 && length(b) == r$N_Group2,
             march_direction = dir_of(g1, g2, median(a), median(b)))
}))
print(old[, .(march_figure, pair, march_n1, march_n2, march_p, march_p_recomputed, n_reproduced, march_direction)])
stopifnot("January values do not reproduce the March n" = all(old$n_reproduced),
          "January values do not reproduce the March p" = isTRUE(all.equal(old$march_p, old$march_p_recomputed, tolerance = 1e-4)))
tt <- rbindlist(tests)
tt[, `:=`(comparison = paste(group1, "vs", group2), pair = pair_key(group1, group2),
          direction = dir_of(group1, group2, median1, median2))]
out <- merge(tt, old[, !"march_p_recomputed"], by = c("panel", "pair"), all.x = TRUE, sort = FALSE)
out[, direction_changed := ifelse(is.na(march_direction), NA, march_direction != direction)]
# The xlsx B1 test used ALL cells, the March B1 plot melanoma cells + melanocytes (the Fig2e cell set): give that one too.
b1p <- jc[Malignant == "Malignant" | majority_celltype == "Melanocytes", .(v = mean(CytoTRACE2_Score, na.rm = TRUE), n = .N),
          by = .(Sample, g = sub("^healthy$", "non-malignant", origin))][n >= 10]
out[panel == "Fig2e_cytotrace_by_origin_location", march_plot_cells_check := mapply(function(a, b) {
  x <- b1p[g == a, v]; y <- b1p[g == b, v]
  sprintf("March plot cells (melanoma + melanocytes, %d vs %d Samples): %s, p = %.2g", length(x), length(y),
          dir_of(a, b, median(x), median(y)), suppressWarnings(wilcox.test(x, y, exact = FALSE)$p.value))
}, sub(" vs .*", "", comparison), sub(".* vs ", "", comparison))]
out <- out[, .(panel, score, comparison, n_patients1, n_patients2, median_patient_mean1 = median1,
               median_patient_mean2 = median2, p_patient_level = p, direction, march_figure, march_test,
               march_n1, march_n2, march_p, march_direction, direction_changed, march_plot_cells_check)]
fwrite(out, file.path(OUTT, "tests_cytotrace_mp.csv"))
print(out[, .(panel, comparison, n_patients1, n_patients2, p = signif(p_patient_level, 2), direction,
              march_p = signif(march_p, 2), march_direction, direction_changed)])
cat("done:", length(list.files(FIG, "^(Fig2|Fig3|SFig4|SFig6).*\\.png$")), "PNGs in", FIG, "\n")
