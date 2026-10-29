#!/usr/bin/env Rscript
# 23b: the two R2 L157 plots of 23_pertumor_cytotrace_cnv.R (January object) re-run on the October 2026 atlas
# (Aizhan, revision doc comment, 2026-10-27: "could you re run those two plots with newest dataset?").
#   F1  CytoTRACE2 per tumour (intratumour spread), malignant cells
#   F2  within-tumour Spearman, number of aberrant chromosomes (has_cnv_chr1-22) vs CytoTRACE2, raw and
#       adjusted for genes detected per cell (log1p nFeature_RNA)
#   F3  the other reading of "genomic complexity" (Aizhan, same doc: CytoTRACE1-style "many unique genes"):
#       within-tumour Spearman, genes detected vs CytoTRACE2
# Inputs: new_data_2026b/full_meta_full.csv.gz (iCNV calls, Malignant, nFeature_RNA) and
# Outputs/CytoTRACE2_2026b/cytotrace2_scores.qs2. Unit = tumour (orig.ident) with >= 50 malignant cells; the
# between-tumour line uses one value per patient. Run: Rscript 23b_pertumor_cytotrace_cnv_2026b.R (login node, ~2 min)
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
A   <- ".."
OUT <- file.path(A, "PerTumor_CytoTRACE2_CNV_2026b")
dir.create(OUT, showWarnings = FALSE)
MIN_CELLS <- 50
chrs <- paste0("chr", 1:22)

m <- fread(file.path(A, "new_data_2026b/full_meta_full.csv.gz"),
           select = c("cell", "orig.ident", "Patient_nr", "origin", "Malignant", "nFeature_RNA", paste0("has_cnv_", chrs)))
ct <- qs2::qs_read(file.path(A, "Outputs/CytoTRACE2_2026b/cytotrace2_scores.qs2"), nthreads = 1)
m <- merge(m[Malignant == "Malignant" & origin %in% c("primary", "metastasis")], ct[, .(cell, ct2 = CytoTRACE2_Score)], by = "cell")
stopifnot(nrow(m) > 1e5, !anyNA(m$ct2))
hc <- paste0("has_cnv_", chrs)
m[, n_aberrant := rowSums(as.matrix(.SD[, lapply(.SD, function(x) x %in% c(TRUE, "TRUE", 1, "1"))])), .SDcols = hc]
stopifnot("no aberrant chromosome called at all - check the has_cnv_ encoding" = m[, mean(n_aberrant > 0)] > 0.5)
m <- m[, if (.N >= MIN_CELLS) .SD, by = orig.ident]
m[, origin := factor(origin, c("primary", "metastasis"))]
cat("malignant cells:", nrow(m), "| tumours:", uniqueN(m$orig.ident), "| patients:", uniqueN(m$Patient_nr), "\n")

res_d <- function(y, z) resid(lm(y ~ log1p(z)))
per <- m[, {
  ok <- sd(n_aberrant) > 0
  .(n = .N, ct2_med = median(ct2), naberr_med = median(n_aberrant),
    rho_cnv = if (ok) cor(n_aberrant, ct2, method = "spearman") else NA_real_,
    rho_cnv_adj = if (ok) cor(res_d(n_aberrant, nFeature_RNA), res_d(ct2, nFeature_RNA), method = "spearman") else NA_real_,
    rho_genes = cor(nFeature_RNA, ct2, method = "spearman"))
}, by = .(orig.ident, Patient_nr, origin)]
fwrite(per, file.path(OUT, "per_tumour_summary.csv"))

# between tumours: one value per patient (patient means), primaries and metastases separately
pt <- m[, .(ct2 = mean(ct2), n_aberrant = mean(n_aberrant)), by = .(Patient_nr, origin)]
btw <- pt[, { t <- suppressWarnings(cor.test(ct2, n_aberrant, method = "spearman")); .(patients = .N, rho = t$estimate, p = t$p.value) }, by = origin]
fwrite(btw, file.path(OUT, "between_patient_cnv_vs_cytotrace2.csv"))
print(per[, .(tumours = .N, rho_cnv = median(rho_cnv, na.rm = TRUE), rho_cnv_adj = median(rho_cnv_adj, na.rm = TRUE),
              rho_genes = median(rho_genes))]); print(btw)

sign_p <- function(x) { x <- x[!is.na(x)]; binom.test(sum(x > 0), length(x))$p.value }
fmt <- function(x) formatC(x, digits = 2, format = "g")
org_cols <- c(primary = "#4363D8", metastasis = "#911D51")
save2 <- function(p, f, w, h) { ggsave(file.path(OUT, paste0(f, ".png")), p, width = w, height = h, dpi = 300, bg = "white")
  ggsave(file.path(OUT, paste0(f, ".pdf")), p, width = w, height = h, device = cairo_pdf) }

# F1 violin, tumours ordered by median CytoTRACE2
m[, tum := factor(orig.ident, per[order(ct2_med), orig.ident])]
p1 <- ggplot(m, aes(tum, ct2, fill = origin)) +
  geom_violin(scale = "width", colour = NA, alpha = 0.85) + geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white", alpha = 0.6) +
  scale_fill_manual(values = org_cols, name = NULL) +
  labs(x = NULL, y = "CytoTRACE2 score", title = "Intratumour heterogeneity of CytoTRACE2 (malignant cells, per tumour)",
       caption = sprintf("October 2026 atlas: %s malignant cells, %d tumours with >= %d cells, %d patients. Ordered by median.",
                         format(nrow(m), big.mark = ","), uniqueN(m$orig.ident), MIN_CELLS, uniqueN(m$Patient_nr))) +
  theme_classic(base_size = 11) + theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 8), legend.position = "top",
                                        plot.caption = element_text(hjust = 0))
save2(p1, "F1_cytotrace2_per_tumour_violin", 13, 5.5)

# F2 forest, raw and depth-adjusted, one row per tumour
fl <- melt(per, id.vars = c("orig.ident", "origin", "n"), measure.vars = c("rho_cnv", "rho_cnv_adj"), variable.name = "what", value.name = "rho")
fl[, what := factor(c(rho_cnv = "raw", rho_cnv_adj = "adjusted for genes detected")[as.character(what)],
                    c("raw", "adjusted for genes detected"))]
fl[, tum := factor(orig.ident, per[order(rho_cnv), orig.ident])]
sub2 <- sprintf("Median rho %s raw, %s adjusted (%d tumours; sign test p = %s raw, %s adjusted). Between patients: primaries rho %s (p %s), metastases rho %s (p %s).",
                fmt(median(per$rho_cnv, na.rm = TRUE)), fmt(median(per$rho_cnv_adj, na.rm = TRUE)), sum(!is.na(per$rho_cnv)),
                fmt(sign_p(per$rho_cnv)), fmt(sign_p(per$rho_cnv_adj)),
                fmt(btw[origin == "primary", rho]), fmt(btw[origin == "primary", p]), fmt(btw[origin == "metastasis", rho]), fmt(btw[origin == "metastasis", p]))
p2 <- ggplot(fl[!is.na(rho)], aes(rho, tum)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey50") +
  geom_point(aes(colour = origin, shape = what, size = n), alpha = 0.85) +
  scale_colour_manual(values = org_cols, name = NULL) + scale_shape_manual(values = c(16, 1), name = NULL) +
  scale_size_area(max_size = 4, name = "malignant cells", breaks = c(500, 2000, 8000)) +
  labs(x = "Within-tumour Spearman rho, aberrant chromosomes vs CytoTRACE2", y = NULL,
       title = "Genomic complexity (aberrant chromosomes) vs stemness, per tumour",
       caption = paste(strwrap(sub2, 115), collapse = "\n")) +
  theme_classic(base_size = 11) + theme(axis.text.y = element_text(size = 7), plot.caption = element_text(hjust = 0))
save2(p2, "F2_genomic_complexity_vs_stemness_per_tumour", 9, 9)

# F3 genes detected vs CytoTRACE2 per tumour
per[, tum := factor(orig.ident, per[order(rho_genes), orig.ident])]
p3 <- ggplot(per, aes(rho_genes, tum, colour = origin, size = n)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey50") + geom_point(alpha = 0.85) +
  scale_colour_manual(values = org_cols, name = NULL) + scale_size_area(max_size = 4, name = "malignant cells", breaks = c(500, 2000, 8000)) +
  coord_cartesian(xlim = c(min(0, per$rho_genes), 1)) +
  labs(x = "Within-tumour Spearman rho, genes detected per cell vs CytoTRACE2", y = NULL,
       title = "Transcriptional complexity (genes detected) vs stemness, per tumour",
       caption = sprintf("Median rho %s over %d tumours (sign test p = %s). CytoTRACE1 used genes detected as its stemness proxy.",
                         fmt(median(per$rho_genes)), nrow(per), fmt(sign_p(per$rho_genes)))) +
  theme_classic(base_size = 11) + theme(axis.text.y = element_text(size = 7), plot.caption = element_text(hjust = 0))
save2(p3, "F3_genes_detected_vs_stemness_per_tumour", 9, 9)
cat("wrote", OUT, "\n")
