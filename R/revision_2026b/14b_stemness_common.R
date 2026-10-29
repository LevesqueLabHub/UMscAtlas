# Shared by 14b (Boston) and 14c (pUM26): gene sets for the October 2026 stemness scoring and the
# figure helpers. Sourced by 14b_score_xenium_boston_2026b.R, 14c_score_xenium_pum26_2026b.R and both Rmds.
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(patchwork) })

# ---- Gene sets -------------------------------------------------------------------------------------------
stemness_sig_dir <- "../Outputs/CytoTRACE2_2026b"
# January signature of report 14 (top/bottom 50 malignant-correlated PANEL genes); copy kept by 14b.
jan_sig_dir <- "../Outputs/Xenium_Stemness_Boston_2026b/inputs"

# Cell cycle: Seurat cc.genes (Tirosh 2015) + the 2019 update, the list family report 14 used.
cycling_genes <- unique(c(Seurat::cc.genes.updated.2019$s.genes, Seurat::cc.genes.updated.2019$g2m.genes,
                          Seurat::cc.genes$s.genes, Seurat::cc.genes$g2m.genes))
# Heat-shock chaperones and immediate-early / dissociation-stress genes (after van den Brink et al. 2017).
heat_shock_genes <- c("HSPA1A", "HSPA1B", "HSPA2", "HSPA4", "HSPA4L", "HSPA5", "HSPA6", "HSPA8", "HSPA9", "HSPB1",
                      "HSPB8", "HSPD1", "HSPE1", "HSPH1", "HSP90AA1", "HSP90AB1", "HSP90B1", "DNAJA1", "DNAJA4",
                      "DNAJB1", "DNAJB4", "DNAJB6", "BAG3", "CRYAB", "CHORDC1", "STIP1", "SERPINH1", "HSF1")
ieg_genes <- c("FOS", "FOSB", "FOSL1", "FOSL2", "JUN", "JUNB", "JUND", "EGR1", "EGR2", "EGR3", "ATF3", "IER2", "IER3",
               "IER5", "NR4A1", "NR4A2", "NR4A3", "DUSP1", "KLF2", "KLF4", "KLF6", "ZFP36", "ZFP36L1", "PPP1R15A",
               "CCN1", "CCN2", "BTG2", "GADD45B", "GADD45G", "SOCS3", "MYC", "CEBPB", "CEBPD", "NFKBIA", "SGK1", "UBC")

# October panel-matched malignant signature, built the way the January one was: the 01b malignant gene-CytoTRACE2
# correlations (gene_cor_malignant.csv, already restricted to genes detected in >= 1% of malignant cells),
# restricted to the Xenium panel, top 50 by correlation = Top, bottom 50 = Bottom. Report 14's January files
# (cytotrace2_{top,bottom}_genes_xenium.csv) are reproduced exactly - same genes, order and values - by this rule
# applied to the January correlations (Outputs/cytotrace2_malignant_gene_correlations.csv); no other filter.
panel_sig_file <- file.path(stemness_sig_dir, "stemness_signature_malignant_xenium_panel_top50.csv")
build_panel_signature <- function(panel, n = 50) {
  cr <- fread(file.path(stemness_sig_dir, "gene_cor_malignant.csv"))[gene %in% panel & !is.na(correlation)]
  stopifnot(nrow(cr) >= 2 * n)
  rbind(cr[order(-correlation)][1:n, .(gene, correlation, direction = "Top")],
        cr[order(correlation)][1:n, .(gene, correlation, direction = "Bottom")])
}

stemness_sets <- function(panel) {
  new_sig <- rbindlist(lapply(c("malignant", "healthy_melanocyte", "all_cells"), function(s)
    fread(file.path(stemness_sig_dir, sprintf("stemness_signature_%s_top50.csv", s)))[, signature := s]))
  stopifnot(new_sig[, .N, by = .(signature, direction)]$N == 50)
  panel_sig <- build_panel_signature(panel)
  if (file.exists(panel_sig_file))   # the saved file must be what this panel gives (Boston and pUM26 share it)
    stopifnot(all.equal(as.data.frame(fread(panel_sig_file)), as.data.frame(panel_sig)))
  jan_top <- fread(file.path(jan_sig_dir, "cytotrace2_top_genes_xenium.csv"))$gene
  jan_bot <- fread(file.path(jan_sig_dir, "cytotrace2_bottom_genes_xenium.csv"))$gene
  hs_ieg <- unique(c(heat_shock_genes, ieg_genes))
  on_panel <- function(g) intersect(g, panel)
  sets <- list()
  for (s in unique(new_sig$signature)) for (d in c("Top", "Bottom"))
    sets[[paste0("new_", s, "_", d)]] <- on_panel(new_sig[signature == s & direction == d, gene])
  sets$panel_malignant_Top    <- panel_sig[direction == "Top", gene]      # PRIMARY score
  sets$panel_malignant_Bottom <- panel_sig[direction == "Bottom", gene]
  sets$jan_malignant_Top    <- on_panel(jan_top)
  sets$jan_malignant_Bottom <- on_panel(jan_bot)
  sets$Cycling       <- on_panel(cycling_genes)
  sets$HeatShock_IEG <- on_panel(hs_ieg)
  # Controls: the malignant Top set without the genes it shares with each reference set, so a Top-vs-heat-shock
  # (or cycling) correlation is not just shared genes correlating with themselves.
  sets$panel_malignant_Top_noHS <- setdiff(sets$panel_malignant_Top, hs_ieg)
  sets$panel_malignant_Top_noCC <- setdiff(sets$panel_malignant_Top, cycling_genes)
  sets$new_malignant_Top_noHS <- setdiff(sets$new_malignant_Top, hs_ieg)
  stopifnot(lengths(sets) > 0)
  cov_row <- function(sig, dir, g) data.table(signature = sig, direction = dir, n_signature = length(g),
    n_on_panel = sum(g %in% panel), genes_on_panel = paste(on_panel(g), collapse = ", "),
    genes_missing = paste(setdiff(g, panel), collapse = ", "))
  coverage <- rbind(new_sig[, cov_row(signature, direction, gene), by = .(signature, direction)][, -(1:2)],
                    cov_row("january_malignant", "Top", jan_top), cov_row("january_malignant", "Bottom", jan_bot))
  list(sets = sets, coverage = coverage, new_sig = new_sig, panel_sig = panel_sig, jan_top = jan_top, jan_bot = jan_bot,
       cycling = cycling_genes, heat_shock = heat_shock_genes, ieg = ieg_genes, panel = panel)
}

# ---- Figures ---------------------------------------------------------------------------------------------
cytotrace_colors <- c("#000004FF", "#3B0F70FF", "#8C2981FF", "#DE4968FF", "#FE9F6DFF", "#FCFDBFFF")
rho <- function(a, b) suppressWarnings(cor(a, b, method = "spearman", use = "complete.obs"))

# Three side-by-side maps of one section: cell types | stemness | CNV (layout of report 19).
# f: cells of one section with x, y, celltype, stem (score) and cnv columns; cnv_scale: a ggplot colour scale.
# Raster size and point radius for one section's map: the short side gets `short_px` pixels and the long side
# follows the tissue aspect ratio (so cells are not stretched under coord_fixed), capped at `max_px`. The point
# radius is `frac` of the median nearest-neighbour distance between cell centroids, in raster pixels: 0.5 leaves
# white gaps that dominate at display size and 0.75 still looked pale; 0.9 reads as filled tissue, a zoom still resolves cells.
spatial_raster <- function(x, y, short_px = 4000, max_px = 20000, frac = 0.9) {
  xr <- diff(range(x)) * 1.1; yr <- diff(range(y)) * 1.1     # ggplot's default 5% expansion on each side
  ratio <- max(xr, yr) / min(xr, yr)
  short_px <- min(short_px, max_px / ratio)                   # a cap shrinks BOTH sides, so the aspect is kept
  px <- if (yr >= xr) c(short_px, short_px * ratio) else c(short_px * ratio, short_px)
  nn <- median(RANN::nn2(cbind(x, y), k = 2)$nn.dists[, 2])
  list(px = round(px), nn_um = nn, pointsize = max(1, frac * nn * px[1] / xr))
}

spatial_3panel <- function(f, ct_cols, stem_lim, stem_title, cnv_scale, cnv_title, title, Q = c(0.01, 0.99),
                           ras = spatial_raster(f$x, f$y)) {
  # Manuscript panels (Victoria, 2026-10-15): legend text and titles 3x the theme default (8.8 -> 26.4 pt,
  # 11 -> 33 pt), x / y coordinate axes kept on every map.
  th <- theme_classic(base_size = 22) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5), legend.title = element_text(size = 33),
          legend.text = element_text(size = 26.4), legend.key.size = unit(0.55, "in"))
  ax <- labs(x = "x (µm)", y = "y (µm)")
  # scattermore rasterises itself; ggrastr::rasterise() inside a res-150 png() corrupted the heap (munmap_chunk abort)
  pt <- function() scattermore::geom_scattermore(pointsize = ras$pointsize, pixels = ras$px)
  p1 <- ggplot(f, aes(x, y, colour = celltype)) + pt() +
    scale_colour_manual(values = ct_cols, drop = TRUE, name = "Cell type") +
    guides(colour = guide_legend(override.aes = list(size = 9))) +
    coord_fixed() + th + ax + labs(title = "Cell types")
  p2 <- ggplot(f[order(stem)], aes(x, y, colour = stem)) + pt() +
    scale_colour_gradientn(colours = cytotrace_colors, limits = stem_lim, oob = scales::squish,
                           name = sprintf("Stemness\n(PLAID z)\ncapped at\n%g-%gth pct", 100 * Q[1], 100 * Q[2])) +
    coord_fixed() + th + ax + labs(title = stem_title)
  p3 <- ggplot(f[order(!is.na(cnv), abs(cnv))], aes(x, y, colour = cnv)) + pt() + cnv_scale +
    coord_fixed() + th + ax + labs(title = cnv_title)
  (p1 | p2 | p3) + plot_annotation(title = title, theme = theme(plot.title = element_text(size = 30, face = "bold", hjust = 0.5)))
}

# Figure height for a 34-inch-wide 3-panel map from the tissue aspect ratio (~6.2 in per panel next to the
# ~4 in wide 3x legends); the 9.5 in floor is what the cell-type legend (up to 12 levels at 26 pt) needs.
spatial_width <- 34
spatial_height <- function(x, y) round(min(26, max(9.5, 6.2 * diff(range(y)) / diff(range(x)) + 3)), 1)

# Stemness vs per-cell CNV RMS (14d_cnv_rms_per_cell.py) in melanoma cells, the design of the submitted Fig. 2h
# (fix_b7_correlation.R, B7_cnv_stemness_correlation_melanoma). m: sample, stem, cnv_rms.
rms_csv <- "../Outputs/Xenium_Stemness_Boston_2026b/cnv_rms_per_cell.csv"
# stat_cor's default label, except that a p between 2.2e-16 and 1e-4 gets 2 digits (ggpubr 1.0.0 printed 9.959955e-14)
cor_lab <- aes(label = after_stat(sprintf("italic(R)~`=`~%s*`,`~italic(p)~%s", signif(r, 2),
                 ifelse(p < 2.2e-16, "`<`~2.2e-16", paste0("`=`~", signif(p, 2))))))
stem_cnv_scatter <- function(m, xlab, ncol = 3, labeller = "label_value", ylab = "CNV score (per-cell RMS)") {
  # Grey points as the submitted panel (Victoria, 2026-10-16; the density colouring of 2026-10-15 is reverted).
  ggplot(as.data.frame(m), aes(stem, cnv_rms)) +
    geom_point(colour = "grey50", alpha = 0.2, size = 0.5) +
    geom_smooth(method = "lm", formula = y ~ x, colour = "black", se = TRUE, linewidth = 0.8) +
    # white box: bare text over the dense points was unreadable (2026-10-16)
    ggpubr::stat_cor(cor_lab, method = "pearson", label.x.npc = 0.05, label.y.npc = 0.95, size = 3.5,
                     geom = "label", fill = "white", border.colour = NA) +
    facet_wrap(~sample, scales = "free", ncol = ncol, labeller = labeller) +
    theme_classic(base_size = 12) + theme(strip.text = element_text(face = "bold")) +
    labs(x = xlab, y = ylab, title = "Melanoma cells: CNV score vs stemness score (per section)")
}

# Manuscript panels (PDF + PNG); returns the plot visibly, so a chunk ending in save_panel() also shows it
panel_dir <- "../Figures_2026b/Manuscript_panels"
save_panel <- function(p, name, w, h, dpi = 300) {
  dir.create(panel_dir, recursive = TRUE, showWarnings = FALSE)
  ggsave(file.path(panel_dir, paste0(name, ".png")), p, width = w, height = h, dpi = dpi, bg = "white")
  ggsave(file.path(panel_dir, paste0(name, ".pdf")), p, width = w, height = h, device = cairo_pdf)
  p
}

# Cell-level density of two variables, one facet per section, Spearman rho printed
hex_by_sample <- function(dat, xv, yv, xlab, ylab, title, who = "melanoma cells only", ncol = 3) {
  r <- dat[, .(rho = rho(get(xv), get(yv)), n = .N), by = sample]
  r[, lab := sprintf("rho = %.2f\nn = %s cells", rho, format(n, big.mark = ","))]
  ggplot(dat, aes(.data[[xv]], .data[[yv]])) +
    geom_hex(bins = 50) + scale_fill_viridis_c(trans = "log10", name = "Cells") +
    geom_smooth(method = "lm", se = FALSE, colour = "red", linewidth = 0.6) +
    geom_text(data = r, aes(x = -Inf, y = Inf, label = lab), hjust = -0.05, vjust = 1.2, size = 3.2, inherit.aes = FALSE) +
    facet_wrap(~sample, scales = "free", ncol = ncol) + theme_bw(base_size = 11) +
    labs(x = xlab, y = ylab, title = title, subtitle = paste("Cell-level Spearman rho within each section,", who))
}

# Per-section rho lollipop. ps: sample, score, covariate, rho; pairs: score, covariate, label
lolli <- function(ps, pairs, title, sample_levels) {
  x <- merge(ps, pairs, by = c("score", "covariate"))
  x[, `:=`(sample = factor(sample, levels = rev(sample_levels)), label = factor(label, levels = pairs$label))]
  ggplot(x, aes(rho, sample)) + geom_vline(xintercept = 0, colour = "grey50") +
    geom_segment(aes(x = 0, xend = rho, yend = sample), colour = "grey60") + geom_point(size = 3, colour = "#B00068") +
    geom_text(aes(label = sprintf("%.2f", rho)), nudge_y = 0.3, size = 3.2) +
    facet_wrap(~label, ncol = 3) + coord_cartesian(xlim = c(min(-0.1, x$rho) - 0.05, max(x$rho) + 0.08)) +
    theme_bw(base_size = 11) + theme(plot.margin = margin(5.5, 18, 5.5, 5.5)) +
    labs(x = "Spearman rho (cell-level, within section)", y = NULL, title = title)
}
