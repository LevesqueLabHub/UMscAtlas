# UM Atlas GeneNMF test (Paul, 2026-10-08): rerun multiNMF with ribosomal, mitochondrial and
# immediate-early / heat-shock genes blocked from the HVG pool, and test MP enrichment with goseq
# (gene length as the bias covariate) next to the plain clusterProfiler ORA, for the new and the
# original MPs. Same samples, seed and per-sample downsampling as 17_GeneNMF_Cancer_Gene_Modules.Rmd,
# so the blocklist is the only change. BLOCK_SET (env): "ribo_mito_stress" (default, also blocks
# immediate-early / heat-shock genes) or "ribo_mito" (Victoria and Aizhan, 2026-10-08: keep the stress
# genes, a known UM metaprogram). Outputs: Outputs/GeneNMF_2026b/blocklist_test[_ribo_mito]/.
suppressPackageStartupMessages({ library(Seurat); library(dplyr); library(clusterProfiler); library(org.Hs.eg.db) })

base <- ".."
nmf_dir <- file.path(base, "Outputs/GeneNMF_2026b")
block_set <- Sys.getenv("BLOCK_SET", "ribo_mito_stress")
stopifnot(block_set %in% c("ribo_mito_stress", "ribo_mito"))
out <- file.path(nmf_dir, if (block_set == "ribo_mito") "blocklist_test_ribo_mito" else "blocklist_test"); dir.create(out, showWarnings = FALSE)
atlas_file <- file.path(base, "new_data_2026b/Objects/full Dataset/final_Subtyping_iCNV_reclustered.qs2")
universe <- qs2::qs_read(file.path(nmf_dir, "atlas_inputs_for_report.qs2"), nthreads = 1)$universe

# ---- blocklist ----
ieg_hs <- c("FOS", "FOSB", "FOSL1", "FOSL2", "JUN", "JUNB", "JUND", "EGR1", "EGR2", "EGR3", "IER2", "IER3", "IER5",
            "ATF3", "ZFP36", "ZFP36L1", "ZFP36L2", "DUSP1", "DUSP6", "NR4A1", "NR4A2", "NR4A3", "KLF2", "KLF4", "KLF6",
            "PPP1R15A", "GADD45B", "GADD45G", "BTG1", "BTG2", "CCN1", "CCN2", "SOCS3", "NFKBIA", "MAFF", "HES1", "DDIT3",
            "DDIT4", "RHOB", "MCL1", "SGK1", "CEBPD", "CEBPB", "TXNIP", "BAG3", "CRYAB", "HSPH1", "UBC", "IRF1", "MYC")
blocked <- function(g) grepl("^RP[SL][0-9]|^RPLP|^MT-|^MRP[SL][0-9]", g) |
  (block_set == "ribo_mito_stress" & (grepl("^HSP[ABDEH]|^HSP90|^DNAJ", g) | g %in% ieg_hs))
block <- universe[blocked(universe)]
cat("blocklist:", length(block), "genes (ribosomal", sum(grepl("^RP[SL][0-9]|^RPLP|^MRP[SL]", block)),
    ", mitochondrial", sum(grepl("^MT-", block)), ", IEG / heat shock", sum(!grepl("^RP[SL][0-9]|^RPLP|^MRP[SL]|^MT-", block)), ")\n")
writeLines(block, file.path(out, "blocklist.txt"))

# ---- program pool (cached) ----
prog_file <- file.path(out, "geneNMF_programs_blocklist.qs2")
if (!file.exists(prog_file)) {
  um_atlas <- qs2::qs_read(atlas_file, nthreads = 1)
  stopifnot(ncol(um_atlas) == 349696)
  malignant <- subset(um_atlas, subset = majority_celltype == "Melanoma cells"); rm(um_atlas); invisible(gc())
  malig_counts <- table(malignant$orig.ident)
  samples_for_nmf <- names(malig_counts[malig_counts >= 300])
  set.seed(42)  # identical to the report
  sampleObj.list <- Seurat::SplitObject(subset(malignant, subset = orig.ident %in% samples_for_nmf), split.by = "orig.ident")
  sampleObj.list <- lapply(sampleObj.list, function(obj) if (ncol(obj) > 500) subset(obj, cells = sample(colnames(obj), 500)) else obj)
  rm(malignant); invisible(gc())
  programs <- GeneNMF::multiNMF(sampleObj.list, assay = "RNA", k = 4:9, min.exp = 0.05, hvg.blocklist = block)
  qs2::qs_save(programs, prog_file, nthreads = 1)
}
programs <- qs2::qs_read(prog_file, nthreads = 1)
orig_programs <- qs2::qs_read(file.path(nmf_dir, "geneNMF_metaprograms.qs2"), nthreads = 1)$programs
stopifnot("different sample/k set from the report's pool" = setequal(names(programs), names(orig_programs)))

# ---- nMP sweep: metrics as in 17b ----
sweep <- bind_rows(lapply(3:10, function(n) {
  m <- GeneNMF::getMetaPrograms(programs, nMP = n)$metaprograms.metrics
  data.frame(nMP = n, MP = rownames(m), lapply(m, function(x) as.numeric(as.character(x))))
}))
write.csv(sweep, file.path(out, "nMP_sweep_blocklist.csv"), row.names = FALSE)

# ---- gene sets and lengths ----
gw <- GenomicFeatures::genes(EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86)
gene_len <- tapply(BiocGenerics::width(gw), gw$gene_name, max) / 1000
msig <- function(...) msigdbr::msigdbr(species = "Homo sapiens", ...) %>% dplyr::select(gs_name, gene_symbol)
sets <- list(`GO BP` = msig(collection = "C5", subcollection = "GO:BP"), Hallmark = msig(collection = "H"),
             Reactome = msig(collection = "C2", subcollection = "CP:REACTOME"))
w3 <- readxl::read_excel("../reference/Weizman_3CA_meta_programs_2025-01-29.xlsx", sheet = "Malignant")
w3 <- lapply(as.list(w3), function(g) unique(na.omit(g)))

# ORA (clusterProfiler::enricher, as the report) and goseq (Wallenius, length bias) on one MP
enrich_mp <- function(genes, univ) {
  univ_len <- univ[!is.na(gene_len[univ])]
  pwf <- goseq::nullp(setNames(as.integer(univ_len %in% genes), univ_len), bias.data = unname(gene_len[univ_len]), plot.fit = FALSE)
  bind_rows(lapply(names(sets), function(src) {
    s <- sets[[src]]
    ora <- enricher(genes, TERM2GENE = s, universe = univ, pvalueCutoff = 1, qvalueCutoff = 1, minGSSize = 10, maxGSSize = 500)
    ora <- if (is.null(ora)) tibble() else as_tibble(ora@result) %>% transmute(term = ID, overlap = Count, p_ora = pvalue, padj_ora = p.adjust)
    s2 <- dplyr::filter(s, gene_symbol %in% univ_len) %>% group_by(gs_name) %>% dplyr::filter(between(n(), 10, 500)) %>% ungroup()
    gs <- suppressMessages(goseq::goseq(pwf, gene2cat = as.data.frame(s2[, c("gene_symbol", "gs_name")]), method = "Wallenius", use_genes_without_cat = TRUE))
    gs <- as_tibble(gs) %>% transmute(term = category, p_goseq = over_represented_pvalue) %>% mutate(padj_goseq = p.adjust(p_goseq, "BH"))
    full_join(ora, gs, by = "term") %>% mutate(source = src)
  }))
}

describe_cut <- function(mp, label, univ) {
  genes <- mp$metaprograms.genes
  m <- mp$metaprograms.metrics; m[] <- lapply(m, function(x) as.numeric(as.character(x)))
  summ <- bind_rows(lapply(names(genes), function(k) { g <- genes[[k]]
    j <- which.max(sapply(w3, function(x) length(intersect(g, x)) / length(x)))
    tibble(run = label, MP = k, n_genes = length(g), top_genes = paste(head(g, 12), collapse = ", "),
           median_kb = median(gene_len[g], na.rm = TRUE), share_over_100kb = mean(gene_len[g] > 100, na.rm = TRUE),
           share_blocked = mean(blocked(g)), sampleCoverage = m[k, "sampleCoverage"], silhouette = m[k, "silhouette"],
           best_3CA = names(w3)[j], best_3CA_shared = length(intersect(g, w3[[j]])), best_3CA_size = length(w3[[j]])) }))
  enr <- bind_rows(lapply(names(genes), function(k) enrich_mp(genes[[k]], univ) %>% mutate(run = label, MP = k)))
  list(summary = summ, enrichment = enr)
}

univ_new <- setdiff(universe, block)
res <- list()
for (n in c(3, 6)) {
  mp_file <- file.path(out, sprintf("geneNMF_metaprograms_blocklist_nMP%d.qs2", n))
  if (!file.exists(mp_file)) qs2::qs_save(GeneNMF::getMetaPrograms(programs, nMP = n), mp_file, nthreads = 1)
  res[[paste0("blocklist_nMP", n)]] <- describe_cut(qs2::qs_read(mp_file, nthreads = 1), sprintf("blocklist nMP%d", n), univ_new)
  res[[paste0("original_nMP", n)]] <- describe_cut(qs2::qs_read(file.path(nmf_dir, sprintf("geneNMF_metaprograms_nMP%d.qs2", n)), nthreads = 1),
                                                   sprintf("original nMP%d", n), universe)
}
summ <- bind_rows(lapply(res, `[[`, "summary")); enr <- bind_rows(lapply(res, `[[`, "enrichment"))
write.csv(summ, file.path(out, "mp_summary.csv"), row.names = FALSE)
write.csv(enr, file.path(out, "mp_enrichment_ora_vs_goseq.csv"), row.names = FALSE)
genes_xlsx <- unlist(lapply(c(3, 6), function(n) { g <- qs2::qs_read(file.path(out, sprintf("geneNMF_metaprograms_blocklist_nMP%d.qs2", n)), nthreads = 1)$metaprograms.genes
  setNames(lapply(g, function(x) data.frame(gene = x)), sprintf("nMP%d_%s", n, names(g))) }), recursive = FALSE)
writexl::write_xlsx(genes_xlsx, file.path(out, "GeneNMF_blocklist_GeneLists.xlsx"))

# ---- console digest ----
options(width = 220)
print(as.data.frame(summ %>% mutate(across(where(is.numeric), ~ round(.x, 2)))), right = FALSE)
cat("\nnMP sweep (blocklist): mean silhouette, MPs with negative silhouette, min coverage\n")
print(sweep %>% group_by(nMP) %>% summarise(mean_sil = round(mean(silhouette), 3), n_neg = sum(silhouette < 0), min_cov = round(min(sampleCoverage), 2)))
cat("\nTop terms per MP (ORA padj < 0.05 or goseq padj < 0.05), best 4 per MP and source by goseq\n")
print(as.data.frame(enr %>% dplyr::filter(padj_ora < 0.05 | padj_goseq < 0.05) %>% group_by(run, MP, source) %>%
  slice_min(padj_goseq, n = 4, with_ties = FALSE) %>% ungroup() %>%
  transmute(run, MP, source, term = substr(term, 1, 55), overlap, padj_ora = signif(padj_ora, 2), padj_goseq = signif(padj_goseq, 2))), right = FALSE)
cat("\ndone\n")
