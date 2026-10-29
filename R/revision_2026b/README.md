# Revision analyses on the October 2026 atlas

Code for the analyses re-run for the revision on the reprocessed atlas object
(`final_Subtyping_iCNV_reclustered.qs2`: 349,696 cells, 63 libraries, 45 patients; TME subset
`TMESubset_Revision.qs2`, 149,062 cells). The scripts one folder up are the original submission and are
left unchanged.

Paths are relative: `..` is the analysis directory that holds the atlas objects
(`../new_data_2026b/Objects/...`), the metadata dumps (`../new_data_2026b/*.csv.gz`) and the outputs.
Reports `source("config_colors_2026b.R")` from this folder, so render them from here. Outputs and rendered
reports are not in the repository.

| Analysis | Files | What it answers |
|---|---|---|
| GloScope + scECODA | `03b_GloScope_Full_Cohort.Rmd` | Sample-level divergence (GloScope KL, ANOSIM) for origin, secondary mutation, location, primary mutation, organ, sex, age; one-library-per-patient sensitivity; tumour vs own normal with within-patient permutations; scECODA CLR composition vs GloScope on the same samples. |
| scECODA benchmark | `25_ECODA_vs_GloScope_Benchmark.R` | ECODA vs GloScope grouping recovery (July 2026 object, `GloScope_Full_Outputs_2026`). |
| CytoTRACE2 | `run_cytotrace2_2026b.R` (scores, per library), `01b_CytoTRACE2_2026b.Rmd` (report) | Developmental potential by origin, location and mutation; stemness signature genes. |
| sccomp | `04_sccomp_Compositional_Analysis.Rmd` | TME composition by origin, secondary mutation and location, with `(1 \| Patient_nr)`. |
| sccomp, T/NK states | `04c_*` (paired model check, 14 T/NK subtypes), `04d_*` (subtypes combined into lineage and state groups, 7 tumour-normal pairs), `04e_*` (state groups, whole cohort), `04f_*` (discovery in the 7 pairs, validation in samples of other patients) | Which T/NK states differ between normal tissue and the primary tumour of the same patient, and whether that holds in unpaired patients. Each `.R` fits, each `.Rmd` reports. |
| GeneNMF | `17_GeneNMF_Cancer_Gene_Modules.Rmd` (+ `17_nmp_section.Rmd`, `17_extras_*.Rmd` as knitr children), `17b`-`17f` | De novo meta-programs in malignant cells (nMP sweep, follow-up scoring, mutation sweep, ribosomal/mitochondrial blocklist rerun with goseq); `17f_GeneNMF_Blocklist_2026b.Rmd` reports the blocklist rerun. |
| mc-ASTRA (MOFAcell) | `mc_astra/01_pseudobulk.R`, `02_mofacell.py`, `04_stability_summary.py`, `05_network.py`; report `26_mcASTRA_MOFAcell_2026b.Rmd` | Multicellular factors across samples and their association with origin, mutation, location and batch. The stability refits run `02_mofacell.py` with `SEED` 1-5 x `N_FACTORS` 4/6/8 and `OUT_DIR=../mc_astra_2026b/stability/K<k>_seed<s>`. |
| 3CA metaprograms | `score_3CA_ucell_2026b.R` (UCell, used by the figures), `02b_score_3CA_2026b.R` + `02b_Cancer_Metaprograms_2026b.Rmd` (report; PLAID kept as a cross-check) | The 67 Weizmann 3CA malignant metaprograms in melanoma cells and melanocytes, by origin, location and mutation, patient-level. |
| Manuscript panels | `manuscript_panels_2026b_cytotrace_mp.R`, `manuscript_panels_2026b_gloscope.R`, `Fig1h_heatmap_from_csv.R`; report `Manuscript_Panels_2026b.Rmd` | The CytoTRACE2, metaprogram and GloScope figure panels (table below), one PDF + 300 dpi PNG per panel, named by figure ID, with the patient-level tests in `../Outputs/Manuscript_panels_2026b/`. |
| Xenium stemness vs CNV | `14b_score_xenium_boston_2026b.R` + `14b_Xenium_Stemness_Boston_2026b.Rmd` (6 Boston metastases), `14c_extract_pum26_icnv_obs.py`, `14e_pum26_icnv_unsmoothed.py`, `14c_score_xenium_pum26_2026b.R` + `14c_Xenium_Stemness_pUM26_2026b.Rmd` (pUM26 eye and liver), `14d_cnv_rms_per_cell.py` (per-cell CNV RMS), shared helpers `14b_stemness_common.R` | Spatial maps of cell type, CNV and stemness, and the per-cell stemness vs CNV correlation. |
| Xenium iCNV | `xenium_icnv_supp_figure.py`, `xenium_icnv_figures_only.py` | infercnvpy on all cells of each Boston section (macrophages as reference): chromosome profile per lineage, chr3 loss and chr8 gain maps of the melanoma cells. |
| Per-tumour CytoTRACE2 vs CNV | `23b_pertumor_cytotrace_cnv_2026b.R` | Within each tumour, does the number of aberrant chromosomes go with CytoTRACE2? |

Helpers: `config_colors_2026b.R` (palettes), `plaid_score.R` (gene-set scorer used for the GeneNMF
programs and the Xenium stemness scores: log1p CP10k, per-gene z-score over the scored cells, mean z per set).

Xenium inputs: `../Xenium_Boston_Sullivan/` (merged Boston object and the InSituCNV outputs
`06_InSituCNV_Outputs/`, `07_InSituCNV_Breast_Outputs/`) and `../Xenium_pUM26/SPLIT_Reannotation/` (pUM26
SPLIT-purified regions and iCNV h5ads).

## Figure panels and the code that makes them

Panel letters follow the revised manuscript. Panels not listed come from the scripts one folder up, or from
the reports in the table above (the October multicellular factors of Fig. 1l-m and Suppl. Fig. 3e-g:
`26_mcASTRA_MOFAcell_2026b.Rmd`).

| Panel | File name (`Figures_2026b/Manuscript_panels/`) | Script |
|---|---|---|
| Fig. 1h | `Fig1h_gloscope_anosim_summary` (table: `Fig1h_source/anosim_manuscript_panels.csv`) | `manuscript_panels_2026b_gloscope.R`; restyle from the table with `Fig1h_heatmap_from_csv.R` |
| Fig. 1i, 1j, 1k | `Fig1i_mds_origin`, `Fig1j_mds_pUM_secondary_mutation`, `Fig1k_mds_pUM_location` | `manuscript_panels_2026b_gloscope.R` |
| Suppl. Fig. 3b, 3c, 3d | `SFig3b_kl_heatmap`, `SFig3c_mds_samples_labelled`, `Fig1i_mds_primary_mutation` | `manuscript_panels_2026b_gloscope.R` |
| Fig. 2c, 2d, 2e, 2i | `Fig2c_cytotrace_potency_umap`, `Fig2d_cytotrace_violin_celltype`, `Fig2e_cytotrace_by_origin_location`, `Fig2i_cytotrace_pUM_location` | `manuscript_panels_2026b_cytotrace_mp.R` |
| Fig. 2g, 2h | `Fig2g_spatial_<section>`, `Fig2h_stemness_vs_cnv_boston` | `14b_Xenium_Stemness_Boston_2026b.Rmd` |
| Suppl. Fig. 4a, 4b, 4c | `SFig4a_cytotrace_origin_age`, `SFig4b_cytotrace_origin_gender`, `SFig4c_cytotrace_secondary_mutation_gender`; tested within each stage: `SFig4a_cytotrace_age_by_stage`, `SFig4b_cytotrace_gender_by_stage` | `manuscript_panels_2026b_cytotrace_mp.R` |
| Suppl. Fig. 5a, 5b | `SFig5a_spatial_pUM26_eye`, `SFig5a_spatial_pUM26_liver`, `SFig5b_stemness_vs_cnv_pUM26` | `14c_Xenium_Stemness_pUM26_2026b.Rmd` |
| Fig. 3a, 3c-j | `Fig3a_mp_top30_correlation`, `Fig3c/3e/3g_umap_MP45/MP17/MP9`, `Fig3d/3f/3h_box_MP45/MP17/MP9`, `Fig3i_*_MP22_24_origin`, `Fig3j_*_MP22_24_treatment` | `manuscript_panels_2026b_cytotrace_mp.R` |
| Suppl. Fig. 6a, 6b | `SFig6a_mp_top20_mean_ucell`, `SFig6b_box_MP45_pUM_location` | `manuscript_panels_2026b_cytotrace_mp.R` |
| Reviewer, fibroblasts by site | `Rev_L169_fibroblast_cytotrace_site` | `manuscript_panels_2026b_cytotrace_mp.R` |

Boxplot style, shared by every boxplot panel (`box_panel()` and `thm` in
`manuscript_panels_2026b_cytotrace_mp.R`): `theme_classic(base_size = 12)`; `geom_boxplot(width = 0.6,
alpha = 0.45, outlier.shape = NA)`; one point per tumour or tissue (`orig.ident` with at least 10 cells) drawn
with `ggbeeswarm::geom_quasirandom(width = 0.22, alpha = 0.85, size = 2.2)`; brackets with `ggsignif`. Tests are
two-sided Wilcoxon on patient means (`Patient_nr` is the unit), unadjusted; a group with fewer than 3 patients is
drawn but not tested.

Software: R 4.6 with Seurat 5, GloScope, vegan, scECODA, sccomp 2.1 (cmdstan, `inference_method = "hmc"`),
GeneNMF, CytoTRACE2 1.1.0, UCell, ggsignif, ggbeeswarm, scattermore; mc-ASTRA 0.1.0 with MOFA-FLEX 0.2.0rc2
(Python 3.13); infercnvpy 0.6.0 for the Xenium iCNV. The Weizmann 3CA table is read from
`../reference/Weizman_3CA_meta_programs_2025-01-29.xlsx`.
