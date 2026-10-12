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

Helpers: `config_colors_2026b.R` (palettes), `plaid_score.R` (gene-set scorer used for the GeneNMF
programs: log1p CP10k, per-gene z-score over the scored cells, mean z per set).

Software: R 4.6 with Seurat 5, GloScope, vegan, scECODA, sccomp 2.1 (cmdstan, `inference_method = "hmc"`),
GeneNMF, CytoTRACE2 1.1.0; mc-ASTRA 0.1.0 with MOFA-FLEX 0.2.0rc2 (Python 3.13). The Weizmann 3CA table is
read from `../reference/Weizman_3CA_meta_programs_2025-01-29.xlsx`.
