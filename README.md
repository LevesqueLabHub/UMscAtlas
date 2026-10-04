# UMscAtlas
Analysis of healthy, primary, and metastatic uveal melanoma using scRNA-seq and Xenium spatial transcriptomics. This repository includes workflows for preprocessing, integration, cell-type annotation, and characterization of malignant and TME states.


## Metadata variables documentation

| Metadata variable | Name in figure | Explanation |
| --- | --- | --- |
| primary_location  | Intraocular Location  | Location within the eye where the primary tumor arose. NA for samples that are not primary UM |
| origin | Disease Progression | Indicates if cells come from a healthy sample, primary tumor or metastasis |
| orig.ident | Sample | Identifies physical tissue sample cells originate from |
| location | Location | Location (intraocular location and location of metastasis) for all samples in the data set |
| samplename | Sample Name | Identifies biological replicates sequenced in different libraries (use instead of Sample column (that shows projectnumber and Index from FGCZ |
| healthy_location | Intraocular Location Healthy | Identifies intraocular location from healthy tissue |
| Condition | Batch | Indicates sequencing batch for each sample |
| primary_mutation | Primary Mutation | Identifies primary mutation (GNAQ, GNA11, PLCB4 or CYSLTR2) for each sample (assessed by WES or targeted sequencing) |
| secondary_mutation | Secondary Mutation | Identifies secondary mutation (BAP1, SF3B1 or EIF1AX) for each sample (assessed by WES or targeted sequencing) |
| Patient_nr | Patient Number | assigns a patient number (P1-P36) to each sequenced sample and shows which patients have multiple samples sequenced |
| organ | Organ | Identifies organ from which sample is originating |
| Gender | Gender | Identifies gender of patient sample was derived from |
| TreatmentStage_ofProcessedSample | Treatment | Indicates last treatment of patient before sampling of the sequenced sample |
| Stage_atDiagnosis | Stage at Diagnosis | Disease stage at first diagnosis |
| Mets_atDiagnosis | Metastasis at Diagnosis | Identifies if metastasis was present at disease prognosis |
| Age_atDiagnosisY | Age at Diagnosis | Identifies age at diagnosis of tumor in years |
| RiskCategory | Risk Category | Identifies risk category assigned by collaborators in Australia for primary tumors |
| Time_sinceDiagnosis | Time since Diagnosis | Shows calculated time since diagnosis in years until October 2025 |
| PrimaryTreatment | Primary Treatment | Identifies treatment of primary tumor |
| LocalProgressionPrimary | Local Progression Primary Tumor | Indicates wheter primary tumor showed disease progression or not |
| LocalProgressionTreatment | Treatment of Local Progression | Identifies which treatment was used if didease was locally progressive |
| DevelopmentOfMetastasis | Metastasis Development | Indicates whether disease developed into metastatic disease from primary tumor |
| Time_Primary-Met.Y | Time to Metastasis | Indicates time in years from primary tumor to metastasis development |
| MetastaticSite | Metastatic Site | Identifies site of metastasis development |
| MetastasisTreatment | Treatment Metastasis | Treatment of metastatic tumor |
| Time_sinceMetastasisDiagnosisY | Time since Diagnosis of Metastasis | Identifies time since diagnosis of metastasis in years until october 2025 |
| PatientStatus | Status Patient | Indicates whether patient is deceased or alive |
| LocationPrimaryTumor | Location of Primary Tumor | Identifies intaocular location of primary tumor also for metastases samples |

## Symbols to Use in Plots for h-p-m
- **healthy**: ● Circle (nr. 1)
- **primary**: ▲ Triangle (nr. 2)
- **metastases**: ◆ Diamond (nr. 5)

## Symbols to Use in Plots for primary mutation
- **GNAQ**: ● Circle (nr. 1)
- **GNA11**: ▲ Triangle (nr. 2)
- **PLCB4**: ◆ Diamond (nr. 5)
- **none**: ▼ Triangle, point down (nr. 6)
- **unknown**: ✵ Star (nr. 8)

## Available Datasets

- **Full Dataset**: All cells across all samples.
- **TME Subset**: Tumor microenvironment cells, subsetted from Full Dataset (all cells except Melanocytes and Melanoma cells).
- **Melanoma/Melanocytes Subset**: Melanoma cells and melanocytes from healthy - primary - metastases subsetted from Full Dataset.
- **Melanoma Subset**: Subset of melanoma cells from p-m UM samples subsetted from Full Dataset.

> **Notes:**
> - Melanoma cells were defined using inferred copy number variation (iCNV) profiles to distinguish malignant cells from non-malignant Melanocytes populations.
> - All subsets are **Harmony-corrected**, except the **Full Dataset**.

```mermaid
graph TD
  A[<b>Full Dataset</b><br/>13_final_Dataset_reclustered_28012025.qs<br/>*not integrated*] --> B[<b>TME Subset</b><br/>01_integrated_multicorr_TME.qs<br/>*harmony corrected*]
  A --> C[<b>Melanoma/Melanocytes Subset</b><br/>01_integrated_melanoma_melanocyte_subset.qs<br/>*harmony corrected*]
  B --> E[<b>T/NK cell Subset</b><br/>01_integrated_TNK_subset.qs<br/>*harmony corrected*]
  C --> D[<b>Melanoma Subset</b><br/>01_integrated_melanoma_subset.qs<br/>*harmony corrected*]
```

## Reference genome and annotation

**Pinned in [`R/genome_annotation.R`](R/genome_annotation.R).** `source()` it in any script that
filters, renames or joins genes.

```
Homo_sapiens/GENCODE/GRCh38.p13/Annotation/Release_42-2023-01-30
```

GENCODE 42 / Ensembl 108, 62,684 genes. All libraries in the atlas were aligned and quantified
against this one reference, so genome build and annotation are not a batch or cohort difference.

Two things the pin exists to catch:

- **The chromatin modalities are on a different build.** The scATAC and multiome (ARC) runs used
  `GRCh38.p14 / Release_48-2025-07-03`. Release 42 -> 48 changes the gene set and some symbols, so any
  ATAC or multiome join to the atlas must be checked for alias drift rather than assumed to line up.
- **Annotation packages are not part of the reference but do move results.** `org.Hs.eg.db` **3.21.0** was
  used for every GO/KEGG mapping (`Augur.R`, `Milo.R`, `CytoTRACE2.Rmd`,
  `DEA_degPatterns_GO_KEGG_of_location.R`). `check_annotation_packages()` warns when the installed
  version differs. The `GWASTools` version (`centromeres.hg38` in `CNV_Arm_Analysis.Rmd`) is not asserted.

`check_features_against_reference(object, features_file)` verifies a Seurat object's features against
the reference's gene annotation table and stops if fewer than 95 % are found, which is what a build
mismatch looks like.

## Metadata Structure in Objects

Each Seurat object contains the following components:

#### 1. meta.data
Cell-level annotations and sample information:
- **`majority_celltype`**: Broad cell type annotation.
- **`SubTyping`**: Detailed immune cell subtypes.
- **`orig.ident`**: Original sample identifier (use for all analyses except GloScope).
- **`samplename`**: Library identifier.

> **Tip:** Refer to the summary table at the top and the provided R script for recommended naming conventions and color palettes.

---

#### 2. reductions
Dimensionality reduction results for visualization and integration:
- **`umap`**: UMAP coordinates for clustering and visualization.
- **`pca`**: PCA embeddings for exploratory analysis.
- **Integration**: All subsets are Harmony-corrected (except the Full Dataset).
-   use reduction: **`umap`** for **Full Dataset**
-   use reduction: **`umap_Condition_and_orig.ident_50PC_theta2`** for **Harmony-corrected Datasets**

---

#### 3. assays
Gene expression data:
- **`RNA`**: Normalized gene expression matrix (default assay for most analyses).

