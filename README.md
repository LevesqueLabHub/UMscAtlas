# UMscAtlas
Analysis of healthy, primary, and metastatic uveal melanoma using scRNA-seq and Xenium spatial transcriptomics. This repository includes workflows for preprocessing, integration, cell-type annotation, and characterization of malignant and TME states.

## Metadata variables documentation

| Metadata variable | Name in figure | Explanation |
| --- | --- | --- |
| primary_location  | Intraocular location  | Location within the eye where the primary tumor arose. NA for samples that are not primary UM |
| origin | Disease progression | Indicates if cells come from a healthy sample, primary tumor or metastasis |
| orig.ident | Sample | Identifies physical tissue sample cells originate from |
| location | Location | Location (intraocular location and location of metastasis) for all samples in the data set |
| samplename | Sample Name | Identifies biological replicates sequenced in different libraries (use instead of Sample column (that shows projectnumber and Index from FGCZ |
| healthy_location | Intraocular location Healthy | Identifies intraocular location from healthy tissue |
| Condition | Batch | Indicated Sequencing Batch for each Sample |
| primary_mutation | Primary Mutation | Idenfies  primary mutation (GNAQ, GNA11, PLCB4 or CYSLTR2) for each sample (assesed by WES or targeted sequencing) |
| secondary_mutation | Secondary Mutation | Idenfies  primary mutation (BAP1, SF3B1 or EIF1AX) for each sample (assesed by WES or targeted sequencing) |
| Patient_nr | Patient Number | assigns a Patient number (P1-P36) to each sequenced sample and shows which Patient has multiple samples sequenced |
| organ | Organ | Identifies organ from which sample is originating |
| Gender | Gender | Identifies gender of patient sample was derived from |
| TreatmentStage_ofProcessedSample | Treatment | Indicates last treatment of Patient before sampling of the sequenced Sample |
| Stage_atDiagnosis | Stage at Diagnosis | Disease stage at first Diagnosis |
| Mets_atDiagnosis | Metastasis at Diagnosis | Identifies if metastasis was present at disease prognosis |
| Age_atDiagnosisY | Age at Diagnosis | Identifies age at diagnosis of tumor in Years |
| RiskCategory | Risk Category | Identifies Risk Category assigned by collaborators in Australia for primary tumors |
| Time_sinceDiagnosis | Time since Diagnosis | Shows calculated time since Diagnosis in Years until October 2025 |
| PrimaryTreatment | Primary Treatment | Identifies treatment of primary tumor |
| LocalProgressionPrimary | Local Progression Primary Tumor | Indicates wheter primary tumor showed disease progression or not |
| LocalProgressionTreatment | Treatment of Local Progression | Identifies which treatment was used if didease was localy progressive |
| DevelopmentOfMetastasis | Metastasis Development | Indicates wheter diseade develops into metastatic disease from primary tumor |
| Time_Primary-Met.Y | Indicates time in years from primary tumor to metastasis development |
| MetastaticSite | Metastatic Site | Identifies site of Metastasis development |
| MetastasisTreatment | Treatment Metastasis | Treatment of metastastatic tumor |
| Time_sinceMetastasisDiagnosisY | Time since diagnosis of Metastasis | Identifies time since diagnosis of metastasis in years until october 2025 |
| PatientStatus | Status Patient | Indicates wheter patient is deceased or alive |
| LocationPrimaryTumor | Location of Primary Tumor | Identifies intaocular location of primary tumor also for metastases samples |
