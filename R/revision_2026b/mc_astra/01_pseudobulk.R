# UM Atlas -> Sample x majority_celltype pseudobulk counts for mc-ASTRA (MOFAcell).
# Rerun of mc_astra/01_pseudobulk.R (July object) on the October 2026 object. 2026-10-05.
# Output (pb/): counts.mtx (genes x groups), genes.tsv, obs.csv (one row per group), samples.csv.
# Changes vs July: qs2 input, Condition (sequencing batch) + orig.ident + secondary_mut_sequencing
# kept as sample covariates; short_name added so samples can be matched to the July run.
library(Seurat)
library(Matrix)

out <- "../mc_astra_2026b/pb"
dir.create(out, showWarnings = FALSE)

# nthreads = 1: multi-threaded qs2 reads of multi-GB Seurat objects can deadlock silently.
obj <- qs2::qs_read("../new_data_2026b/Objects/full Dataset/final_Subtyping_iCNV_reclustered.qs2",
                    nthreads = 1)
stopifnot("not the October 2026 object" = ncol(obj) == 349696 && length(unique(obj$Sample)) == 63)
DefaultAssay(obj) <- "RNA"
obj <- JoinLayers(obj)
md <- obj@meta.data
cat("cells", ncol(obj), "samples", length(unique(md$Sample)), "patients", length(unique(md$Patient_nr)), "\n")
print(table(md$majority_celltype, useNA = "ifany"))

# GNA1Q was fixed at source in October; kept as a no-op guard.
md$primary_mutation <- ifelse(md$primary_mutation %in% c("GANQ", "GNA1Q"), "GNAQ", md$primary_mutation)

counts <- LayerData(obj, assay = "RNA", layer = "counts")
grp <- factor(paste(md$Sample, md$majority_celltype, sep = "||"))
pb <- counts %*% t(fac2sparse(grp))           # genes x groups, summed raw counts
rm(obj, counts); gc()

key <- do.call(rbind, strsplit(levels(grp), "||", fixed = TRUE))
obs <- data.frame(group = levels(grp), Sample = key[, 1], celltype = key[, 2],
                  psbulk_cells = as.integer(table(grp)[levels(grp)]),
                  psbulk_counts = Matrix::colSums(pb))

covars <- c("Sample", "orig.ident", "Condition", "Patient_nr", "origin", "organ", "primary_mutation",
            "secondary_mutation", "secondary_mut_sequencing", "primary_location", "healthy_location",
            "metastatic_location", "Gender", "Age_atDiagnosisY")
stopifnot(all(covars %in% colnames(md)))
samples <- unique(md[, covars])
stopifnot(!anyDuplicated(samples$Sample))      # one metadata row per Sample
samples$n_cells <- as.integer(table(md$Sample)[samples$Sample])
# Library prefix and _GEX suffix stripped: unique for all 63, and equal to the July name for every library.
samples$short_name <- sub("_GEX.*$", "", sub("^\\d+_\\d+-", "", samples$Sample))
stopifnot(!anyDuplicated(samples$short_name))

writeMM(pb, file.path(out, "counts.mtx"))
writeLines(rownames(pb), file.path(out, "genes.tsv"))
write.csv(obs, file.path(out, "obs.csv"), row.names = FALSE)
write.csv(samples, file.path(out, "samples.csv"), row.names = FALSE)
cat("groups", ncol(pb), "genes", nrow(pb), "\n")
