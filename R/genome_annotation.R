# Pinned reference genome and annotation for the UM single-cell atlas.
#
# Every one of the 63 libraries in the atlas was aligned and quantified against ONE reference. This
# file records it so that a re-run, a new batch, or anything joining external features to the atlas
# uses the same gene universe, and fails loudly if it does not.
#
# Audited 2026-08-19: `refBuild` read from parameters.tsv / dataset.tsv of every CellRanger,
# CellBender, ScSeurat and ScMultiOmics run in FGCZ projects p31662 and p28409. All three stages of
# all five batches' parent chains agree.
#
# source() this at the top of any script that filters, renames or joins genes.

GENOME_BUILD      <- "GRCh38.p13"
ANNOTATION_SOURCE <- "GENCODE"
ANNOTATION_REL    <- "Release_42-2023-01-30"   # GENCODE 42, Ensembl 108
REFBUILD          <- file.path("Homo_sapiens", ANNOTATION_SOURCE, GENOME_BUILD, "Annotation", ANNOTATION_REL)

# FGCZ-internal path. Not reachable from outside FGCZ; the constants above are the portable part.
FGCZ_REFERENCE_DIR <- file.path("/srv/GT/reference", REFBUILD)
FGCZ_FEATURES_FILE <- file.path(FGCZ_REFERENCE_DIR, "Genes", "features_annotation_byGene.txt")
N_GENES_IN_REFERENCE <- 62684L   # data rows in features_annotation_byGene.txt (62,685 lines incl. header)

# Annotation-dependent packages. These are NOT part of the reference above but they change results,
# and nothing else in this repo records them. Versions used for the published figures:
ANNOTATION_PACKAGES <- c(
  "org.Hs.eg.db" = "3.21.0"    # GO/KEGG id mapping in Augur.R, Milo.R, CytoTRACE2.Rmd, DEA_*.R
)
# GWASTools (centromeres.hg38 in CNV_Arm_Analysis.Rmd) is NOT installed in the FGCZ R library, so its
# version could not be measured here and is deliberately not asserted. Whoever runs that Rmd should
# add the version they used rather than have this file claim one.

# NOT the same build: the chromatin modalities. o39391_CellRangerATACCount_2025-08-26 and
# o39575_CellRangerARCCount_2025-09-16 were run on GRCh38.p14 / Release_48-2025-07-03, while
# o37268_CellRangerATACCount_2025-01-07 used the atlas build. Release 42 -> 48 changes the gene set and
# some symbols, so any ATAC/multiome join to the atlas must be checked for alias drift, not assumed.
ATAC_ARC_BUILD <- "Homo_sapiens/GENCODE/GRCh38.p14/Annotation/Release_48-2025-07-03"

#' Assert the loaded annotation packages match the pinned versions.
#' Warns rather than stops by default: an upgrade should be noticed, not silently tolerated, but it
#' should also not block someone re-running a figure on a newer install.
check_annotation_packages <- function(strict = FALSE) {
  for (pkg in names(ANNOTATION_PACKAGES)) {
    want <- ANNOTATION_PACKAGES[[pkg]]
    have <- tryCatch(as.character(packageVersion(pkg)), error = function(e) NA_character_)
    if (is.na(have)) next
    if (!identical(have, want)) {
      msg <- sprintf("%s is %s, figures were produced with %s - enrichment results can move", pkg, have, want)
      if (strict) stop(msg) else warning(msg, call. = FALSE)
    }
  }
  invisible(TRUE)
}

#' Check a Seurat object's features against the pinned reference. FGCZ-only (needs /srv/GT).
#' Returns the fraction of features found, and stops if it is implausibly low - which is what an
#' accidental build mismatch looks like.
check_features_against_reference <- function(object, min_frac = 0.95) {
  if (!file.exists(FGCZ_FEATURES_FILE)) {
    warning("reference not reachable, skipping feature check: ", FGCZ_FEATURES_FILE, call. = FALSE)
    return(invisible(NA_real_))
  }
  ref  <- data.table::fread(FGCZ_FEATURES_FILE, select = "gene_name", showProgress = FALSE)
  feat <- rownames(object)
  if (!length(feat)) stop("object has no rownames - nothing to check against the reference")
  frac <- mean(feat %in% unique(ref$gene_name))
  message(sprintf("%d of %d features (%.2f%%) present in %s",
                  sum(feat %in% unique(ref$gene_name)), length(feat), 100 * frac, ANNOTATION_REL))
  if (frac < min_frac)
    stop(sprintf("only %.2f%% of features match %s - wrong reference build, or symbols were renamed",
                 100 * frac, REFBUILD))
  invisible(frac)
}
