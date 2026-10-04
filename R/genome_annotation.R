# Pinned reference genome and annotation for the UM single-cell atlas.
#
# Every library in the atlas was aligned and quantified against ONE reference. This file records it
# so that a re-run, a new batch, or anything joining external features to the atlas uses the same
# gene universe, and fails loudly if it does not.
#
# source() this at the top of any script that filters, renames or joins genes.

GENOME_BUILD      <- "GRCh38.p13"
ANNOTATION_SOURCE <- "GENCODE"
ANNOTATION_REL    <- "Release_42-2023-01-30"   # GENCODE 42, Ensembl 108
REFBUILD          <- file.path("Homo_sapiens", ANNOTATION_SOURCE, GENOME_BUILD, "Annotation", ANNOTATION_REL)

N_GENES_IN_REFERENCE <- 62684L   # genes in the GENCODE 42 annotation used

# Annotation-dependent packages. These are NOT part of the reference above but they change results,
# and nothing else in this repo records them. Versions used for the published figures:
ANNOTATION_PACKAGES <- c(
  "org.Hs.eg.db" = "3.21.0"    # GO/KEGG id mapping in Augur.R, Milo.R, CytoTRACE2.Rmd, DEA_*.R
)
# GWASTools (centromeres.hg38 in CNV_Arm_Analysis.Rmd) version is deliberately not asserted. Whoever
# runs that Rmd should add the version they used rather than have this file claim one.

# NOT the same build: the chromatin modalities (scATAC, multiome ARC) were run on
# GRCh38.p14 / Release_48-2025-07-03. Release 42 -> 48 changes the gene set and some symbols, so any
# ATAC/multiome join to the atlas must be checked for alias drift, not assumed.
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

#' Check a Seurat object's features against the pinned reference.
#' `features_file`: a gene annotation table of the reference with a `gene_name` column.
#' Returns the fraction of features found, and stops if it is implausibly low - which is what an
#' accidental build mismatch looks like.
check_features_against_reference <- function(object, features_file, min_frac = 0.95) {
  if (!file.exists(features_file)) {
    warning("reference not reachable, skipping feature check: ", features_file, call. = FALSE)
    return(invisible(NA_real_))
  }
  ref  <- data.table::fread(features_file, select = "gene_name", showProgress = FALSE)
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
