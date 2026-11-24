#!/usr/bin/env python
"""pUM26 eye / liver iCNV without neighbour smoothing, per-cell RMS for the Suppl. Fig. 5c scatter (Xenium_Stemness.Rmd).

Same as the iCNV run of these sections (infercnvpy 0.6, normalize_total + log1p, macrophage reference, window 100,
step 10, dynamic threshold 1.5) but without the k = 100 neighbour smoothing, like the Boston run of Fig. 2i.
Writes only the per-cell scores. Long-running (infercnvpy on all cells); run once per section before the Rmd.

Usage: python Xenium_Stemness_unsmoothed_icnv.py <xenium_dir> <outputs_dir> <pUM26_eye|pUM26_liver>
"""
import sys

import infercnvpy as cnv
import numpy as np
import scanpy as sc
from scipy.sparse import issparse

xenium_dir, out_root, pair = sys.argv[1:4]

adata = sc.read_h5ad(f"{xenium_dir}/pUM26/{pair}/icnv_input.h5ad")
adata = adata[:, adata.var[["chromosome", "start", "end"]].notna().all(axis=1)].copy()
adata.var["start"] = adata.var["start"].astype(int)
adata.var["end"] = adata.var["end"].astype(int)
adata.var["chromosome"] = ("chr" + adata.var["chromosome"].astype(str)).astype("category")
ref = adata.uns["reference_cat"]
ref = [ref] if isinstance(ref, str) else list(ref)

sc.pp.normalize_total(adata)
sc.pp.log1p(adata)
cnv.tl.infercnv(adata, reference_key="cell_type", reference_cat=ref, window_size=100, step=10, dynamic_threshold=1.5)

X = adata.obsm["X_cnv"]
X = X.toarray() if issparse(X) else X
obs = adata.obs[["cell_type"]].copy()
obs.index.name = "barcode"
obs["cnv_rms_unsmoothed"] = np.sqrt(np.mean(X ** 2, axis=1))
obs.to_csv(f"{out_root}/Xenium_Stemness_pUM26_2026b/{pair}_icnv_unsmoothed.csv")
