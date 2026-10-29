#!/usr/bin/env python
"""14c step 1: per-cell iCNV columns for pUM26 eye + liver, WITH the cell barcode.

The report_cache/*_obs.csv files written by Xenium_pUM26/SPLIT_Reannotation/21_extract_cnv_obs.py
dropped the obs index (index=False), so they cannot be joined to other per-cell data. This reads only
obs from the 1.8 / 8.6 GB *_cnv.h5ad files (backed mode, X and obsm untouched) and keeps the barcode.
Run: python 14c_extract_pum26_icnv_obs.py
"""
import os
import anndata as ad

src = "../Xenium_pUM26/SPLIT_Reannotation/xenium_icnv"
out = "../Outputs/Xenium_Stemness_pUM26_2026b"
os.makedirs(out, exist_ok=True)
# x/y_centroid are the coordinates 17_xenium_icnv.R recovered from the original (pre-SPLIT) region objects
cols = ["cell_type", "refined_immune_type", "chr3_cnv", "chr8q_cnv", "cnv_rms_score", "x_centroid", "y_centroid"]
for pair in ["pUM26_eye", "pUM26_liver"]:
    a = ad.read_h5ad(f"{src}/{pair}_cnv.h5ad", backed="r")
    df = a.obs[[c for c in cols if c in a.obs.columns]].copy()
    df.index.name = "barcode"
    df.to_csv(f"{out}/{pair}_icnv_obs.csv")
    print(pair, df.shape, "Melanoma:", int((df["cell_type"] == "Melanoma").sum()), "index head:", list(df.index[:2]))
    a.file.close()
