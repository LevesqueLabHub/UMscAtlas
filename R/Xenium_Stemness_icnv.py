#!/usr/bin/env python
"""Per-cell iCNV values for Xenium_Stemness.Rmd, read from the stored infercnvpy results (no CNV recomputation).

Boston sections: mean inferred CNV over the chr3 windows and over the chr8 windows of obsm/X_cnv, with cell type and
centroids (Fig. 2h / Suppl. Fig. 5b maps). pUM26 eye / liver: the obs columns chr3_cnv, chr8q_cnv, cnv_rms_score and
centroids stored by the iCNV run (Suppl. Fig. 5a). Only obs, var and obsm/X_cnv are read; existing outputs are kept.

Usage: python -I Xenium_Stemness_icnv.py <xenium_dir> <outputs_dir>
"""
import glob
import os
import sys

import anndata as ad
import h5py
import numpy as np
import pandas as pd

xenium_dir, out_root = sys.argv[1], sys.argv[2]


def window_bounds(chrom, window=100, step=10):
    """Column range of each chromosome in X_cnv. infercnvpy (window_size 100, step 10, chrX/chrY excluded) keeps one
    window every `step` genes per chromosome, or a single window if it has <= `window` genes, chr1-22 in order.
    The stored uns chr_pos is in gene space, so it cannot index X_cnv."""
    n_genes = pd.Series(chrom).value_counts()
    chrs = [f"chr{i}" for i in range(1, 23) if f"chr{i}" in n_genes]
    n_win = [int(np.ceil((n_genes[c] - window + 1) / step)) if n_genes[c] > window else 1 for c in chrs]
    end = np.cumsum(n_win)
    return {c: (int(e - n), int(e)) for c, n, e in zip(chrs, n_win, end)}


def chr_means(g, bounds, chunk=50_000):
    """Per-cell mean of a CSR matrix over each column range (implicit zeros count as 0)."""
    assert g.attrs["encoding-type"] == "csr_matrix"
    n = int(g.attrs["shape"][0])
    indptr = g["indptr"][:].astype(np.int64)
    out = {c: np.empty(n) for c in bounds}
    for a in range(0, n, chunk):
        b = min(a + chunk, n)
        data, idx = g["data"][indptr[a]:indptr[b]], g["indices"][indptr[a]:indptr[b]]
        row = np.repeat(np.arange(b - a), np.diff(indptr[a:b + 1]))
        for c, (s, e) in bounds.items():
            m = (idx >= s) & (idx < e)
            out[c][a:b] = np.bincount(row[m], weights=data[m], minlength=b - a) / (e - s)
    return out


# Boston: chr3 and chr8 per cell
boston_csv = f"{out_root}/Xenium_Stemness_Boston_2026b/cnv_chr3_chr8_per_cell.csv"
if not os.path.exists(boston_csv):
    rows = []
    for p in sorted(glob.glob(f"{xenium_dir}/boston/cnv/*_cnv.h5ad")):
        with h5py.File(p, "r") as h:
            obs, var = ad.io.read_elem(h["obs"]), ad.io.read_elem(h["var"])
            wb = window_bounds(var["chromosome"])
            assert wb["chr22"][1] == h["obsm/X_cnv"].attrs["shape"][1], "window count does not match X_cnv"
            m = chr_means(h["obsm/X_cnv"], {c: wb[c] for c in ("chr3", "chr8")})
        rows.append(pd.DataFrame({"sample": os.path.basename(p).removesuffix("_cnv.h5ad"), "barcode": obs.index,
                                  "cell_type": obs["cell_type"].astype(str).values, "x": obs["x_centroid"].values,
                                  "y": obs["y_centroid"].values, "chr3_cnv": m["chr3"], "chr8_cnv": m["chr8"]}))
        print(rows[-1]["sample"].iloc[0], len(obs), "cells")
    os.makedirs(os.path.dirname(boston_csv), exist_ok=True)
    pd.concat(rows, ignore_index=True).to_csv(boston_csv, index=False)
    print("wrote", boston_csv)

# pUM26: stored per-cell iCNV columns with the barcode
cols = ["cell_type", "refined_immune_type", "chr3_cnv", "chr8q_cnv", "cnv_rms_score", "x_centroid", "y_centroid"]
for pair in ["pUM26_eye", "pUM26_liver"]:
    out = f"{out_root}/Xenium_Stemness_pUM26_2026b/{pair}_icnv_obs.csv"
    if os.path.exists(out):
        continue
    with h5py.File(f"{xenium_dir}/pUM26/{pair}/icnv_cnv.h5ad", "r") as h:
        obs = ad.io.read_elem(h["obs"])
    obs = obs[cols]
    obs.index.name = "barcode"
    os.makedirs(os.path.dirname(out), exist_ok=True)
    obs.to_csv(out)
    print("wrote", out)
