#!/usr/bin/env python
"""14d: per-cell CNV RMS = sqrt(mean(X_cnv^2)) over the CNV windows, as in the submitted Fig. 2h (B7).

The 09 score (cnv_scores_09method_all_samples.csv) is one value per CNV cluster, so a per-cell scatter
against it shows vertical stripes; the manuscript plotted a per-cell RMS. obsm/X_cnv is CSR in every
file, read here in row chunks with h5py (X and the rest of the 0.6-8.6 GB files are never loaded).
Implicit zeros count as 0 in the mean (denominator = number of windows).
Positive controls, printed and asserted:
  - Boston: per-cell RMS averaged per CNV cluster (= unique 09 score within a section, melanoma cells)
    must correlate with the 09 per-cluster mean |CNV|.
  - pUM26: must equal obs cnv_rms_score, which script 18 stored with the iCNV run.
Run: python -I 14d_cnv_rms_per_cell.py   (~2 min)
"""
import h5py
import numpy as np
import pandas as pd

boston = "../Xenium_Boston_Sullivan"
liver, breast = f"{boston}/06_InSituCNV_Outputs", f"{boston}/07_InSituCNV_Breast_Outputs"
pum26 = "../Xenium_pUM26/SPLIT_Reannotation/xenium_icnv"
out_csv = "../Outputs/Xenium_Stemness_Boston_2026b/cnv_rms_per_cell.csv"
paths = {  # the files compute_cnv_scores_09method.py reads, plus pUM26
    "S22-19XXX-A1_Slide12": f"{liver}/S22-19XXX-A1_Slide12_cnv.h5ad",
    "S22-40XXX-A1_Slide12": f"{liver}/S22-40XXX-A1_Slide12_cnv.h5ad",
    "S22-87XX-A2_Slide11": f"{liver}/S22-87XX-A2_Slide11_cnv.h5ad",
    "S23-32XXX-A1_Slide11": f"{breast}/S23-32XXX-A1_Slide11_cnv.h5ad",
    "S23-70XXX-A1_Slide12": f"{liver}/S23-70XXX-A1_Slide12_cnv.h5ad",
    "S23-89X-A_Slide11": f"{liver}/S23-89X-A_Slide11_cnv.h5ad",
    "pUM26_eye": f"{pum26}/pUM26_eye_cnv.h5ad",
    "pUM26_liver": f"{pum26}/pUM26_liver_cnv.h5ad",
}


def rms_rows(g, chunk=50_000):
    assert g.attrs["encoding-type"] == "csr_matrix", g.attrs["encoding-type"]
    n, ncol = (int(v) for v in g.attrs["shape"])
    indptr = g["indptr"][:].astype(np.int64)
    out = np.empty(n)
    for a in range(0, n, chunk):
        b = min(a + chunk, n)
        sq = g["data"][indptr[a]:indptr[b]] ** 2
        row = np.repeat(np.arange(b - a), np.diff(indptr[a:b + 1]))
        out[a:b] = np.sqrt(np.bincount(row, weights=sq, minlength=b - a) / ncol)
    return out


rows, ref = [], {}
for s, p in paths.items():
    with h5py.File(p, "r") as h:
        bc = h["obs/_index"].asstr()[:]
        rms = rms_rows(h["obsm/X_cnv"])
        if "cnv_rms_score" in h["obs"]:
            ref[s] = h["obs/cnv_rms_score"][:]
    assert len(bc) == len(rms) and len(set(bc)) == len(bc)
    rows.append(pd.DataFrame({"sample": s, "barcode": bc, "cnv_rms": rms}))
    print(s, len(bc), "cells, RMS median", round(float(np.median(rms)), 5))
res = pd.concat(rows, ignore_index=True)

# Control 1 (pUM26): identical to the RMS stored by the iCNV run
for s, v in ref.items():
    d = float(np.max(np.abs(res.loc[res["sample"] == s, "cnv_rms"].to_numpy() - v)))
    print(f"control pUM26 {s}: max |RMS - obs cnv_rms_score| = {d:.2e}")
    assert d < 1e-6

# Control 2 (Boston): cluster mean of per-cell RMS vs the 09 per-cluster score
c09 = pd.read_csv(f"{liver}/cnv_scores_09method_all_samples.csv").rename(columns={"cell_barcode": "barcode"})
m = res.merge(c09, on=["sample", "barcode"], how="inner")
assert len(m) == len(c09), (len(m), len(c09))
cl = m[m["cnv_score"] > 0].groupby(["sample", "cnv_score"], as_index=False)["cnv_rms"].mean()
per = cl.groupby("sample")[["cnv_score", "cnv_rms"]].apply(lambda x: pd.Series({"k": len(x), "r": np.corrcoef(x["cnv_score"], x["cnv_rms"])[0, 1]}))
r_all = np.corrcoef(cl["cnv_score"], cl["cnv_rms"])[0, 1]
print(per.round(3).to_string())
print(f"control Boston: Pearson r over all {len(cl)} CNV clusters = {r_all:.3f} (per section {per['r'].min():.2f} to {per['r'].max():.2f})")
assert r_all > 0.8

res.to_csv(out_csv, index=False)
print("Saved", out_csv, len(res), "rows")
