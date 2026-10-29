#!/usr/bin/env python
"""14e: pUM26 eye + liver iCNV WITHOUT neighbour smoothing, the way the January Boston run was made.

Paul, 2026-10-20: Victoria wants Fig. 2h on the January (unsmoothed) InSituCNV run; pUM26 (Suppl. Fig. 5b)
had only ever been run with the Jensen et al. k = 100 smoothing (Xenium_pUM26/SPLIT_Reannotation/
18_xenium_icnv.py, 2026-07-02). This is that script minus the smoothing block: same input h5ad, same genes,
normalize_total + log1p, same reference label, window_size 100, step 10, dynamic_threshold 1.5, same RMS and
chr3 / chr8q (second half of chr8 windows) scores. Writes only the per-cell scores (not an 8 GB h5ad).
Run: 14e_pum26_icnv_unsmoothed.sbatch (gi_insitucnv env).
"""
import sys
import numpy as np
import scanpy as sc
import infercnvpy as cnv
from scipy.sparse import issparse

pair = sys.argv[1]
src = "../Xenium_pUM26/SPLIT_Reannotation/xenium_icnv"
out = "../Outputs/Xenium_Stemness_pUM26_2026b"

adata = sc.read_h5ad(f"{src}/{pair}.h5ad")
keep = adata.var[["chromosome", "start", "end"]].notna().all(axis=1)
adata = adata[:, keep].copy()
adata.var["start"] = adata.var["start"].astype(int)
adata.var["end"] = adata.var["end"].astype(int)
adata.var["chromosome"] = ("chr" + adata.var["chromosome"].astype(str)).astype("category")
ref = adata.uns["reference_cat"]
ref = [ref] if isinstance(ref, str) else list(ref)
print(f"{pair}: {adata.n_obs} cells, {adata.n_vars} genes, reference={ref}, n_ref={int(adata.obs['cell_type'].isin(ref).sum())}")

sc.pp.normalize_total(adata)
sc.pp.log1p(adata)
cnv.tl.infercnv(adata, reference_key="cell_type", reference_cat=ref, window_size=100, step=10, dynamic_threshold=1.5)

X = adata.obsm["X_cnv"]
X = X.toarray() if issparse(X) else X
chr_pos = adata.uns["cnv"]["chr_pos"]
n_win = X.shape[1]
assert max(chr_pos.values()) < n_win, "chr_pos is not in window space"
order = sorted(chr_pos.items(), key=lambda kv: kv[1])
win_chr = np.empty(n_win, dtype=object)
for i, (c, s) in enumerate(order):
    win_chr[s:(order[i + 1][1] if i + 1 < len(order) else n_win)] = c


def arm(c, a=0.0, b=1.0):
    idx = np.where(win_chr == c)[0]
    return X[:, idx[int(a * len(idx)):int(np.ceil(b * len(idx)))]].mean(axis=1)


obs = adata.obs[["cell_type"]].copy()
obs.index.name = "barcode"
obs["cnv_rms_unsmoothed"] = np.sqrt(np.mean(X ** 2, axis=1))
obs["chr3_cnv_unsmoothed"] = arm("chr3")
obs["chr8q_cnv_unsmoothed"] = arm("chr8", 0.5)
mel = obs["cell_type"] == "Melanoma"
print(f"{pair}: melanoma {int(mel.sum())}; mean chr3 {obs.loc[mel, 'chr3_cnv_unsmoothed'].mean():.4f} "
      f"(others {obs.loc[~mel, 'chr3_cnv_unsmoothed'].mean():.4f}); mean chr8q {obs.loc[mel, 'chr8q_cnv_unsmoothed'].mean():.4f} "
      f"(others {obs.loc[~mel, 'chr8q_cnv_unsmoothed'].mean():.4f})")
obs.to_csv(f"{out}/{pair}_icnv_unsmoothed.csv")
print("wrote", f"{out}/{pair}_icnv_unsmoothed.csv")
