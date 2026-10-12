# Is the reference origin/BAP1 factor (out/, K=6 seed 42; REF_FACTOR env var) recovered across seeds 1-5 x K in {4,6,8}?
# Per refit: the factor best matching the reference by sample scores, its |r|, its covariate padj,
# the within-primary BAP1 test, and gene-loading agreement per view. October rerun of mc_astra/, 2026-10-05.
import os
from pathlib import Path
import numpy as np, pandas as pd, anndata as ad
from scipy.stats import mannwhitneyu

D = Path("../mc_astra_2026b")
REF = os.environ["REF_FACTOR"]

def loadings(run_dir):
    a = ad.read_h5ad(run_dir / "amodel.h5ad")
    return pd.DataFrame(a.varm["gene_loadings"], columns=a.uns["gene_loadings_columns"], index=a.var_names)

def primary_bap1_p(f, col):
    p = f[f.origin == "primary"]
    b = p.loc[p.secondary_mutation == "BAP1", col]
    o = p.loc[p.secondary_mutation.isin(["SF3B1", "EIF1AX", "none"]), col]
    return mannwhitneyu(b, o).pvalue

ref_f = pd.read_csv(D / "out/factor_scores.csv", index_col=0)
# July Factor2 (the July origin/BAP1 axis), keyed on the short sample name: does any refit recover it?
import re
jul = pd.read_csv(D.parent / "mc_astra/out/factor_scores.csv", index_col=0)["Factor2"]
jul.index = [re.sub(r"_GEX.*$", "", re.sub(r"^\d+_\d+-", "", s)) for s in jul.index]
ref_w = loadings(D / "out")
# loading columns look like "<view>:<gene>" (append_view_to_var); views compared separately
views = sorted({c.split(":")[0] for c in ref_w.columns})

rows = []
for run in sorted((D / "stability").glob("K*_seed*")):
    if not (run / "factor_covariate_padj.csv").exists():
        rows.append({"run": run.name, "status": "missing"}); continue
    f = pd.read_csv(run / "factor_scores.csv", index_col=0)
    fac = [c for c in f.columns if c.startswith("Factor")]
    shared = ref_f.index.intersection(f.index)
    r = f.loc[shared, fac].corrwith(ref_f.loc[shared, REF])
    best = r.abs().idxmax()
    sign = np.sign(r[best])
    pm = pd.read_csv(run / "factor_covariate_padj.csv", index_col=0)
    w = loadings(run)
    row = {"run": run.name, "status": "ok", "n_samples": len(f), "best": best, "r_scores": round(r[best], 3),
           "second_best_abs_r": round(r.abs().drop(best).max(), 3),
           "padj_origin": pm.loc[best, "origin"], "padj_BAP1": pm.loc[best, "secondary_BAP1"],
           "padj_batch": pm.loc[best, "batch"], "padj_smut_pUM": pm.loc[best, "smut_pUM"],
           "padj_log_ncells": pm.loc[best, "log_ncells"], "p_primary_BAP1": primary_bap1_p(f, best)}
    for v in views:
        cols = [c for c in ref_w.columns if c.startswith(v + ":") and c in w.columns]
        row[f"r_load_{v}"] = round(np.corrcoef(ref_w.loc[REF, cols], sign * w.loc[best, cols])[0, 1], 3)
    fs = f.set_index("short_name").loc[lambda d: d.index.isin(jul.index), fac]
    rj = fs.corrwith(jul.loc[fs.index])
    row.update({"july_F2_best": rj.abs().idxmax(), "july_F2_abs_r": round(rj.abs().max(), 3)})
    rows.append(row)

res = pd.DataFrame(rows)
res.to_csv(D / "stability_summary.csv", index=False)
pd.set_option("display.width", 250); pd.set_option("display.max_columns", 30)
print(res.round(4).to_string())
ok = res[res.status == "ok"]
print(f"\n{len(ok)}/{len(res)} refits ok; |r_scores| >= 0.8 in {(ok.r_scores.abs() >= 0.8).sum()}; "
      f"origin padj < 0.05 in {(ok.padj_origin < 0.05).sum()}; BAP1 padj < 0.05 in {(ok.padj_BAP1 < 0.05).sum()}; "
      f"within-primary BAP1 p < 0.05 in {(ok.p_primary_BAP1 < 0.05).sum()}")
