# UM Atlas - mc-ASTRA (MOFAcell) fit on the October 2026 object, following docs/notebooks/GetStarted_MOFA.ipynb.
# Rerun of mc_astra/02_mofacell.py (July) with the same parameters. Input: pb/ from 01_pseudobulk.R. 2026-10-05.
# Changes vs July: batch = Condition (Batch01-05) replaces the library-prefix "run" (checked 1:1 below);
# scoped covariates added (secondary mutation within pUM / choroidal pUM, location within pUM / non-malignant);
# per-view R2 saved as variance_explained.csv.
import warnings
warnings.filterwarnings("ignore")
import math, os, re
from pathlib import Path
import numpy as np, pandas as pd, scipy.io, anndata as ad, mudata as md
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import scanpy as sc, decoupler as dc, mofaflex as mf, mc_astra as mca

# ponytail: monkeypatch MOFA-FLEX (main @ 0.2.0rc2) post-fit mean imputation. IterativeImputer over a
# 1-group x p-feature matrix fits p BayesianRidge models x 10 iters even when nothing is missing
# (~40 of the 53 min of the first full run), and is the identity then. With one group and NaNs it
# silently DROPS those columns, so fill them with the row mean instead. Drop this if fixed upstream.
from mofaflex._core.likelihoods import base as _mf_lik
_orig_impute = _mf_lik.Likelihood._impute_feature_summary_statistics
def _fast_impute(self, stats, missing_values=np.nan):
    m = np.stack(tuple(stats.values()), axis=0)
    if not np.isnan(m).any():
        return dict(stats)
    if m.shape[0] == 1:
        return {k: np.where(np.isnan(v), np.nanmean(v), v) for k, v in stats.items()}
    return _orig_impute(self, stats, missing_values)
_mf_lik.Likelihood._impute_feature_summary_statistics = _fast_impute

D = Path("../mc_astra_2026b")
# Defaults = the reference fit in out/; 03_stability.sh overrides all three per array task.
OUT = Path(os.environ.get("OUT_DIR", D / "out")); OUT.mkdir(parents=True, exist_ok=True)
N_FACTORS = int(os.environ.get("N_FACTORS", 6))
SEED = int(os.environ.get("SEED", 42))
MAX_EPOCHS = int(os.environ.get("MAX_EPOCHS", 10_000))   # smoke test: MAX_EPOCHS=50

# ---- pseudobulk AnnData (groups x genes) ----
# Dense on purpose: mc_astra.up.filt.filter_genes_byexpr breaks on sparse X (np.matrix lib sizes);
# dc.pp.pseudobulk in the tutorial returns dense too.
X = scipy.io.mmread(D / "pb/counts.mtx").T.toarray().astype(np.float32)
obs = pd.read_csv(D / "pb/obs.csv")
obs["celltype"] = obs.celltype.str.replace("/", "_")   # h5ad keys cannot hold "/" (B/Plasma, T/NK)
genes = (D / "pb/genes.tsv").read_text().splitlines()
pdata = ad.AnnData(X=X, obs=obs.set_index("group"), var=pd.DataFrame(index=genes))
pdata.obs.index.name = None

meta = pd.read_csv(D / "pb/samples.csv").set_index("Sample", drop=False)
meta.index.name = None
# Sequencing batch. July used the order prefix of the library name ("run"); in October every library
# carries a prefix and Condition (Batch01-05) is the lab's own batch column. Fail if they disagree.
meta["run"] = [m.group(1) if (m := re.match(r"^(o?\d+)_", s)) else "none" for s in meta.Sample]
meta["batch"] = meta.Condition
assert (meta.groupby("run").batch.nunique() == 1).all() and (meta.groupby("batch").run.nunique() == 1).all(), \
    meta.groupby(["run", "batch"]).size()
meta["enrich16S"] = meta.Sample.str.contains("_16S_").map({True: "16S", False: "std"})
meta["secondary_BAP1"] = np.where(meta.secondary_mutation == "BAP1", "BAP1", "other")
meta["log_ncells"] = np.log10(meta.n_cells)
# Scoped covariates: NaN outside their scope, so get_associations tests only the samples they apply to.
prim = meta.origin == "primary"
called = meta.secondary_mutation.isin(["BAP1", "SF3B1", "EIF1AX"])
meta["smut_pUM"] = meta.secondary_mutation.where(prim & called)
meta["smut_choroid_pUM"] = meta.secondary_mutation.where(prim & called & (meta.primary_location == "choroid"))
meta["location_pUM"] = meta.primary_location.where(prim & meta.primary_location.isin(["choroid", "iris", "ciliary body"]))
meta["location_nonmal"] = meta.healthy_location.where(meta.origin == "non-malignant")
print(meta.groupby(["batch", "origin"]).size().unstack(fill_value=0))

# ---- views: one per cell type, obs = Sample ----
contaminants = mca.up.get_contaminant_genes(pdata, view_group="celltype", min_count=10, min_total_count=20,
                                            large_n=10, min_prop=0.4, lfc_threshold=3, p_threshold=0.001)
pdata.obs = pdata.obs.set_index("Sample", drop=False); pdata.obs.index.name = None
views = mca.up.split_anndata_by_celltype(pdata=pdata, grouping="celltype")
for v in views.values():
    v.var["total_counts"] = np.asarray(v.X.sum(axis=0)).ravel()

min_samples = math.trunc(len(meta) * 0.4)   # a view must cover >= 40 % of samples (25/63)
mca.up.filt.filter_anndata_by_ncells(views, min_cells=20)
mca.up.filt.filter_views_by_samples(views, min_rows=min_samples)
mca.up.filt.filter_genes_byexpr(views, min_count=20, min_prop=0.4)
mca.up.filt.filter_views_by_genes(views, min_genes_per_view=100)
mca.up.filt.filter_samples_by_coverage(views, threshold=0, min_prop=0.90)
mca.up.filt.filter_views_by_samples(views, min_rows=min_samples)
mca.up.utils.save_raw_counts(views, layer_name="raw_counts")
mca.up.norm_log(views, target_sum=1e6, exclude_highly_expressed=False, max_value=None, center=True)
mca.up.filter_hvgs(views, groupby=None, ngroups_cut=None)
mca.up.filter_smpls_by_nview(views, min_views=4)   # tutorial uses 5; with 6 views that drops 9/16 non-malignant (no Melanoma/B view)
mca.up.filter_genes_by_celltype(views, {k: v for k, v in contaminants.items() if k in views})
summary = pd.DataFrame({k: {"samples": v.n_obs, "genes": v.n_vars} for k, v in views.items()}).T
print(summary); summary.to_csv(OUT / "views_kept.csv")

# ---- fit ----
mca.up.utils.append_view_to_var(views)
mdata = md.MuData(views)
model = mf.terms.MofaFlex(n_factors=N_FACTORS, weight_prior="SpikeSlab", init_factors="pca")
model.fit(mdata, seed=SEED, save_path=str(OUT / "mofaflex_model.h5"), lr=0.01, early_stopper_patience=1000,
          likelihoods="Normal", max_epochs=MAX_EPOCHS, plot_data_overview=False)
amodel = mca.down.model_to_anndata(anndata_dict=views, metadata=meta, model=model)
amodel.write_h5ad(OUT / "amodel.h5ad")
pd.DataFrame(amodel.X, index=amodel.obs_names, columns=[f"Factor{i+1}" for i in range(amodel.n_vars)]) \
  .join(meta).to_csv(OUT / "factor_scores.csv")

# Per-view R2 of every factor (amodel.var), for the report.
amodel.var.rename(columns=lambda c: c.split(":")[0]).to_csv(OUT / "variance_explained.csv")
try:
    ve = mf.pl.variance_explained(model)
    (ve.save(OUT / "variance_explained.png", dpi=200) if hasattr(ve, "save") else plt.savefig(OUT / "variance_explained.png", dpi=200, bbox_inches="tight"))
except Exception as e:
    print("variance_explained plot failed:", e)
plt.close("all")

# ---- which covariates do factors track? technical first, then biology ----
cov = ["batch", "enrich16S", "log_ncells", "origin", "organ", "primary_mutation", "secondary_BAP1", "smut_pUM",
       "smut_choroid_pUM", "location_pUM", "location_nonmal", "Gender", "Age_atDiagnosisY"]
pmat = mca.down.get_pval_matrix(amodel, cov)
pmat.to_csv(OUT / "factor_covariate_padj.csv"); print(pmat)
mca.pl.plot_pval_tiles(pmat, title=f"MOFA {N_FACTORS} factors - adj. p", star_threshold=0.05)
plt.savefig(OUT / "factor_covariate_padj.png", dpi=200, bbox_inches="tight"); plt.close("all")

for c in ["origin", "batch"]:
    st = mca.down.get_associations(amodel, test_variable=c, test_type="categorical", random_effect=None)
    st.to_csv(OUT / f"assoc_{c}.csv")
    try:
        pd.Series(mca.down.calc_total_variance(amodel, st, pval_thrs=0.05)).to_csv(OUT / f"total_var_{c}.csv")
    except Exception as e:
        print(f"calc_total_variance {c} failed:", e)
    sc.pl.violin(amodel, [f"Factor{i+1}" for i in range(N_FACTORS)], groupby=c, rotation=90, size=3, show=False)
    plt.savefig(OUT / f"violin_{c}.png", dpi=150, bbox_inches="tight"); plt.close("all")

# ---- PROGENy on gene loadings, per view ----
gl = pd.DataFrame(amodel.varm["gene_loadings"], columns=amodel.uns["gene_loadings_columns"], index=amodel.var_names)
gl = mca.down.split_by_view(gl)
pws = mca.down.run_ulm_per_view(view_dict=gl, net=dc.op.progeny(organism="human", top=500))
for k, v in pws.items():
    for kk, vv in v.items():
        if isinstance(vv, pd.DataFrame): vv.to_csv(OUT / f"progeny_{k.replace('/', '-').replace(' ', '_')}_{kk}.csv")
mca.pl.plot_mcell_funcomics(result_dict=pws, result_key="pw_acts", pval_key="pw_padj", p_threshold=0.05, top_n=5,
                            cmap="coolwarm", figsize=(14, 3), ytick_rotation=0, use_var=True, share_color_scale=True, center=0.0)
plt.savefig(OUT / "progeny_per_view.png", dpi=200, bbox_inches="tight"); plt.close("all")
print("done")
