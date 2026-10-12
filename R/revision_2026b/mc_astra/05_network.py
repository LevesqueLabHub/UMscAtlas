# Origin/BAP1 factor (FACTOR env var) of the reference fit (out/): multicellular coordination network and the
# ligand-receptor pairs that could underlie it, following the tail of GetStarted_MOFA.ipynb.
# October rerun of mc_astra/05_network.py, 2026-10-05. The sign of "positive" follows the factor's own scores.
import warnings
warnings.filterwarnings("ignore")
import os
from pathlib import Path
import pandas as pd, anndata as ad
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import liana as li, mc_astra as mca

D = Path("../mc_astra_2026b")
FAC = os.environ["FACTOR"]
# Network OLS uses complete cases only. In July all 6 views left 14 samples (B/Plasma the bottleneck) and
# dropping B/Plasma left 31. DROP_VIEWS (comma-separated) removes views; VIEWSET names the output dir.
VIEWSET = os.environ.get("VIEWSET", "all")
DROP = [v for v in os.environ.get("DROP_VIEWS", "").split(",") if v]
OUT = D / f"out/network_{FAC}_{VIEWSET}"; OUT.mkdir(parents=True, exist_ok=True)

amodel = ad.read_h5ad(D / "out/amodel.h5ad")
all_views = sorted(amodel.obsm.keys())
VIEWS = None if not DROP else [v for v in all_views if v not in DROP]
assert all(v in all_views for v in DROP), (DROP, all_views)
used = VIEWS or all_views
cc = amodel.obs[[f"{v}_n_cells" for v in used]].notna().all(axis=1)
pd.DataFrame({"viewset": [VIEWSET], "views": [";".join(used)], "n_complete": [int(cc.sum())],
              **{f"n_{o}": [int((amodel.obs.origin[cc] == o).sum())] for o in ["non-malignant", "primary", "metastasis"]}}
             ).to_csv(OUT / "complete_samples.csv", index=False)
print("complete samples:", int(cc.sum()), amodel.obs.origin[cc].value_counts().to_dict())
gl = pd.DataFrame(amodel.varm["gene_loadings"], columns=amodel.uns["gene_loadings_columns"], index=amodel.var_names)
gl = mca.down.split_by_view(gl)
if VIEWS: gl = {k: v for k, v in gl.items() if k in VIEWS}

mins = mca.down.get_multicell_net(test_model=amodel, sel_factor=FAC, views=VIEWS, random_effect=None,
                                  standardize=True, drop_na=True, verbose=True)
resource = li.rs.select_resource("consensus"); resource["lr_exists"] = 1

for side in ["positive", "negative"]:
    net = mins[side]
    net.to_csv(OUT / f"network_{side}.csv", index=False)
    print(f"\n== {side}: {len(net)} edges\n", net.sort_values("cor_estimate", key=abs, ascending=False).head(12).to_string())
    mca.pl.plot_mcell_network(df=net, weight_col="cor_estimate", abs_cutoff=0.2, keep_negative=False,
                              title=f"{FAC} {side} - coordination", reciprocal_curvature=0.3, default_curvature=0.4)
    plt.savefig(OUT / f"network_{side}.png", dpi=200, bbox_inches="tight"); plt.close("all")
    mca.pl.plot_interaction_tileplot(df=net, value_col="cor_estimate", row_col="target", col_col="predictor",
                                     cluster_by="rows", same_order_for_rows_cols=True, fill_missing_for_clustering=0.0,
                                     cmap="RdBu_r", center=0.0, show_values=True, value_decimals=2, figsize=(8, 6))
    plt.savefig(OUT / f"tiles_{side}.png", dpi=200, bbox_inches="tight"); plt.close("all")

    lr = mca.down.lr_usage(gene_loadings=gl, resource=resource, network_df=net.copy(), sel_factor=FAC,
                           weight_col="cor_estimate", abs_cutoff=0, keep_negative=False, loading_type=side, n_top=None)
    lr.to_csv(OUT / f"lr_{side}.csv", index=False)
    print(f"-- {side}: {len(lr)} L-R events\n", lr.head(20).to_string())
    if len(lr):
        mca.pl.plot_comm_overview(plot_df=lr, tile_width=0.6, tile_height=0.6, text_size=5, figsize=(8, 5))
        plt.savefig(OUT / f"lr_overview_{side}.png", dpi=200, bbox_inches="tight"); plt.close("all")
        mca.pl.plot_lr_circos(plot_df=lr, link_width=2, arrow_height=8, arrow_width=8, gene_label_size=7,
                              sector_gap=10, figsize=(6, 6), sector_r_lim=(80, 100))
        plt.savefig(OUT / f"lr_circos_{side}.png", dpi=200, bbox_inches="tight"); plt.close("all")
print("done")
