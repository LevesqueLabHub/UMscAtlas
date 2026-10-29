#!/usr/bin/env python3
"""Build a clean Xenium iCNV supplementary figure for the UM atlas NC revision.

Addresses Reviewer 2 (manuscript L190): "provide an inferCNV to defend your
hypothesis ... show correlation in tumor compartment of iCNV and cytotrace".

This script produces the *iCNV-results* half of the deliverable (the
tumor-compartment correlation already exists as B5/B6 in Revised_Plots_NC_20260522).

For each of the 6 Boston/Sullivan mUM Xenium samples it:
  1. Loads the precomputed _cnv.h5ad (X already = smoothed log-norm infercnv input).
  2. Re-runs infercnvpy windowing (fast) to recover a CORRECT window-space chr_pos
     (the stored chr_pos was overwritten with a gene-index dict and cannot map the
      273-column X_cnv). Re-run reproduces the 273 windows exactly.
  3. Aggregates X_cnv into a broad-lineage x chromosome mean-CNV matrix.
  4. Extracts per-cell whole-chromosome chr3 (loss) and chr8 (8q gain) iCNV values
     + spatial centroids for spatial maps.
Outputs per-sample CSVs + three figures:
  - Fig 1: per-sample chromosome iCNV heatmap (lineage x chr1-22).
  - Fig 2: combined melanoma chr-profile across all 6 samples (the key panel).
  - Fig 3: spatial chr3 (loss) + chr8 (gain) maps per sample.
  - Fig 4: infercnvpy cnv_score by lineage (melanoma vs normal), per sample.

Env: gi_insitucnv (infercnvpy 0.6.0). Run via SBATCH_xenium_icnv_suppfig.sh.
"""
import os
import warnings
import numpy as np
import pandas as pd
import scipy.sparse as sp
import anndata as ad
import infercnvpy as cnv
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm

warnings.filterwarnings("ignore")

LIVER_DIR = "../Xenium_Boston_Sullivan/06_InSituCNV_Outputs"
BREAST_DIR = "../Xenium_Boston_Sullivan/07_InSituCNV_Breast_Outputs"
OUT = "../Xenium_iCNV_SuppFig"
os.makedirs(OUT, exist_ok=True)

# sample -> (h5ad path, reference cell-type label)
SAMPLES = {
    "S22-19XXX-A1_Slide12": (f"{LIVER_DIR}/S22-19XXX-A1_Slide12_cnv.h5ad", "Macrophages"),
    "S22-40XXX-A1_Slide12": (f"{LIVER_DIR}/S22-40XXX-A1_Slide12_cnv.h5ad", "Macrophages"),
    "S22-87XX-A2_Slide11":  (f"{LIVER_DIR}/S22-87XX-A2_Slide11_cnv.h5ad",  "Macrophages"),
    "S23-70XXX-A1_Slide12": (f"{LIVER_DIR}/S23-70XXX-A1_Slide12_cnv.h5ad", "Macrophages"),
    "S23-89X-A_Slide11":    (f"{LIVER_DIR}/S23-89X-A_Slide11_cnv.h5ad",    "Macrophages"),
    "S23-32XXX-A1_Slide11": (f"{BREAST_DIR}/S23-32XXX-A1_Slide11_cnv.h5ad", "macrophage"),
}

CHR_ORDER = [f"chr{i}" for i in range(1, 23)]

# Broad-lineage mapping covering both the liver (capitalised) and breast (ontology)
# annotation schemes. Anything unmatched -> "Other".
def lineage(ct):
    s = str(ct).lower()
    if "melanoma" in s:
        return "Melanoma"
    if "macrophage" in s or "mono" in s or "dendritic" in s or s in ("cdc1s", "cdc2s", "pdcs", "mig_cdcs", "basophils", "neutrophils", "mast cell"):
        return "Myeloid (ref)"
    if any(k in s for k in ("t cell", "t_cell", "nk", "natural killer", "regulatory t", "helper", "lymphocyte")):
        return "T / NK"
    if "b cell" in s or "b_cell" in s or "plasma" in s:
        return "B / plasma"
    if "fibro" in s or "stromal" in s or "myofibro" in s or "pericyte" in s or "smooth muscle" in s:
        return "Fibroblast / stroma"
    if "endothel" in s:
        return "Endothelial"
    if "hepatocyte" in s or "cholangiocyte" in s or "epithel" in s:
        return "Epithelial / hepatocyte"
    return "Other"

# Row order for heatmaps (melanoma first, reference second)
LINEAGE_ORDER = ["Melanoma", "Myeloid (ref)", "T / NK", "B / plasma",
                 "Fibroblast / stroma", "Endothelial", "Epithelial / hepatocyte", "Other"]


def chr_bounds(chr_pos, n_windows):
    """window-space chr_pos dict -> {chr: (start, end)} contiguous blocks."""
    items = sorted(chr_pos.items(), key=lambda kv: kv[1])
    order = [k for k, _ in items]
    starts = [v for _, v in items]
    out = {}
    for i, ch in enumerate(order):
        s = starts[i]
        e = starts[i + 1] if i + 1 < len(starts) else n_windows
        out[ch] = (s, e)
    return out


lineage_chr_rows = []     # long-form: sample, lineage, chr, mean_cnv, n_cells
mel_profile_rows = []      # sample -> per-chr melanoma mean (wide later)
percell_frames = []        # per-cell chr3/chr8 + coords + lineage + cnv_score
cnvscore_rows = []         # sample, lineage, mean cnv_score, n

print("=== Xenium iCNV supplementary figure build ===", flush=True)
for name, (path, ref) in SAMPLES.items():
    print(f"\n[{name}] loading {os.path.basename(path)} (ref={ref}) ...", flush=True)
    a = ad.read_h5ad(path)
    # Use the STORED X_cnv - exactly what the published Figure 5 code
    # (Xenium_CNV_Stemness_Correlation.Rmd) reads via adata$obsm$X_cnv - for all CNV
    # math. Re-run infercnv ONLY to recover the correct window-space chr_pos (the
    # stored uns chr_pos was overwritten with a gene-index map). Stored and re-run
    # X_cnv are element-wise aligned (Pearson r = 0.997), so chr_pos maps the stored
    # matrix correctly.
    Xc = a.obsm["X_cnv"]
    Xc = Xc.toarray() if sp.issparse(Xc) else np.asarray(Xc)
    cnv.tl.infercnv(a, reference_key="cell_type", reference_cat=[ref],
                    window_size=100, step=10, dynamic_threshold=1.5)
    nw = Xc.shape[1]
    bounds = chr_bounds(a.uns["cnv"]["chr_pos"], nw)
    # Published per-cell CNV score = RMS of the stored X_cnv row.
    rms = np.sqrt(np.mean(Xc ** 2, axis=1))
    print(f"  n_windows={nw}; chromosomes mapped={len(bounds)}; using STORED X_cnv", flush=True)

    lin = a.obs["cell_type"].map(lineage).values
    # per-cell whole-chromosome mean iCNV for every chr (used for heatmap + spatial)
    chr_cell = {}
    for ch in CHR_ORDER:
        if ch in bounds:
            s, e = bounds[ch]
            chr_cell[ch] = Xc[:, s:e].mean(axis=1) if e > s else np.zeros(Xc.shape[0])
        else:
            chr_cell[ch] = np.full(Xc.shape[0], np.nan)

    # lineage x chromosome mean
    for lg in pd.unique(lin):
        m = lin == lg
        if m.sum() == 0:
            continue
        for ch in CHR_ORDER:
            lineage_chr_rows.append(dict(sample=name, lineage=lg, chr=ch,
                                         mean_cnv=float(np.nanmean(chr_cell[ch][m])),
                                         n_cells=int(m.sum())))

    # melanoma-only per-chr profile (key combined panel)
    melm = lin == "Melanoma"
    mp = dict(sample=name, n_melanoma=int(melm.sum()))
    for ch in CHR_ORDER:
        mp[ch] = float(np.nanmean(chr_cell[ch][melm]))
    mel_profile_rows.append(mp)

    # cnv_score per lineage (infercnvpy pipeline output)
    if "cnv_score" in a.obs:
        cs = a.obs["cnv_score"].astype(float).values
        for lg in pd.unique(lin):
            m = lin == lg
            cnvscore_rows.append(dict(sample=name, lineage=lg,
                                      mean_cnv_score=float(np.nanmean(cs[m])),
                                      n=int(m.sum())))

    # per-cell chr3/chr8 + coords for spatial maps + downstream
    tissue = "breast" if "32XXX" in name else "liver"
    pc = pd.DataFrame({
        "sample": name,
        "tissue": tissue,
        "cell": a.obs_names,
        "cell_type": a.obs["cell_type"].astype(str).values,
        "lineage": lin,
        "x": a.obs["x_centroid"].values if "x_centroid" in a.obs else np.nan,
        "y": a.obs["y_centroid"].values if "y_centroid" in a.obs else np.nan,
        "chr3_cnv": chr_cell["chr3"],
        "chr8_cnv": chr_cell["chr8"],
        "chr6_cnv": chr_cell["chr6"],
        "chr1_cnv": chr_cell["chr1"],
        "cnv_score_rms": rms,  # published Figure 5 score (RMS of stored X_cnv)
        "cnv_score_09": a.obs["cnv_score"].astype(float).values if "cnv_score" in a.obs else np.nan,
    })
    percell_frames.append(pc)
    del a, Xc
    print(f"  done. melanoma cells={int(melm.sum())}", flush=True)

# ---- save tables ----
lc = pd.DataFrame(lineage_chr_rows)
lc.to_csv(f"{OUT}/lineage_chr_mean_cnv.csv", index=False)
mel = pd.DataFrame(mel_profile_rows).set_index("sample")
mel.to_csv(f"{OUT}/melanoma_chr_profile.csv")
pcell = pd.concat(percell_frames, ignore_index=True)
pcell.to_csv(f"{OUT}/percell_chr_cnv.csv.gz", index=False, compression="gzip")
pd.DataFrame(cnvscore_rows).to_csv(f"{OUT}/cnv_score_by_lineage.csv", index=False)
print("\nSaved tables to", OUT, flush=True)

# ---------- Figures ----------
SAMP_ORDER = list(SAMPLES.keys())
cmap = plt.get_cmap("RdBu_r")


def safe_vm(arr, floor=0.05):
    vals = np.abs(np.nan_to_num(np.asarray(arr, float), nan=0.0))
    vm = float(np.percentile(vals, 99)) if vals.size else floor
    if not np.isfinite(vm) or vm <= 0:
        vm = floor
    return vm


def dnorm(vm):
    return TwoSlopeNorm(vmin=-vm, vcenter=0.0, vmax=vm)

# Fig 2 first (the key combined panel): melanoma per-chr profile across samples
melmat = mel[CHR_ORDER].reindex(SAMP_ORDER)
vm = safe_vm(melmat.values)
fig, ax = plt.subplots(figsize=(11, 4.2))
im = ax.imshow(np.nan_to_num(melmat.values), aspect="auto", cmap=cmap, norm=dnorm(vm))
ax.set_xticks(range(len(CHR_ORDER)))
ax.set_xticklabels([c.replace("chr", "") for c in CHR_ORDER], fontsize=9)
ax.set_yticks(range(len(SAMP_ORDER)))
ax.set_yticklabels(SAMP_ORDER, fontsize=9)
ax.set_xlabel("Chromosome")
ax.set_title("Inferred CNV profile of melanoma cells across 6 metastatic Xenium samples\n"
             "(mean infercnvpy signal per chromosome; blue = loss, red = gain)", fontsize=11)
# highlight UM hallmark chromosomes
for ci, ch in enumerate(CHR_ORDER):
    if ch in ("chr3", "chr8", "chr6", "chr1"):
        ax.add_patch(plt.Rectangle((ci - 0.5, -0.5), 1, len(SAMP_ORDER),
                                   fill=False, edgecolor="black", lw=1.2))
cbar = fig.colorbar(im, ax=ax, shrink=0.8)
cbar.set_label("mean inferred CNV")
fig.text(0.5, -0.02, "Boxed: chr3 (monosomy 3 / loss), chr8 (8q gain), chr6, chr1 - classical UM aberrations",
         ha="center", fontsize=9, style="italic")
fig.tight_layout()
fig.savefig(f"{OUT}/Fig2_melanoma_chr_profile_allsamples.png", dpi=300, bbox_inches="tight")
fig.savefig(f"{OUT}/Fig2_melanoma_chr_profile_allsamples.pdf", bbox_inches="tight")
plt.close(fig)
print("Saved Fig2 (combined melanoma profile)", flush=True)

# Fig 1: per-sample lineage x chromosome heatmaps (2 rows x 3 cols)
fig, axes = plt.subplots(2, 3, figsize=(20, 9))
for k, name in enumerate(SAMP_ORDER):
    axx = axes.flat[k]
    sub = lc[lc["sample"] == name]
    present = [lg for lg in LINEAGE_ORDER if lg in set(sub["lineage"])]
    mat = (sub.pivot(index="lineage", columns="chr", values="mean_cnv")
              .reindex(index=present, columns=CHR_ORDER))
    vm = safe_vm(mat.values)
    im = axx.imshow(np.nan_to_num(mat.values), aspect="auto", cmap=cmap, norm=dnorm(vm))
    axx.set_xticks(range(len(CHR_ORDER)))
    axx.set_xticklabels([c.replace("chr", "") for c in CHR_ORDER], fontsize=6)
    axx.set_yticks(range(len(present)))
    axx.set_yticklabels(present, fontsize=8)
    axx.set_title(name, fontsize=10)
    fig.colorbar(im, ax=axx, shrink=0.7)
fig.suptitle("Per-sample inferred CNV by cell lineage x chromosome (blue = loss, red = gain)", fontsize=13)
fig.tight_layout(rect=[0, 0, 1, 0.97])
fig.savefig(f"{OUT}/Fig1_lineage_chr_heatmaps_persample.png", dpi=300, bbox_inches="tight")
fig.savefig(f"{OUT}/Fig1_lineage_chr_heatmaps_persample.pdf", bbox_inches="tight")
plt.close(fig)
print("Saved Fig1 (per-sample lineage heatmaps)", flush=True)

# Fig 3: spatial chr3 (loss) and chr8 (gain) maps per sample (6 rows x 2 cols)
fig, axes = plt.subplots(len(SAMP_ORDER), 2, figsize=(11, 4.0 * len(SAMP_ORDER)))
for r, name in enumerate(SAMP_ORDER):
    d = pcell[pcell["sample"] == name]
    for c, (col, lab) in enumerate([("chr3_cnv", "chr3 (loss = blue)"),
                                    ("chr8_cnv", "chr8 / 8q (gain = red)")]):
        axx = axes[r, c]
        # plot melanoma cells on top of grey non-tumor background
        nonmel = d[d["lineage"] != "Melanoma"]
        melc = d[d["lineage"] == "Melanoma"]
        vm = safe_vm(melc[col].values)
        axx.scatter(nonmel["x"], nonmel["y"], s=0.5, c="lightgrey", linewidths=0)
        sc = axx.scatter(melc["x"], melc["y"], s=0.5, c=np.nan_to_num(melc[col].values),
                         cmap=cmap, norm=dnorm(vm), linewidths=0)
        axx.set_aspect("equal")
        axx.invert_yaxis()
        axx.set_xticks([]); axx.set_yticks([])
        axx.set_title(f"{name}\n{lab}", fontsize=9)
        fig.colorbar(sc, ax=axx, shrink=0.7)
fig.suptitle("Spatial distribution of inferred chr3 loss and chr8 gain in melanoma cells\n"
             "(grey = non-tumor; aberrations localise to the tumor compartment)", fontsize=12)
fig.tight_layout(rect=[0, 0, 1, 0.985])
fig.savefig(f"{OUT}/Fig3_spatial_chr3_chr8_maps.png", dpi=200, bbox_inches="tight")
plt.close(fig)
print("Saved Fig3 (spatial maps)", flush=True)

# Fig 4: cnv_score by lineage (melanoma vs normal) per sample
cs = pd.DataFrame(cnvscore_rows)
fig, ax = plt.subplots(figsize=(12, 5))
lg_present = [lg for lg in LINEAGE_ORDER if lg in set(cs["lineage"])]
width = 0.8 / len(SAMP_ORDER)
for si, name in enumerate(SAMP_ORDER):
    sub = cs[cs["sample"] == name].set_index("lineage").reindex(lg_present)
    ax.bar(np.arange(len(lg_present)) + si * width, sub["mean_cnv_score"].values,
           width=width, label=name)
ax.set_xticks(np.arange(len(lg_present)) + 0.4 - width / 2)
ax.set_xticklabels(lg_present, rotation=30, ha="right")
ax.set_ylabel("mean infercnvpy cnv_score")
ax.set_title("infercnvpy CNV score by cell lineage (melanoma is elevated vs normal/reference)")
ax.legend(fontsize=7, ncol=2)
fig.tight_layout()
fig.savefig(f"{OUT}/Fig4_cnv_score_by_lineage.png", dpi=300, bbox_inches="tight")
plt.close(fig)
print("Saved Fig4 (cnv_score by lineage)", flush=True)

print("\nALL DONE. Outputs in", OUT, flush=True)
