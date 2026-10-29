#!/usr/bin/env python3
"""Render the Xenium iCNV figures from the cached CSVs produced by
xenium_icnv_supp_figure.py (re-render step; no infercnv recompute).
Robust to all-zero / degenerate colour ranges."""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm

OUT = "../Xenium_iCNV_SuppFig"
CHR_ORDER = [f"chr{i}" for i in range(1, 23)]
LINEAGE_ORDER = ["Melanoma", "Myeloid (ref)", "T / NK", "B / plasma",
                 "Fibroblast / stroma", "Endothelial", "Epithelial / hepatocyte", "Other"]
cmap = plt.get_cmap("RdBu_r")


def safe_vm(arr, floor=0.05):
    vals = np.abs(np.nan_to_num(np.asarray(arr, float), nan=0.0))
    vm = float(np.percentile(vals, 99)) if vals.size else floor
    if not np.isfinite(vm) or vm <= 0:
        vm = floor
    return vm


def dnorm(vm):
    return TwoSlopeNorm(vmin=-vm, vcenter=0.0, vmax=vm)


lc = pd.read_csv(f"{OUT}/lineage_chr_mean_cnv.csv")
mel = pd.read_csv(f"{OUT}/melanoma_chr_profile.csv").set_index("sample")
cs = pd.read_csv(f"{OUT}/cnv_score_by_lineage.csv")
pcell = pd.read_csv(f"{OUT}/percell_chr_cnv.csv.gz")
SAMP_ORDER = list(mel.index)

# ---- Fig 2: combined melanoma chr profile across samples ----
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
for ci, ch in enumerate(CHR_ORDER):
    if ch in ("chr3", "chr8", "chr6", "chr1"):
        ax.add_patch(plt.Rectangle((ci - 0.5, -0.5), 1, len(SAMP_ORDER),
                                   fill=False, edgecolor="black", lw=1.2))
fig.colorbar(im, ax=ax, shrink=0.8).set_label("mean inferred CNV")
fig.text(0.5, -0.02, "Boxed: chr3 (monosomy 3 / loss), chr8 (8q gain), chr6, chr1 - classical UM aberrations",
         ha="center", fontsize=9, style="italic")
fig.tight_layout()
fig.savefig(f"{OUT}/Fig2_melanoma_chr_profile_allsamples.png", dpi=300, bbox_inches="tight")
fig.savefig(f"{OUT}/Fig2_melanoma_chr_profile_allsamples.pdf", bbox_inches="tight")
plt.close(fig)
print("Fig2 OK")

# ---- Fig 1: per-sample lineage x chromosome ----
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
print("Fig1 OK")

# ---- Fig 3: spatial chr3 / chr8 maps ----
fig, axes = plt.subplots(len(SAMP_ORDER), 2, figsize=(11, 4.0 * len(SAMP_ORDER)))
for r, name in enumerate(SAMP_ORDER):
    d = pcell[pcell["sample"] == name]
    for c, (col, lab) in enumerate([("chr3_cnv", "chr3 (loss = blue)"),
                                    ("chr8_cnv", "chr8 / 8q (gain = red)")]):
        axx = axes[r, c]
        melc = d[d["lineage"] == "Melanoma"]
        nonmel = d[d["lineage"] != "Melanoma"]
        vm = safe_vm(melc[col].values)
        axx.scatter(nonmel["x"], nonmel["y"], s=0.5, c="lightgrey", linewidths=0)
        sc = axx.scatter(melc["x"], melc["y"], s=0.5, c=np.nan_to_num(melc[col].values),
                         cmap=cmap, norm=dnorm(vm), linewidths=0)
        axx.set_aspect("equal"); axx.invert_yaxis()
        axx.set_xticks([]); axx.set_yticks([])
        axx.set_title(f"{name}\n{lab}", fontsize=9)
        fig.colorbar(sc, ax=axx, shrink=0.7)
fig.suptitle("Spatial distribution of inferred chr3 loss and chr8 gain in melanoma cells\n"
             "(grey = non-tumor; aberrations localise to the tumor compartment)", fontsize=12)
fig.tight_layout(rect=[0, 0, 1, 0.985])
fig.savefig(f"{OUT}/Fig3_spatial_chr3_chr8_maps.png", dpi=200, bbox_inches="tight")
plt.close(fig)
print("Fig3 OK")

# ---- Fig 4: cnv_score by lineage ----
fig, ax = plt.subplots(figsize=(12, 5))
lg_present = [lg for lg in LINEAGE_ORDER if lg in set(cs["lineage"])]
width = 0.8 / len(SAMP_ORDER)
for si, name in enumerate(SAMP_ORDER):
    sub = cs[cs["sample"] == name].set_index("lineage").reindex(lg_present)
    ax.bar(np.arange(len(lg_present)) + si * width, np.nan_to_num(sub["mean_cnv_score"].values),
           width=width, label=name)
ax.set_xticks(np.arange(len(lg_present)) + 0.4 - width / 2)
ax.set_xticklabels(lg_present, rotation=30, ha="right")
ax.set_ylabel("mean infercnvpy cnv_score")
ax.set_title("infercnvpy CNV score by cell lineage (melanoma is elevated vs normal/reference)")
ax.legend(fontsize=7, ncol=2)
fig.tight_layout()
fig.savefig(f"{OUT}/Fig4_cnv_score_by_lineage.png", dpi=300, bbox_inches="tight")
plt.close(fig)
print("Fig4 OK")
print("ALL FIGURES DONE")
