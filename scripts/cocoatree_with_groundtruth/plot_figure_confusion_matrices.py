import matplotlib.pyplot as plt
import os
import itertools
import numpy as np
import pandas as pd
from utils.vis import sectors_cm
from utils.sectors import compute_IOU_metric
from plotmastery.utils_subfigure import add_letter_and_title
from plotmastery import utils_heatmap
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("dataset")
args = parser.parse_args()

dataset = args.dataset


def compute_all_vs_all(sectors_per_method, comp1=0, comp2=0):
    IOU_metric = np.zeros(
        (len(sectors_per_method), len(sectors_per_method)))
    for i, met1 in enumerate(sectors_per_method):
        for j, met2 in enumerate(sectors_per_method):
            IOU_metric[i, j] = compute_IOU_metric(met1[comp1], met2[comp2])
    return IOU_metric


# Load precomputed results
sca_res = pd.read_csv(
    f"results/cocoatree_gt/{dataset}/cocoatree_SCA_none.csv")
sca_res = sca_res.loc[~sca_res["filtered_msa_pos"].isna()]

mi_res = pd.read_csv(
    f"results/cocoatree_gt/{dataset}/cocoatree_MI_none.csv")
mi_res = mi_res.loc[~mi_res["filtered_msa_pos"].isna()]

nmi_res = pd.read_csv(
    f"results/cocoatree_gt/{dataset}/cocoatree_NMI_none.csv")
nmi_res = nmi_res.loc[~nmi_res["filtered_msa_pos"].isna()]

miapc_res = pd.read_csv(
    f"results/cocoatree_gt/{dataset}/cocoatree_MI_APC.csv")
miapc_res = miapc_res.loc[~miapc_res["filtered_msa_pos"].isna()]

sector_cols = [
    c for c in sca_res.columns if c.startswith("xcor")]
sector_cols.sort()

orig_sector_cols = [
    c for c in sca_res.columns if c.startswith("orig_sector")]

# Get original sectors in the original order
orig_sectors = [np.where(sca_res[c])[0] for c in orig_sector_cols]
sca_sectors = [np.where(sca_res[c])[0] for c in sector_cols]
mi_sectors = [np.where(mi_res[c])[0] for c in sector_cols]
nmi_sectors = [np.where(nmi_res[c])[0] for c in sector_cols]
miapc_sectors = [np.where(miapc_res[c])[0] for c in sector_cols]

all_sectors = [orig_sectors, sca_sectors, mi_sectors, nmi_sectors,
               miapc_sectors]

n_comp = len(sca_sectors)

letters = ["A", "B", "C", "D", "E", "F"]
sector_names = ["Red", "Green", "Blue", "Purple", "Crimson", "Orange"]

cmaps = list(sectors_cm.values())
fig, axes = plt.subplots(figsize=(8, 8), ncols=n_comp, nrows=n_comp)
for i, j in itertools.product(range(n_comp), range(n_comp)):
    if i != j:
        cmap = sectors_cm["others"]
    else:
        cmap = cmaps[i]

    ax = axes[i, j]
    if i > j:
        ax.set_axis_off()
        continue

    IOU_metric_comp1 = compute_all_vs_all(all_sectors, comp1=i, comp2=j)
    m = ax.matshow(IOU_metric_comp1, vmin=0, cmap=cmap)
    utils_heatmap.annotate_heatmap(
        m, valfmt="{x:1.0f}",
        fontsize="x-small")
    if i == 0:
        ax.set_xticklabels(
            ["", "ori.", "SCA", "MI", "NMI", "MI+APC"],
            fontsize="x-small")
        ax.set_title(sector_names[j], fontweight="bold", fontsize="medium")
    else:
        ax.set_xticks([])
    if i == j:
        ax.set_yticklabels(
            ["", "ori.", "SCA", "MI", "NMI", "MI+APC"],
            fontsize="x-small")
        ax.set_ylabel(sector_names[j], fontweight="bold")
    else:
        ax.set_yticks([])

os.makedirs(f"figures/{dataset}", exist_ok=True)
fig.savefig(f"figures/{dataset}/confusion_all.pdf")
fig.savefig(f"figures/{dataset}/confusion_all.png")

if dataset == "rivoire":
    fig, axes = plt.subplots(figsize=(16, 3), ncols=n_comp,
                             tight_layout=True)

else:
    fig, axes = plt.subplots(figsize=(8, 3), ncols=n_comp,
                             tight_layout=True)
for i in range(n_comp):

    ax = axes[i]
    IOU_metric_comp1 = compute_all_vs_all(
        all_sectors, comp1=i,
        comp2=i)
    m = ax.matshow(IOU_metric_comp1, vmin=0, cmap=cmaps[i])
    utils_heatmap.annotate_heatmap(
        m, valfmt="{x:1.0f}",
        fontsize="x-small")
    ax.set_xticklabels(
            ["", "ori.", "SCA", "MI", "NMI", "MI+APC"],
            fontsize="x-small")
    ax.set_yticklabels(
            ["", "ori.", "SCA", "MI", "NMI", "MI+APC"],
            fontsize="x-small")

    ax.tick_params(axis='both', which='both', labelsize='x-small',
                   bottom=True, top=False, labeltop=False, labelbottom=True,
                   left=True, right=False, labelright=False, labelleft=True)

    add_letter_and_title(ax, letters[i], sector_names[i] + " XCoR",
                         fontsize="small")


fig.savefig(f"figures/{dataset}/confusion_row.pdf")
fig.savefig(f"figures/{dataset}/confusion_row.png")
