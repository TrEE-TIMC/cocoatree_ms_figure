import matplotlib.pyplot as plt
import os
import itertools
import numpy as np
import pandas as pd
from utils.vis import sectors_cm
from plotmastery.utils_subfigure import add_letter_and_title
from plotmastery import utils_heatmap
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("dataset")
args = parser.parse_args()

dataset = args.dataset

# Load precomputed results
sca_res = pd.read_csv(f"results/cocoatree_gt/{dataset}/cocoatree_SCA_none.csv")
sca_res = sca_res.loc[~sca_res["filtered_msa_pos"].isna()]

mi_res = pd.read_csv(f"results/cocoatree_gt/{dataset}/cocoatree_MI_none.csv")
mi_res = mi_res.loc[~mi_res["filtered_msa_pos"].isna()]

nmi_res = pd.read_csv(f"results/cocoatree_gt/{dataset}/cocoatree_NMI_none.csv")
nmi_res = nmi_res.loc[~nmi_res["filtered_msa_pos"].isna()]

miapc_res = pd.read_csv(f"results/cocoatree_gt/{dataset}/cocoatree_MI_APC.csv")
miapc_res = miapc_res.loc[~miapc_res["filtered_msa_pos"].isna()]


sector_cols = [c for c in sca_res.columns if c.startswith("sector")]
orig_sector_cols = [c for c in sca_res.columns if c.startswith("orig_sector")]

# Get original sectors in the original order
orig_sectors = [np.where(sca_res[c])[0] for c in orig_sector_cols]


# Now create some form of confusion matrix
def compute_IOU_metric(set1, set2, metric="intersection"):
    union = len(set(set1).union(set2))
    intersection = len(set(set1).intersection(set2))
    if metric == "IOU":
        return intersection / union
    else:
        return intersection


def compute_all_vs_all(sectors_per_method, comp1=0, comp2=0):
    IOU_metric = np.zeros(
        (len(sectors_per_method), len(sectors_per_method)))
    for i, met1 in enumerate(sectors_per_method):
        for j, met2 in enumerate(sectors_per_method):
            IOU_metric[i, j] = compute_IOU_metric(met1[comp1], met2[comp2])
    return IOU_metric


def get_best_ordered_sectors(res, dataset="halabi", type="SCA"):
    sector_cols = [c for c in res.columns if c.startswith("sector")]

    sectors = []
    # It's not the best strategy, but just get iteratively the "best
    # sector"

    all_scores = np.zeros((len(sector_cols), len(sector_cols)))
    for i, s in enumerate(orig_sectors):
        for j, s1 in enumerate(sector_cols):
            all_scores[i, j] = compute_IOU_metric(
                    s,
                    np.where(res[s1])[0])
    order = all_scores.argmax(axis=1)
    if len(np.unique(order)) != len(sector_cols):
        order = np.arange(len(sector_cols))
        if dataset == "halabi" and type == "NMI":
            order = [2, 1, 0]
        if dataset == "halabi" and type == "MI":
            order = [1, 0, 2]
        if dataset == "rhomboid" and type == "MI":
            order = [0, 2, 1]

    sectors = [np.where(res[sector_cols[o]])[0] for o in order]
    return sectors


sca_sectors = get_best_ordered_sectors(sca_res, dataset=dataset, type="SCA")
mi_sectors = get_best_ordered_sectors(mi_res, dataset=dataset, type="MI")
nmi_sectors = get_best_ordered_sectors(nmi_res, dataset=dataset, type="NMI")
miapc_sectors = get_best_ordered_sectors(
    miapc_res, dataset=dataset,
    type="MIAPC")
all_sectors = [orig_sectors, sca_sectors, mi_sectors, nmi_sectors,
               miapc_sectors]

n_comp = len(sca_sectors)

letters = ["A", "B", "C", "D", "E", "F"]
sector_names = ["Green", "Red", "Blue", "Purple", "Crimson", "Orange"]

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

    add_letter_and_title(ax, letters[i], sector_names[i] + " sector",
                         fontsize="small")


fig.savefig(f"figures/{dataset}/confusion_row.pdf")
fig.savefig(f"figures/{dataset}/confusion_row.png")
