import matplotlib.pyplot as plt
import itertools
import numpy as np
import pandas as pd
from plotmastery.utils_subfigure import add_letter_and_title
from plotmastery import utils_heatmap

dataset = "rhomboid"

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


def get_best_ordered_sectors(res):
    sector_cols = [c for c in res.columns if c.startswith("sector")]

    sectors = []
    # It's not the best strategy, but just get iteratively the "best
    # sector"
    for s in orig_sectors:
        all_scores = []
        for s1 in sector_cols:
            all_scores.append(
                compute_IOU_metric(
                    s,
                    np.where(res[s1])[0]))

        s = sector_cols.pop(np.argmax(all_scores))
        sectors.append(np.where(res[s])[0])
    return sectors


sca_sectors = get_best_ordered_sectors(sca_res)
mi_sectors = get_best_ordered_sectors(mi_res)
nmi_sectors = get_best_ordered_sectors(nmi_res)
miapc_sectors = get_best_ordered_sectors(miapc_res)
all_sectors = [orig_sectors, sca_sectors, mi_sectors, nmi_sectors,
               miapc_sectors]

n_comp = len(sca_sectors)

letters = ["A", "B", "C"]
sector_names = ["Red", "Blue", "Green"]

fig, axes = plt.subplots(figsize=(8, 8), ncols=n_comp, nrows=n_comp)
for i, j in itertools.product(range(3), range(3)):

    ax = axes[i, j]
    if i > j:
        ax.set_axis_off()
        continue

    IOU_metric_comp1 = compute_all_vs_all(all_sectors, comp1=i, comp2=j)
    m = ax.matshow(IOU_metric_comp1, vmin=0, cmap="Oranges")
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


cmaps = ["Reds", "Blues", "Greens"]
fig, axes = plt.subplots(figsize=(8, 3), ncols=n_comp,
                         tight_layout=True)
for i in range(n_comp):

    ax = axes[i]
    IOU_metric_comp1 = compute_all_vs_all(all_sectors, comp1=i,
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

    add_letter_and_title(ax, letters[i], sector_names[i])
