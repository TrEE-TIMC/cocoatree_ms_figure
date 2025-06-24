import matplotlib.pyplot as plt
import itertools
import numpy as np
import pandas as pd
from utils.vis import sectors_cm
from utils.sectors import get_best_ordered_sectors
from plotmastery.utils_subfigure import add_letter_and_title
from plotmastery import utils_heatmap
from utils.sectors import compute_IOU_metric

datasets = ["halabi", "rhomboid", "DHFR"]

def compute_all_vs_all(sectors_per_method, comp1=0, comp2=0):
    IOU_metric = np.zeros(
        (len(sectors_per_method), len(sectors_per_method)))
    for i, met1 in enumerate(sectors_per_method):
        for j, met2 in enumerate(sectors_per_method):
            IOU_metric[i, j] = compute_IOU_metric(met1[comp1], met2[comp2])
    return IOU_metric


cmaps = list(sectors_cm.values())

fig  = plt.figure(figsize=(7.5, 7.5), 
                  tight_layout=True)
axes = []
for j, dataset in enumerate(datasets):
    dataset_ax = []
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
    sector_cols.sort()
    orig_sector_cols = [c for c in sca_res.columns if c.startswith("orig_sector")]

    # Get original sectors in the original order
    orig_sectors = [np.where(sca_res[c])[0] for c in orig_sector_cols]
    sca_sectors = [np.where(sca_res[c])[0] for c in sector_cols]
    mi_sectors = [np.where(mi_res[c])[0] for c in sector_cols]
    nmi_sectors = [np.where(nmi_res[c])[0] for c in sector_cols]
    miapc_sectors = [np.where(miapc_res[c])[0] for c in sector_cols]


    all_sectors = [orig_sectors, sca_sectors, mi_sectors, nmi_sectors,
                miapc_sectors]

    n_comp = len(sca_sectors)

    letters = ["A", "B", "C", "D"]
    sector_names = ["Green", "Red", "Blue", "Purple"]


    for i in range(n_comp):

        ax = fig.add_subplot(len(datasets), n_comp, (j*n_comp+i+1))
        dataset_ax.append(ax)
        IOU_metric_comp1 = compute_all_vs_all(all_sectors, comp1=i, comp2=i)
        m = ax.matshow(IOU_metric_comp1, vmin=0, cmap=cmaps[i])
        utils_heatmap.annotate_heatmap(
            m, valfmt="{x:1.0f}",
            fontsize="x-small")
        if i == 0:
            ax.set_yticklabels(
                    ["", "ori.", "SCA", "MI", "NMI", "MI+APC"],
                    fontsize="x-small")
        else:
            ax.set_yticklabels([])

        if j == 2:
            ax.set_xticklabels(
                    ["", "ori.", "SCA", "MI", "NMI", "MI+APC"],
                    fontsize="x-small", rotation=90)
        else:
            ax.set_xticklabels([])
        ax.tick_params(axis='both', which='both', labelsize='x-small',
                bottom=True, top=False, labeltop=False, labelbottom=True,
                left=True, right=False, labelright=False, labelleft=True)
    axes.append(dataset_ax)

add_letter_and_title(axes[0][0], letters[0], "Serine protease")
add_letter_and_title(axes[1][0], letters[1], "Rhomboid")
add_letter_and_title(axes[2][0], letters[2], "DHFR")

fig.savefig(f"figures/figure_4.pdf")
