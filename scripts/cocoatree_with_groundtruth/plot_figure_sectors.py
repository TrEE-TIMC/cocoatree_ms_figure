import matplotlib.pyplot as plt
import os
import numpy as np
import pandas as pd
from utils.vis import plot_sectors
from matplotlib.gridspec import GridSpec
from plotmastery.utils_subfigure import add_letter_and_title
from utils.postprocessing import annotate_results
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("dataset")
args = parser.parse_args()

dataset = args.dataset


fig = plt.figure(figsize=(7.5, 6))
gs = GridSpec(40, 40, figure=fig, top=0.9, left=0.1)

extensions = ["cocoatree_SCA_none.csv",
              "cocoatree_MI_none.csv",
              "cocoatree_NMI_none.csv",
              "cocoatree_MI_APC.csv"]

start_i = 0
titles = ["SCA", "MI", "NMI", "MI+APC"]
letters = ["A.", "B.", "C.", "D."]

for i, ext in enumerate(extensions):
    results = pd.read_csv(
        f"results/cocoatree_gt/{dataset}/{ext}")
    results = annotate_results(results)
    results = results.loc[~results["pdb_pos"].isna()]

    columns = [col for col in results.columns if col.startswith("orig_sector")]
    order = [col for col in results.columns if
             col.startswith("sector")]
    order.sort()
    
    # Cocoatree sectors
    ax = fig.add_subplot(gs[start_i, :-2])
    add_letter_and_title(ax, letters[i], title=titles[i])
    plot_sectors(ax, results, columns=order, title="cocoatree")

    # Original sectors
    ax = fig.add_subplot(gs[start_i+1, :-2])
    plot_sectors(ax, results,
                 columns=columns,
                 title="orig")

    # Both
    ax = fig.add_subplot(gs[start_i+2, :-2])
    diff = results["is_both"].values[np.newaxis, :].astype(float)
    ax.matshow(diff, aspect="auto", cmap="Greys")
    ax.set_yticks([0])
    ax.set_yticklabels(["both"], fontweight="bold")
    ax.tick_params(axis='both', which='both', labelsize='x-small',
                   bottom=True, top=False, labeltop=False, labelbottom=True,
                   left=False, right=False, labelleft=False, labelright=True)
    ax.xaxis.set_major_locator(plt.MaxNLocator(5))
    ax.set_xlabel("Position in PDB", fontweight="bold", fontsize="small",
                  labelpad=2)

    start_i += 8

ax = fig.add_subplot(gs[start_i:, :-2])

ax1 = ax.twinx()
ax.plot(np.arange(results.shape[0]), results["msa_conservation"],
        color="black", label="conservation")
ax1.plot(np.arange(results.shape[0]), results["msa_perc_gap"], color="gray",
        label="Percentage gaps")
ax.spines["top"].set_linewidth(0)
ax.spines["right"].set_linewidth(0)
ax1.spines["top"].set_linewidth(0)

ax.tick_params(axis='both', which='both', labelsize='x-small',
               bottom=True, top=False, labeltop=False, labelbottom=True,
               left=True, right=False, labelleft=True, labelright=False)

ax1.tick_params(axis='both', which='both', labelsize='x-small',
                bottom=True, top=False, labeltop=False, labelbottom=True,
                left=False, right=False, labelleft=False, labelright=True,
                labelcolor="gray")
ax1.set_ylabel("% gaps", fontweight="bold", fontsize="small", color="gray")
ax.set_ylabel("Conservation", fontweight="bold", fontsize="small")


ax.xaxis.set_major_locator(plt.MaxNLocator(5))

ax.set_xlabel("Position in PDB", fontweight="bold", fontsize="small",
                labelpad=2)
ax.set_xlim((0, len(results)))

add_letter_and_title(ax,"E.", title="Conservation and gaps")

os.makedirs(f"figures/{dataset}", exist_ok=True)
fig.savefig(f"figures/{dataset}/sectors_pdb.png")
fig.savefig(f"figures/{dataset}/sectors_pdb.pdf")
