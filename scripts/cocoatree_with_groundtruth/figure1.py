import matplotlib.pyplot as plt
import os
import numpy as np
import pandas as pd
from utils.vis import plot_scatter_sectors, create_legend
from matplotlib.gridspec import GridSpec
from plotmastery.utils_subfigure import add_letter_and_title
from utils.colors import sectors_cm


def annotate_results(results):
    results = results.loc[~results["filtered_msa_pos"].isna()]

    sca_sectors = [c for c in results.columns if c.startswith("sector")]
    orig_sectors = [c for c in results.columns if c.startswith("orig_")]

    results["is_cocoatree_sector"] = results[sca_sectors].sum(
        axis=1) > 0
    results["is_orig_sector"] = results[orig_sectors].sum(
        axis=1) > 0
    results["is_both"] = (
        results["is_cocoatree_sector"] & results["is_orig_sector"])
    results["is_only_cocoatree"] = (
        results["is_cocoatree_sector"] & ~results["is_orig_sector"])
    results["is_only_orig"] = (
        ~results["is_cocoatree_sector"] & results["is_orig_sector"])
    return results


fig = plt.figure(figsize=(7.5, 4))
gs = GridSpec(20, 20, figure=fig, top=0.9, left=0.1)

###############################################################################
# Start with halabi results
results = pd.read_csv("results/cocoatree_gt/halabi/cocoatree_SCA_none.csv")
results = annotate_results(results)

ax = fig.add_subplot(gs[:12, :8])
plot_scatter_sectors(ax, results, "IC1", "IC2", annotate=False)
add_letter_and_title(ax, "A.", "Halabi")

###############################################################################
# rhomboid results
results = pd.read_csv("results/cocoatree_gt/rhomboid/cocoatree_SCA_none.csv")
results = annotate_results(results)

ax = fig.add_subplot(gs[:12, 10:-2])
plot_scatter_sectors(ax, results, "IC1", "IC2", annotate=False)
add_letter_and_title(ax, "B.", "Rhomboid")

###############################################################################
# Legend

legend = create_legend()
sec_legend = ax.legend(
    loc=(1, .7), frameon=False, fontsize="x-small",
    handles=legend["sectors"], title="Sectors",
    title_fontproperties={"weight": "bold",
                          "size": "small"},
    alignment="left")
ax.legend(loc=(1, .15), frameon=False, fontsize="x-small",
          handles=legend["methods"], title="Methods",
          title_fontproperties={"weight": "bold",
                                "size": "small"},
          alignment="left")
ax.add_artist(sec_legend)

###############################################################################
# Mapping between all
# results = pd.read_csv("results/cocoatree_gt/halabi/cocoatree_SCA_none.csv")
# results = annotate_results(results)
results = results.loc[~results["pdb_pos"].isna()]

# Cocoatree sectors
ax = fig.add_subplot(gs[17, :-2])
add_letter_and_title(ax, "C.", title="cocoatree vs pysca on rhomboid")
sec1 = results["sector_1"].values[np.newaxis, :].astype(float)
sec1[sec1 == 0] = np.nan
ax.matshow(sec1, aspect="auto", cmap=sectors_cm["sector_1"], vmin=0)

sec2 = results["sector_2"].values[np.newaxis, :].astype(float)
sec2[sec2 == 0] = np.nan
ax.matshow(sec2, aspect="auto", cmap=sectors_cm["sector_2"], vmin=0)

sec3 = results["sector_3"].values[np.newaxis, :].astype(float)
sec3[sec3 == 0] = np.nan
ax.matshow(sec3, aspect="auto", cmap=sectors_cm["sector_3"], vmin=0)

ax.set_xticks([])
ax.set_yticks([])
ax.set_yticks([0])
ax.set_yticklabels(["cocoatree"], fontweight="bold")
ax.tick_params(axis='both', which='both', labelsize='x-small',
               bottom=False, top=False, labeltop=False, labelbottom=False,
               left=False, right=False, labelleft=False, labelright=True)

# Original sectors
ax = fig.add_subplot(gs[18, :-2])
sec1 = results["orig_sector_1"].values[np.newaxis, :].astype(float)
sec1[sec1 == 0] = np.nan
ax.matshow(sec1, aspect="auto", cmap=sectors_cm["sector_3"], vmin=0)

sec2 = results["orig_sector_2"].values[np.newaxis, :].astype(float)
sec2[sec2 == 0] = np.nan
ax.matshow(sec2, aspect="auto", cmap=sectors_cm["sector_1"], vmin=0)

sec3 = results["orig_sector_0"].values[np.newaxis, :].astype(float)
sec3[sec3 == 0] = np.nan
ax.matshow(sec3, aspect="auto", cmap=sectors_cm["sector_2"], vmin=0)
ax.set_yticks([0])
ax.set_yticklabels(["original"], fontweight="bold")
ax.tick_params(axis='both', which='both', labelsize='x-small',
               bottom=False, top=False, labeltop=False, labelbottom=False,
               left=False, right=False, labelleft=False, labelright=True)

# Both
ax = fig.add_subplot(gs[19, :-2])
diff = results["is_both"].values[np.newaxis, :].astype(float)
ax.matshow(diff, aspect="auto", cmap="Greys")
ax.set_yticks([0])
ax.set_yticklabels(["both"], fontweight="bold")
ax.tick_params(axis='both', which='both', labelsize='x-small',
               bottom=True, top=False, labeltop=False, labelbottom=True,
               left=False, right=False, labelleft=False, labelright=True)
ax.set_xticks([0, len(results)])
os.makedirs("figures", exist_ok=True)
fig.savefig("figures/figure_1.pdf")
