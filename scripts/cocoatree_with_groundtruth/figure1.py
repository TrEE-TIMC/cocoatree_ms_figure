import matplotlib.pyplot as plt
import os
import numpy as np
import pandas as pd
from utils.vis import plot_scatter_sectors, create_legend
from utils.vis import plot_sectors
from matplotlib.gridspec import GridSpec
from plotmastery.utils_subfigure import add_letter_and_title
from utils.postprocessing import annotate_results


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
dataset = "rhomboid"
results = pd.read_csv(
    f"results/cocoatree_gt/{dataset}/cocoatree_SCA_none.csv")
results = annotate_results(results)
results = results.loc[~results["pdb_pos"].isna()]

# Cocoatree sectors
ax = fig.add_subplot(gs[17, :-2])
add_letter_and_title(ax, "C.", title=dataset.capitalize())
if dataset == "halabi":
    order = ["sector_3", "sector_1", "sector_2"]
else:
    order = ["sector_2", "sector_3", "sector_1"]
plot_sectors(ax, results, columns=order, title="cocoatree")

# Original sectors
ax = fig.add_subplot(gs[18, :-2])
plot_sectors(ax, results,
             columns=["orig_sector_1", "orig_sector_2", "orig_sector_3"],
             title="orig")

# Both
ax = fig.add_subplot(gs[19, :-2])
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

os.makedirs("figures", exist_ok=True)
fig.savefig("figures/figure_1.pdf")
