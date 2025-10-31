import matplotlib.pyplot as plt
import os
import numpy as np
import pandas as pd
from utils.vis import plot_scatter_sectors, create_legend
from utils.vis import plot_sectors
from matplotlib.gridspec import GridSpec
from plotmastery.utils_subfigure import add_letter_and_title
from utils.postprocessing import annotate_results


fig = plt.figure(figsize=(7.5, 4.2))
gs = GridSpec(38, 30, figure=fig, top=0.94, left=0.07, right=0.93, bottom=0.08)


###############################################################################
# Start with halabi results
results = pd.read_csv("results/cocoatree_gt/halabi/cocoatree_SCA_none.csv")
results = annotate_results(results)

ax = fig.add_subplot(gs[:12, :8])
plot_scatter_sectors(ax, results, "IC1", "IC2", annotate=False,
                     add_labels=False)
ax.set_ylabel("Green", fontweight="bold", fontsize="small", labelpad=2)
ax.set_xlabel("Red", fontweight="bold", fontsize="small", labelpad=2)

add_letter_and_title(ax, "A.", "Serine protease")


###############################################################################
# DHFR results
results = pd.read_csv("results/cocoatree_gt/DHFR/cocoatree_SCA_none.csv")
results = annotate_results(results)

ax = fig.add_subplot(gs[:12, 20:-2])
plot_scatter_sectors(
    ax, results, "IC1", "IC2", annotate=False,
    add_labels=False)
ax.set_ylabel("Green", fontweight="bold", fontsize="small", labelpad=2)
ax.set_xlabel("Red", fontweight="bold", fontsize="small", labelpad=2)

add_letter_and_title(ax, "B.", "DHFR")


###############################################################################
# rhomboid results
results = pd.read_csv("results/cocoatree_gt/rhomboid/cocoatree_SCA_none.csv")
results = annotate_results(results)

ax = fig.add_subplot(gs[:12, 10:18])
plot_scatter_sectors(ax, results, "IC1", "IC2", annotate=False,
                     add_labels=False)
ax.set_ylabel("Green", fontweight="bold", fontsize="small", labelpad=2)
ax.set_xlabel("Red", fontweight="bold", fontsize="small", labelpad=2)

add_letter_and_title(ax, "C.", "Rhomboid")


###############################################################################
# Legend

legend = create_legend()
sec_legend = ax.legend(
    loc=(1.05, .5), frameon=False, fontsize="x-small",
    handles=legend["sectors"], title="Coco group",
    title_fontproperties={"weight": "bold",
                          "size": "small"},
    alignment="left")
ax.legend(loc=(1.05, -0.3), frameon=False, fontsize="x-small",
          handles=legend["methods"], title="Methods",
          title_fontproperties={"weight": "bold",
                                "size": "small"},
          alignment="left")
ax.add_artist(sec_legend)

###############################################################################
# Mapping between all
shift = 8
start_i = 19
dataset = "halabi"
results = pd.read_csv(
    f"results/cocoatree_gt/{dataset}/cocoatree_SCA_none.csv")
results = annotate_results(results)
results = results.loc[~results["pdb_pos"].isna()]

columns = [col for col in results.columns if col.startswith("orig_sector")]
order = [col for col in results.columns if col.startswith("sector")]
# Cocoatree sectors
ax = fig.add_subplot(gs[start_i, :-2])
add_letter_and_title(ax, "D.", title="Serine protease")
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

start_i += shift

dataset = "DHFR"
results = pd.read_csv(
    f"results/cocoatree_gt/{dataset}/cocoatree_SCA_none.csv")
results = annotate_results(results)
results = results.loc[~results["pdb_pos"].isna()]

columns = [col for col in results.columns if col.startswith("orig_sector")]
order = [col for col in results.columns if col.startswith("sector")]
# Cocoatree sectors
ax = fig.add_subplot(gs[start_i, :-2])
add_letter_and_title(ax, "E.", title="DHFR")
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

start_i += shift

dataset = "rhomboid"
results = pd.read_csv(
    f"results/cocoatree_gt/{dataset}/cocoatree_SCA_none.csv")
results = annotate_results(results)
results = results.loc[~results["pdb_pos"].isna()]

columns = [col for col in results.columns if col.startswith("orig_sector")]
order = [col for col in results.columns if col.startswith("sector")]
# Cocoatree sectors
ax = fig.add_subplot(gs[start_i, :-2])
add_letter_and_title(ax, "F.", title="Rhomboid")
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


os.makedirs("figures", exist_ok=True)
fig.savefig("figures/figure_2.pdf")
fig.savefig("figures/figure_2.png", dpi=300)
