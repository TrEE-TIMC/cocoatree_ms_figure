import os
import matplotlib.pyplot as plt
import pandas as pd
import argparse
from plotmastery.utils_subfigure import add_letter_and_title

from utils.postprocessing import annotate_results
from utils.vis import plot_scatter_sectors

parser = argparse.ArgumentParser()
parser.add_argument("filename")
parser.add_argument("--outname", "-o")
parser.add_argument("--pca", action="store_true", default=False)
args = parser.parse_args()

filename = args.filename
outname = args.outname
pca = args.pca


results = pd.read_csv(filename)
results = results.loc[~results["filtered_msa_pos"].isna()]
results = annotate_results(results)

sca_sectors = [c for c in results.columns if c.startswith("xcor")]
n_comp = len(sca_sectors)

fig, axes = plt.subplots(
    nrows=(n_comp-1), ncols=(n_comp-1),
    figsize=(4, 4), tight_layout=True)


coevolution_metric = filename.split("_")[-2]
correction = filename.split("_")[-1].split(".")[0]

# Manually add letters
if coevolution_metric == "SCA":
    letter = "A." if pca else "B."
elif coevolution_metric == "MI" and correction == "none":
    letter = "D." if pca else "E."
elif coevolution_metric == "NMI":
    letter = "G." if pca else "H."
elif coevolution_metric == "MI" and correction == "APC":
    letter = "J." if pca else "K."
else:
    letter = None


if pca:
    title = f"PCA with {coevolution_metric}"
else:
    title = f"ICA with {coevolution_metric}"
if correction != "none":
    title = title + f" ({correction})"

colors = ["Red", "Green", "Blue", "Purple", "Crimson", "Orange"]
for i in range(n_comp):
    for j in range(n_comp):
        if i <= j:
            if i != 0 and j != n_comp - 1:
                fig.delaxes(axes[j, i-1])
            continue
        ax = axes[j, i-1]
        if j == 0 and i-1 == 0 and letter is not None:
            add_letter_and_title(ax, letter, title)
        if pca:
            plot_scatter_sectors(
                ax, results, f"PC{i+1}", f"PC{j+1}",
                annotate=False, add_labels=True)
        else:
            plot_scatter_sectors(
                ax, results, f"IC{i+1}", f"IC{j+1}",
                annotate=False, add_labels=False)
            ax.set_xlabel(colors[i], fontweight="bold", fontsize="small",
                          labelpad=2)
            ax.set_ylabel(colors[j], fontweight="bold", fontsize="small",
                          labelpad=2)

        # Hide the top and right spines
        ax.spines['top'].set_color('none')
        ax.spines['right'].set_color('none')


if outname is not None:
    os.makedirs(os.path.dirname(outname), exist_ok=True)
    fig.savefig(outname, dpi=300)
    fig.savefig(outname.replace(".png", ".pdf"))
