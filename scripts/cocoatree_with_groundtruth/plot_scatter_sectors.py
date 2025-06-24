import os
import matplotlib.pyplot as plt
import pandas as pd
import argparse
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

sca_sectors = [c for c in results.columns if c.startswith("sector")]
n_comp = len(sca_sectors)

fig, axes = plt.subplots(
    nrows=(n_comp-1), ncols=(n_comp-1),
    figsize=(4, 4))


coevolution_metric = filename.split("_")[-2]
correction = filename.split("_")[-1].split(".")[0]

for i in range(n_comp):
    for j in range(n_comp):
        if i <= j:
            if i != 0 and j != n_comp - 1:
                fig.delaxes(axes[j, i-1])
            continue
        ax = axes[j, i-1]
        if pca:
            plot_scatter_sectors(
                ax, results, f"PC{i+1}", f"PC{j+1}",
                annotate=False)
        else:
            plot_scatter_sectors(
                ax, results, f"IC{i+1}", f"IC{j+1}",
                annotate=False)
        ax.legend(frameon=False)
        # Move the left and bottom spines to the center
        ax.spines['left'].set_position('zero')
        ax.spines['bottom'].set_position('zero')

        # Hide the top and right spines
        ax.spines['top'].set_color('none')
        ax.spines['right'].set_color('none')

title = coevolution_metric
if correction != None:
    title = title + f" ({correction})"
fig.suptitle(coevolution_metric, fontweight="bold", fontsize="small")

if outname is not None:
    os.makedirs(os.path.dirname(outname), exist_ok=True)
    fig.savefig(outname)
    fig.savefig(outname.replace(".png", ".pdf"))
