import matplotlib.pyplot as plt
import os
import pandas as pd
from utils.vis import plot_coev_mat_sectors
from plotmastery.utils_subfigure import add_letter_and_title
from utils.postprocessing import annotate_results
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("dataset")
args = parser.parse_args()

dataset = args.dataset

extensions = ["cocoatree_SCA_none.csv",
              "cocoatree_MI_none.csv",
              "cocoatree_NMI_none.csv",
              "cocoatree_MI_APC.csv"]

start_i = 0
titles = ["SCA", "MI", "NMI", "MI+APC"]
letters = ["A.", "B.", "C.", "D."]
axes = [[0, 0], [0, 1], [1, 0], [1, 1]]

# Plot coevolution matrix
fig, axs = plt.subplots(nrows=2, ncols=2, layout='constrained')
for i, ext in enumerate(extensions):

    # Load and process coevolution matrix
    matrix = str(ext.split('.')[0] + '-distance.csv')
    coev_mat = pd.read_csv(f"results/cocoatree_gt/{dataset}/{matrix}")
    coev_mat = coev_mat.to_numpy()
    # Load and process results file
    results = pd.read_csv(f"results/cocoatree_gt/{dataset}/{ext}")
    results = annotate_results(results)
    results = results.loc[~results["pdb_pos"].isna()]
    results.filtered_msa_pos = results.filtered_msa_pos.astype('int')

    add_letter_and_title(axs[axes[i][0], axes[i][1]], letters[start_i],
                         title=titles[start_i])
    plot_coev_mat_sectors(fig, axs[axes[i][0], axes[i][1]], results, coev_mat,
                          title="")
    start_i += 1

os.makedirs(f"figures/{dataset}", exist_ok=True)
fig.savefig(f"figures/{dataset}/coev_matrix_sectors_all_metrics.png", dpi=300)
fig.savefig(f"figures/{dataset}/coev_matrix_sectors_all_metrics.pdf")

# Plot coevolution matrix without global mode (ngm)
start_i = 0
fig, axs = plt.subplots(nrows=2, ncols=2, layout='constrained')
for i, ext in enumerate(extensions):
    # Load and process coevolution matrix without global mode
    matrix_ngm = str(ext.split('.')[0] + '-distance_ngm.csv')
    coev_mat_ngm = pd.read_csv(f"results/cocoatree_gt/{dataset}/{matrix_ngm}")
    coev_mat_ngm = coev_mat_ngm.to_numpy()
    # Load and process results file
    results = pd.read_csv(f"results/cocoatree_gt/{dataset}/{ext}")
    results = annotate_results(results)
    results = results.loc[~results["pdb_pos"].isna()]
    results.filtered_msa_pos = results.filtered_msa_pos.astype('int')

    add_letter_and_title(axs[axes[i][0], axes[i][1]], letters[start_i],
                         title=titles[start_i])
    plot_coev_mat_sectors(fig, axs[axes[i][0], axes[i][1]], results,
                          coev_mat_ngm,
                          title="No global mode")
    start_i += 1

os.makedirs(f"figures/{dataset}", exist_ok=True)
fig.savefig(f"figures/{dataset}/coev_matrix_ngm_sectors_all_metrics.png",
            dpi=300)
fig.savefig(f"figures/{dataset}/coev_matrix_ngm_sectors_all_metrics.pdf")
