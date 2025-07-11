import os

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

from plotmastery.utils_subfigure import add_letter_and_title

from utils.vis import plot_scatter_sectors
from utils.vis import plot_coev_mat_sectors
from utils.vis import plot_image
from utils.postprocessing import annotate_results
from utils.colors_and_labels import labels

import argparse

parser = argparse.ArgumentParser()
parser.add_argument("dataset")
parser.add_argument("--metric", "-m", default="SCA")
parser.add_argument("--correction", "-c", default="none")
args = parser.parse_args()


dataset = args.dataset
metric = args.metric
correction = args.correction
outname = f"figures/{dataset}/panel_{metric}_{correction}.pdf"


results_filename = \
    f"results/cocoatree_gt/{dataset}/cocoatree_{metric}_{correction}.csv"
results = pd.read_csv(results_filename)
results = annotate_results(results)

sca_sectors = [c for c in results.columns if c.startswith("sector")]
n_comp = len(sca_sectors)

fig = plt.figure(figsize=(8.3, 11.7))
gs = GridSpec(1106, 1000, figure=fig, top=0.95, left=0.1, right=0.99,
              bottom=0.05)


###############################################################################
# PC scatter plots

# First subpanel should go on half the page roughly, and then be splitted by
# n_comp

sep = 30
width = int(350 / (n_comp-1))
height = int(225 / (n_comp-1))

for i in range(n_comp):
    for j in range(n_comp):
        if j <= i:
            continue
        ax = fig.add_subplot(
            gs[i*(height+sep):(i+1)*height+i*sep,
               (j-1)*(width+sep):j*width+(j-1)*sep])
        plot_scatter_sectors(
           ax, results, f"PC{j+1}", f"PC{i+1}",
           annotate=False, add_labels=False)

        if i == j-1:
            ax.set_xlabel(f"PC{j+1}", fontweight="bold", fontsize="small",
                         labelpad=2)
            ax.set_ylabel(f"PC{i+1}", fontweight="bold", fontsize="small",
                          labelpad=2)
        else:
            ax.set_xticklabels([])
            ax.set_yticklabels([])

        if i == 0 and j == 1:
            add_letter_and_title(ax, "A.", "Principal components")

###############################################################################
# IC scatter plots

shift = 500

for i in range(n_comp):
    for j in range(n_comp):
        if j <= i:
            continue
        ax = fig.add_subplot(
            gs[i*(height+sep):(i+1)*height+i*sep,
               shift + (j-1)*(width+sep):shift + j*width+(j-1)*sep])
        plot_scatter_sectors(
           ax, results, f"IC{j+1}", f"IC{i+1}",
           annotate=False, add_labels=False)

        if i == j-1:
            ax.set_xlabel(f"{labels[j]}", fontweight="bold", fontsize="small",
                         labelpad=2)
            ax.set_ylabel(f"{labels[i]}", fontweight="bold", fontsize="small",
                          labelpad=2)
        else:
            ax.set_xticklabels([])
            ax.set_yticklabels([])

        if i == 0 and j == 1:
            add_letter_and_title(ax, "B.", "Independant components")

###############################################################################
# Sector coevolution matrix

col_start = (n_comp-1)*(height + sep) + 3*sep
height = 195

ax = fig.add_subplot(gs[col_start:col_start+height, :2*height])

coev_mat = pd.read_csv(results_filename.replace(".csv", "-distance.csv"))
coev_mat.columns = coev_mat.columns.astype(int)

results = pd.read_csv(results_filename)
results = annotate_results(results)

cb = plot_coev_mat_sectors(fig, ax, results, coev_mat)
add_letter_and_title(ax, "C.", "Sector coevolution")


ax = fig.add_subplot(gs[col_start:col_start+height, 500:500+2*height])

coev_mat = pd.read_csv(results_filename.replace(".csv", "-distance_ngm.csv"))
coev_mat.columns = coev_mat.columns.astype(int)

results = pd.read_csv(results_filename)
results = annotate_results(results)

cb = plot_coev_mat_sectors(fig, ax, results, coev_mat)
add_letter_and_title(ax, "D.", "Sector coevolution (without global mode)")

###############################################################################
# 3D structures
col_start = col_start + height + 3*sep

height = 200
width = int((1000-n_comp*sep)/n_comp)
trim = 5


for i, sector in enumerate(labels[:n_comp]):

    ax = fig.add_subplot(
        gs[col_start:col_start+height,
           i*(width+sep):(i+1)*width + i*sep])
    if i == 0:
        add_letter_and_title(ax, "E.", "3D structures")

    path = f"images/cocoatree_gt/{dataset}/cocoatree_{metric}_{correction}_{sector}_sector_1.png"
    plot_image(ax, path, trim=trim)


col_start = col_start + height
trim = 200

for i, sector in enumerate(labels[:n_comp]):

    ax = fig.add_subplot(
        gs[col_start:col_start+height,
           i*(width+sep):(i+1)*width + i*sep])

    path = f"images/cocoatree_gt/{dataset}/cocoatree_{metric}_{correction}_{sector}_sector_2.png"
    plot_image(ax, path, trim=trim)

os.makedirs(os.path.dirname(outname), exist_ok=True)
fig.savefig(outname, dpi=300)
