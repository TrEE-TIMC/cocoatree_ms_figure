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


dataset = "halabi"
metric = "SCA"
correction = "none"


results_filename = \
    f"results/cocoatree_gt/{dataset}/cocoatree_{metric}_{correction}.csv"
results = pd.read_csv(results_filename)
results = annotate_results(results)

sca_sectors = [c for c in results.columns if c.startswith("sector")]
n_comp = len(sca_sectors)


fig = plt.figure(figsize=(8.3, 11.7))
gs = GridSpec(1106, 1000, figure=fig, top=0.95, left=0.1, right=0.95)


###############################################################################
# PC scatter plots

# First subpanel should go on half the page roughly, and then be splitted by
# n_comp

sep = 30
width = int(350 / (n_comp-1))
height = int(225 / (n_comp-1))
ax = fig.add_subplot(gs[:height, :width])
plot_scatter_sectors(
    ax, results, "PC2", "PC1", annotate=False, add_labels=False)
ax.set_xlabel("PC2", fontweight="bold", fontsize="small",
              labelpad=2)
ax.set_ylabel("PC1", fontweight="bold", fontsize="small",
              labelpad=2)
add_letter_and_title(ax, "A.", "Principal components")

ax = fig.add_subplot(gs[:height, width+sep:2*width+sep])
plot_scatter_sectors(
    ax, results, "PC3", "PC1", annotate=False, add_labels=False)
ax.set_xticklabels([])
ax.set_yticklabels([])


ax = fig.add_subplot(gs[height+sep:2*height+sep, width+sep:2*width+sep])
plot_scatter_sectors(
    ax, results, "PC3", "PC2", annotate=False, add_labels=False)
ax.set_xlabel("PC3", fontweight="bold", fontsize="small",
              labelpad=2)
ax.set_ylabel("PC2", fontweight="bold", fontsize="small",
              labelpad=2)

###############################################################################
# IC scatter plots

ax = fig.add_subplot(gs[:height, 500+sep:500+sep+width])
plot_scatter_sectors(
    ax, results, "IC2", "IC1", annotate=False, add_labels=False)
ax.set_xlabel(labels[1], fontweight="bold", fontsize="small",
              labelpad=2)
ax.set_ylabel(labels[0], fontweight="bold", fontsize="small",
              labelpad=2)
add_letter_and_title(ax, "B.", "Independant components")

ax = fig.add_subplot(gs[:height, 500+sep*2+width:500+2*(sep+width)])
plot_scatter_sectors(
    ax, results, "IC3", "IC1", annotate=False, add_labels=False)
ax.set_xticklabels([])
ax.set_yticklabels([])


ax = fig.add_subplot(gs[height+sep:2*height+sep, 500+sep*2+width:500+2*(sep+width)])
plot_scatter_sectors(
    ax, results, "IC3", "IC2", annotate=False, add_labels=False)
ax.set_xlabel(labels[2], fontweight="bold", fontsize="small",
              labelpad=2)
ax.set_ylabel(labels[1], fontweight="bold", fontsize="small",
              labelpad=2)


###############################################################################
# Sector coevolution matrix

col_start = 2*height+4*sep

ax = fig.add_subplot(gs[col_start:col_start+2*height, :2*width+sep])

coev_mat = pd.read_csv(results_filename.replace(".csv", "-distance.csv"))
coev_mat.columns = coev_mat.columns.astype(int)

results = pd.read_csv(results_filename)
results = annotate_results(results)

cb = plot_coev_mat_sectors(fig, ax, results, coev_mat)
add_letter_and_title(ax, "C.", "Sector coevolution")


ax = fig.add_subplot(gs[col_start:col_start+2*height, 500:500+2*(width+sep)])

coev_mat = pd.read_csv(results_filename.replace(".csv", "-distance_ngm.csv"))
coev_mat.columns = coev_mat.columns.astype(int)

results = pd.read_csv(results_filename)
results = annotate_results(results)

cb = plot_coev_mat_sectors(fig, ax, results, coev_mat)
add_letter_and_title(ax, "D.", "Sector coevolution (without global mode)")

###############################################################################
# 3D structures
col_start = col_start + 2*height + 3*sep

width = int((1000-3*sep)/3)
trim = 250

ax = fig.add_subplot(gs[col_start:col_start+2*height, :width])
add_letter_and_title(ax, "E.", "3D structures")

path = f"images/cocoatree_gt/{dataset}/cocoatree_{metric}_{correction}_Red_sector_1.png"
plot_image(ax, path, trim=trim)

ax = fig.add_subplot(gs[col_start:col_start+2*height, width+sep:2*width+sep])
path = f"images/cocoatree_gt/{dataset}/cocoatree_{metric}_{correction}_Green_sector_1.png"
plot_image(ax, path, trim=trim)

ax = fig.add_subplot(
    gs[col_start:col_start+2*height,
       2*(width+sep):3*width+2*sep])
path = f"images/cocoatree_gt/{dataset}/cocoatree_{metric}_{correction}_Blue_sector_1.png"
plot_image(ax, path, trim=trim)


col_start = col_start + 2*height

width = int((1000-3*sep)/3)
trim = 200
ax = fig.add_subplot(gs[col_start:col_start+2*height, :width])
path = f"images/cocoatree_gt/{dataset}/cocoatree_{metric}_{correction}_Red_sector_2.png"
plot_image(ax, path, trim=trim)

ax = fig.add_subplot(gs[col_start:col_start+2*height, width+sep:2*width+sep])
path = f"images/cocoatree_gt/{dataset}/cocoatree_{metric}_{correction}_Green_sector_2.png"
plot_image(ax, path, trim=trim)

ax = fig.add_subplot(
    gs[col_start:col_start+2*height,
       2*(width+sep):3*width+2*sep])
path = f"images/cocoatree_gt/{dataset}/cocoatree_{metric}_{correction}_Blue_sector_2.png"
plot_image(ax, path, trim=trim)
