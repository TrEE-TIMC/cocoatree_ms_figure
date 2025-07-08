import pandas as pd

import argparse
import numpy as np

import matplotlib.pyplot as plt

from plotmastery.utils_subfigure import add_letter_and_title

from utils.postprocessing import annotate_results
from utils.vis import plot_coev_mat_sectors


parser = argparse.ArgumentParser()
parser.add_argument("coev_filename")
parser.add_argument("results_filename")
parser.add_argument("--letter", "-l", default="C.")
parser.add_argument("--title", "-t", default="Sector coevolution matrix")
parser.add_argument("--outname", "-o", default=None)
args = parser.parse_args()

coev_filename = args.coev_filename
results_filename = args.results_filename
letter = args.letter
title = args.title
outname = args.outname

coev_mat = pd.read_csv(coev_filename)
coev_mat.columns = coev_mat.columns.astype(int)

results = pd.read_csv(results_filename)
results = annotate_results(results)

fig, ax = plt.subplots(figsize=(4, 3), tight_layout=True)
cb = plot_coev_mat_sectors(fig, ax, results, coev_mat)
add_letter_and_title(ax, letter, title=title)


if outname is not None:
    fig.savefig(outname, dpi=300)
    fig.savefig(outname.replace(".png", ".pdf"))
