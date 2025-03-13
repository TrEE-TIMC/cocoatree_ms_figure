import os
import matplotlib.pyplot as plt
import pandas as pd
import argparse
from utils.vis import plot_scatter_sectors

parser = argparse.ArgumentParser()
parser.add_argument("filename")
parser.add_argument("--outname", "-o")
args = parser.parse_args()

filename = args.filename
outname = args.outname


results = pd.read_csv(filename)
results = results.loc[~results["filtered_msa_pos"].isna()]

sca_sectors = [c for c in results.columns if c.startswith("sector")]
orig_sectors = [c for c in results.columns if c.startswith("orig_")]

results["is_cocoatree_sector"] = results[sca_sectors].sum(axis=1).astype(bool)
results["is_orig_sector"] = results[orig_sectors].sum(axis=1).astype(bool)
results["is_both"] = results["is_cocoatree_sector"] & results["is_orig_sector"]
results["is_only_cocoatree"] = (
    results["is_cocoatree_sector"] & ~results["is_orig_sector"])
results["is_only_orig"] = (
    ~results["is_cocoatree_sector"] & results["is_orig_sector"])


fig, axes = plt.subplots(nrows=2, figsize=(6, 8))

ax = axes[0]

plot_scatter_sectors(ax, results, "IC1", "IC2")
ax.legend(frameon=False)

ax = axes[1]
plot_scatter_sectors(ax, results, "IC1", "IC3")

if outname is not None:
    os.makedirs(os.path.dirname(outname), exist_ok=True)
    fig.savefig(outname)
