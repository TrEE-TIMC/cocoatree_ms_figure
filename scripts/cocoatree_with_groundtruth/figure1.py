import matplotlib.pyplot as plt
import pandas as pd
from utils.vis import plot_scatter_sectors



def annotate_results(results):
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
    return results



fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(6, 6), tight_layout=True)

###############################################################################
# Start with halabi results
results = pd.read_csv("results/cocoatree_gt/halabi/cocoatree_SCA_none.csv")
results = annotate_results(results)

ax = axes[0, 0]
plot_scatter_sectors(ax, results, "IC1", "IC2", annotate=False)
ax.legend(frameon=False)

ax = axes[1, 0]
plot_scatter_sectors(ax, results, "IC1", "IC3", annotate=False)

###############################################################################
# rhomboid results
results = pd.read_csv("results/cocoatree_gt/rhomboid/cocoatree_SCA_none.csv")
results = annotate_results(results)

ax = axes[0, 1]
plot_scatter_sectors(ax, results, "IC1", "IC2", annotate=False)
ax.legend(frameon=False)

ax = axes[1, 1]
plot_scatter_sectors(ax, results, "IC1", "IC3", annotate=False)
