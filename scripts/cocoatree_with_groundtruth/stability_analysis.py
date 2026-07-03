"""
Perform a stability analysis
"""

import itertools

import pandas as pd
import numpy as np
from cocoatree import datasets
from cocoatree import perform_sca
from cocoatree.statistics import position

import matplotlib.pyplot as plt

from utils.sectors import get_best_ordered_sectors
from utils.sectors import compute_all_vs_all
from utils.vis import sectors_cm
from plotmastery import utils_heatmap
from plotmastery.utils_subfigure import add_letter_and_title

from joblib import Memory

# Here, I'm going to cache some results
mem = Memory(".joblib")


dataset = "DHFR"
coevolution_metric = "SCA"
correction = None

sequence_thresholds = [0.1, 0.2, 0.3]
gap_thresholds = [0.3, 0.4, 0.5]

letters = ["A.", "B.", "C.", "D."]

def compute_sca(dataset, seq_thres=.2, gap_thres=.3):
    if dataset == "rivoire":
        data = datasets.load_S1A_serine_proteases(paper="rivoire")
    elif dataset == "halabi":
        data = datasets.load_S1A_serine_proteases(paper="halabi")
    elif dataset == "rhomboid":
        data = datasets.load_rhomboid_proteases()
    elif dataset == "DHFR":
        data = datasets.load_DHFR()

    n_components = len([k for k in data["sector_positions"].keys()])
    if dataset == "rhomboid":
        n_components = 3

    n_components = len([k for k in data["sector_positions"].keys()])
    # Compute conservation
    conservation = position.compute_conservation(data["alignment"])
    # Perform cocoatree analysis
    coevolution_matrix, coevolution_matrix_ngm, results = perform_sca(
        data["sequence_ids"], data["alignment"],
        n_components=n_components,
        gap_threshold=gap_thres,
        seq_threshold=seq_thres,
        coevolution_metric=coevolution_metric,
        correction=correction)

    # Add original sectors to the results files
    sectors = [
        [str(i) for i in data["sector_positions"][key]]
        for key in data["sector_positions"].keys()]

    pdb_pos = data["pdb_positions"]
    is_mapped = np.array([s != "-" for s in data["alignment"][0]])
    pdb_mapping = [int(val) if f else None
                   for f, val in zip(
                    is_mapped, (is_mapped.cumsum()-1))]
    pdb_pos_mapping = [
        pdb_pos[j]
        if i else None
        for i, j in zip(is_mapped, is_mapped.cumsum()-1)]

    results["pdb_pos"] = pdb_mapping
    results["pdb_named_pos"] = pdb_pos_mapping
    n_components = len([k for k in data["sector_positions"].keys()])

    mask = results["pdb_named_pos"].isna()
    for sector_id in range(n_components):
        results[f"orig_sector_{sector_id+1}"] = np.isin(
            results["pdb_named_pos"],
            sectors[sector_id])
        results.loc[mask, f"orig_sector_{sector_id+1}"] = False


    # Now reorder based on best match to original results

    order = get_best_ordered_sectors(results, dataset=dataset,
                                    type=coevolution_metric,
                                    correction=correction)
    rename = {}
    for i, j in enumerate(order):
        rename[f"IC{j+1}"] = f"IC{i+1}"
        rename[f"xcor_{j+1}"] = f"xcor_{i+1}"
    results.rename(rename, axis=1, inplace=True)

    # Add statistics on the MSA such as conservation / number of gaps
    results["msa_conservation"] = conservation
    al = np.array([[pos for pos in seq] for seq in data["alignment"]])
    perc_gap = (al == "-").sum(axis=0) / len(al) * 100
    results["msa_perc_gap"] = perc_gap
    return results

all_results = []

if dataset == "DHFR":
    sector_cols = ["xcor_1", "xcor_2", "xcor_3", "xcor_4"]
else:
    sector_cols = ["xcor_1", "xcor_2", "xcor_3"]


labels = []
for gap_thres, seq_thres in itertools.product(gap_thresholds, sequence_thresholds):
    res = mem.cache(compute_sca)(dataset, seq_thres, gap_thres)
    all_results.append([np.where(res[c])[0] for c in sector_cols if c in res.columns])
    labels.append(f"gap {gap_thres} - seq {seq_thres}")


fig, axes = plt.subplots(
    ncols=len(sector_cols),
    nrows=len(sector_cols),
    figsize=(7.5, 7.5),
    tight_layout=True)
cmaps = list(sectors_cm.values())

for i, xcor1 in enumerate(sector_cols): 
    for j, xcor2 in enumerate(sector_cols):
        overlap = compute_all_vs_all(all_results, comp1=i, comp2=j)
        ax = axes[i, j]
        if i == j:
            cmap = cmaps[i]
        else:
            cmap = "Greys"
        m = ax.matshow(overlap, vmin=0, cmap=cmap)
        utils_heatmap.annotate_heatmap(
                m, valfmt="{x:1.0f}",
                fontsize="x-small")
        if i == (len(sector_cols) - 1):
            ax.set_xticks(np.arange(len(labels)), labels, fontsize="x-small", rotation=90)
        else:
            ax.set_xticks([])
        if j == 0:
            ax.set_yticks(np.arange(len(labels)), labels, fontsize="x-small")
        else:
            ax.set_yticks([])

        ax.tick_params(
                axis='both', which='both', labelsize='x-small',
                bottom=True, top=False, labeltop=False, labelbottom=True,
                left=True, right=False, labelright=False, labelleft=True)

