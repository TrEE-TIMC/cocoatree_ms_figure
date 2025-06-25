import argparse
import os

import pandas as pd
import numpy as np
from cocoatree import datasets
from cocoatree import perform_sca
from cocoatree.statistics import position

from utils.sectors import get_best_ordered_sectors


parser = argparse.ArgumentParser()
parser.add_argument("dataset")
parser.add_argument("--coevolution-metric", "-m", default="SCA")
parser.add_argument("--correction", "-c", default=None)
parser.add_argument("--outdir", "-o", default=None)
args = parser.parse_args()

dataset = args.dataset
coevolution_metric = args.coevolution_metric
correction = args.correction
correction = correction if correction != "none" else None
outdir = args.outdir

# Load data
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

# Compute conservation
conservation = position.compute_conservation(data["alignment"])

# Perform cocoatree analysis
coevolution_matrix, results = perform_sca(
    data["sequence_ids"], data["alignment"],
    n_components=n_components,
    coevolution_metric=coevolution_metric,
    correction=correction)

# Add original sectors to the results files
sectors = [
    [str(i) for i in data["sector_positions"][key]]
    for key in data["sector_positions"].keys()]

pdb_pos = data["pdb_positions"]
if dataset in ["rhomboid", "DHFR"]:
    from rhomboid_dataset import align_pdb
    msa = align_pdb(data)
    pdb_mapping = msa.indices[0][msa.indices[1] != -1]
    pdb_pos_mapping = np.array(pdb_pos)[pdb_mapping]
else:
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
    rename[f"sector_{j+1}"] = f"sector_{i+1}"
results.rename(rename, axis=1, inplace=True)

# Add statistics on the MSA such as conservation / number of gaps
results["msa_conservation"] = conservation
al = np.array([[pos for pos in seq] for seq in data["alignment"]])
perc_gap = (al == "-").sum(axis=0) / len(al) * 100
results["msa_perc_gap"] = perc_gap

# Write output
if outdir is not None:
    os.makedirs(outdir, exist_ok=True)
    outname = f"cocoatree_{coevolution_metric}"
    if correction is not None:
        outname = outname + f"_{correction}.csv"
    else:
        outname = outname + "_none.csv"
    outname = os.path.join(outdir, outname)
    results.to_csv(outname, index=False)

    outname = outname.replace(".csv", "-distance.csv")
    pd.DataFrame(data=coevolution_matrix).to_csv(outname, index=False)
