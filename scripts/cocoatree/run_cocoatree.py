import argparse
import os
import gzip
from cocoatree.io import load_MSA
from cocoatree import perform_sca


parser = argparse.ArgumentParser()
parser.add_argument("alignment_file")
parser.add_argument("--n-components", "-n", default=3)
parser.add_argument("--coevolution-metric", "-m", default="SCA")
parser.add_argument("--correction", "-c", default=None)
parser.add_argument("--outdir", "-o", default=None)
args = parser.parse_args()

filename = args.alignment_file
n_components = args.n_components
coevolution_metric = args.coevolution_metric
correction = args.correction
correction = correction if correction != "none" else None
outdir = args.outdir

# Load data
with gzip.open(filename, "rt") as f:
    data = load_MSA(f)

# Perform cocoatree analysis
coevolution_matrix, results = perform_sca(
    data["sequence_ids"], data["alignment"],
    n_components=n_components,
    coevolution_metric=coevolution_metric,
    correction=correction)


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
