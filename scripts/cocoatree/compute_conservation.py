import argparse
import pandas as pd
import os
import gzip
from cocoatree.io import load_MSA
from cocoatree.statistics.position import compute_conservation
from cocoatree.msa import filter_sequences


parser = argparse.ArgumentParser()
parser.add_argument("alignment_file")
parser.add_argument("--outdir", "-o", default=None)
args = parser.parse_args()

filename = args.alignment_file
outdir = args.outdir


# Load data
with gzip.open(filename, "rt") as f:
    data = load_MSA(f)

seq_kept, seq_id_kept, pos_kept = filter_sequences(
    data["alignment"],
    data["sequence_ids"])

conservation = compute_conservation(seq_kept)

# Write output
if outdir is not None:
    os.makedirs(outdir, exist_ok=True)
    outname = "conservation.csv"
    outname = os.path.join(outdir, outname)
    conservation = pd.DataFrame(data=conservation, columns=("Conservation", ))
    conservation.to_csv(outname, index=False)
