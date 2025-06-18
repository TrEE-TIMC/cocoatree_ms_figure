import matplotlib.pyplot as plt
import pandas as pd
from joblib import Memory
from cocoatree import datasets
from cocoatree.msa import filter_sequences
from cocoatree.statistics.position import compute_conservation
from cocoatree.statistics import pairwise

mem = Memory(".joblib")

dataset = "halabi"
serine_dataset = datasets.load_S1A_serine_proteases(paper=dataset)
seq_id = serine_dataset["sequence_ids"]
sequences = serine_dataset["alignment"]
n_pos, n_seq = len(sequences[0]), len(sequences)

seq_kept, seq_id_kept, pos_kept = filter_sequences(sequences, seq_id)

# Load precomputed results
res = pd.read_csv(f"results/cocoatree_gt/{dataset}/cocoatree_SCA_none.csv")
res = res.loc[~res["filtered_msa_pos"].isna()] 

sector_cols = [col for col in res.columns if col.startswith("sector")]
is_in_sectors = res[sector_cols].sum(axis=1) > 0

###############################################################################
# Compute matrices and conservation

conservation = mem.cache(compute_conservation)(seq_kept)
sca_matrix = mem.cache(pairwise.compute_sca_matrix)(seq_kept)
mi_matrix = mem.cache(pairwise.compute_mutual_information_matrix)(
    seq_kept,
    normalize=False)
nmi_matrix = mem.cache(pairwise.compute_mutual_information_matrix)(seq_kept)
_, mi_apc_matrix = mem.cache(pairwise.compute_apc)(mi_matrix)

###############################################################################
# plot results
fig, axes = plt.subplots(figsize=(8, 2), nrows=1, ncols=4, squeeze=False,
tight_layout=True)

ax = axes[0, 0]

def plot_conservation_vs_matrix(ax, matrix, colors):
    ax.scatter(conservation, matrix.sum(axis=0),
               c=colors,
               marker=".", linewidth=0)
    ax.spines["right"].set_linewidth(0)
    ax.spines["top"].set_linewidth(0)
    ax.set_xlabel("Conservation", fontweight="bold", fontsize="small")


plot_conservation_vs_matrix(axes[0, 0], sca_matrix, sca_matrix.sum(axis=1))
plot_conservation_vs_matrix(axes[0, 1], mi_matrix, sca_matrix.sum(axis=1))
plot_conservation_vs_matrix(axes[0, 2], nmi_matrix, sca_matrix.sum(axis=1))
plot_conservation_vs_matrix(axes[0, 3], mi_apc_matrix, sca_matrix.sum(axis=1))
