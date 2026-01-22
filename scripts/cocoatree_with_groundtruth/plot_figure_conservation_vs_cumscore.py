import matplotlib.pyplot as plt
import pandas as pd
from joblib import Memory
from cocoatree import datasets
from cocoatree.msa import filter_sequences
from cocoatree.statistics.position import compute_conservation
from cocoatree.statistics import pairwise
from plotmastery.utils_subfigure import add_letter_and_title
from utils.colors_and_labels import colors
from utils.postprocessing import annotate_results
import argparse

mem = Memory(".joblib")

parser = argparse.ArgumentParser()
parser.add_argument("dataset")
args = parser.parse_args()

dataset = args.dataset

if dataset in ["rivoire", "halabi"]:
    data = datasets.load_S1A_serine_proteases(paper=dataset)
elif dataset == "rhomboid":
    data = datasets.load_rhomboid_proteases()
else:
    data = datasets.load_DHFR()

seq_id = data["sequence_ids"]
sequences = data["alignment"]
n_pos, n_seq = len(sequences[0]), len(sequences)

seq_kept, seq_id_kept, pos_kept = filter_sequences(sequences, seq_id)

# Load precomputed results
res = pd.read_csv(f"results/cocoatree_gt/{dataset}/cocoatree_SCA_none.csv")
res = res.loc[~res["filtered_msa_pos"].isna()]
res = annotate_results(res)

sector_cols = [col for col in res.columns if col.startswith("xcor")]
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
fig, axes = plt.subplots(figsize=(8.3, 2.2), nrows=1, ncols=4, squeeze=False,
                         tight_layout=True)


def plot_conservation_vs_matrix(ax, matrix, results=None):
    sector_columns = [c for c in results.columns if c.startswith("xcor")]
    sector_columns.sort()
    mask = ~results["is_only_cocoatree"]
    ax.scatter(conservation[mask], matrix.sum(axis=0)[mask],
               c=colors["default"], marker=".", linewidth=0)
    for c in sector_columns:
        sec_id = c.split("_")[-1]
        mask = res[c]
        c = colors[f"xcor_{sec_id}"] 
        ax.scatter(conservation[mask], matrix.sum(axis=0)[mask],
                   c=c,
                   marker=".", linewidth=0)
    ax.spines["right"].set_linewidth(0)
    ax.spines["top"].set_linewidth(0)
    ax.set_xlabel("Conservation", fontweight="bold", fontsize="small",
                  labelpad=2)
    ax.tick_params(
        axis='both', which='both', labelsize='x-small',
        bottom=True, top=False, labeltop=False, labelbottom=True,
        left=True, right=False, labelleft=True, labelright=False)


plot_conservation_vs_matrix(axes[0, 0], sca_matrix, results=res)
axes[0, 0].set_ylabel("Cum score", fontsize="small", fontweight="bold",
                      labelpad=2)
add_letter_and_title(axes[0, 0], "A.", "SCA")


res = pd.read_csv(f"results/cocoatree_gt/{dataset}/cocoatree_MI_none.csv")
res = res.loc[~res["filtered_msa_pos"].isna()]
res = annotate_results(res)
plot_conservation_vs_matrix(axes[0, 1], mi_matrix, results=res)
add_letter_and_title(axes[0, 1], "B.", "MI")

res = pd.read_csv(f"results/cocoatree_gt/{dataset}/cocoatree_NMI_none.csv")
res = res.loc[~res["filtered_msa_pos"].isna()]
res = annotate_results(res)
plot_conservation_vs_matrix(axes[0, 2], nmi_matrix, results=res)
add_letter_and_title(axes[0, 2], "C.", "NMI")

res = pd.read_csv(f"results/cocoatree_gt/{dataset}/cocoatree_MI_APC.csv")
res = res.loc[~res["filtered_msa_pos"].isna()]
res = annotate_results(res)
plot_conservation_vs_matrix(axes[0, 3], mi_apc_matrix, results=res)
add_letter_and_title(axes[0, 3], "D.", "MI+APC")
fig.savefig(f"figures/{dataset}/conservation_vs_cumscore.pdf")
