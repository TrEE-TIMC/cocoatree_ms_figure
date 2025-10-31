import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import argparse
from cocoatree.io import load_MSA
import logomaker


parser = argparse.ArgumentParser()
parser.add_argument("dataset")
parser.add_argument("coevolution_metric")
parser.add_argument("correction")
args = parser.parse_args()

dataset = args.dataset
coevolution_metric = args.coevolution_metric
correction = args.correction

# Load metadata
if dataset == 'halabi':
    annot_file = '../../data/Trypsin/Halabi/halabi_metadata.csv'
    metadata_list = ['']
elif dataset == 'rivoire':
    annot_file = '../data/Trypsin/Rivoire/rivoire_metadata.csv'
elif dataset == 'DHFR':
    annot_file = '../data/DHFR/DHFR_metadata.csv'
elif dataset == 'rhomboid':
    annot_file = '../data/Rhomboid/rhomboid_metadata_clean.csv'
df_annot = pd.read_csv(annot_file)


# Load results file in order to detect the number of sectors
# TODO: Try to find another way that doesn't require loading a file
results = pd.read_csv(
    f"results/cocoatree_gt/{dataset}/cocoatree_{coevolution_metric}_{correction}.csv")
num_sectors = len([col for col in results.columns if col.startswith("xcor")])
results = results.loc[~results["pdb_pos"].isna()]

cat_name = "prot_specificity"
categories = {

    "prot_specificity": ["chymotrypsin", "kallikrein", "trypsin", "tryptase",
                         "elastase"],

     "Subphylum":  ["vertebrate", "invertebrate", "bacteria", "fungi",
                     "virus"],
     "class":  ["Mammalia", "Actinopterygii", "Amphibia", "Malacostraca",
                "Insecta"]
    }


def get_comp_from_seq(sequences, characters="ACDEFGHIKLMNPQRSTVWY"):
    n_aa = len(sequences[0])
    compDict = {char: [0]*n_aa for char in characters}
    for seq in sequences:
        for aaPos in range(len(seq)):
            aa = seq[aaPos]
            if aa in characters:
                compDict[aa][aaPos] += 1
    return pd.DataFrame.from_dict(compDict)


colors = {
    "A": "#c8c8c8",
    "C": "#e5e501",
    "D": "#e60a0a",
    "E": "#e60a0a",
    "F": "#3232aa",
    "G": "#ebebeb",
    "H": "#8282d2",
    "I": "#0f820f",
    "K": "#145aff",
    "L": "#0f820f",
    "M": "#e5e501",
    "N": "#00dcdc",
    "P": "#dc9682",
    "Q": "#00dcdc",
    "R": "#145aff",
    "S": "#f79400",
    "T": "#f79400",
    "V": "#0f820f",
    "W": "#b45ab4",
    "Y": "#3232aa",
}


for sect in range(1, num_sectors+1):
    # Load sector sequence as fasta file
    sector_file = f"results/cocoatree_gt/{dataset}/cocoatree_xcor_{sect}_{coevolution_metric}_{correction}.fasta"
    sector = load_MSA(sector_file, 'fasta')
    sector_id = np.array(sector["sequence_ids"])
    sector_seq = np.array(sector["alignment"])
    fig, axes = plt.subplots(figsize=(8, 5), nrows=5, ncols=3,
                             tight_layout=True)
    for i, cat_name in enumerate(categories.keys()):
        for j, cat in enumerate(categories[cat_name]):
            ax = axes[j, i]

            cat_sectors = sector_seq[(df_annot[cat_name] == cat)]
            ppm = get_comp_from_seq(cat_sectors)

            crp_logo = logomaker.Logo(ppm, ax=ax, color_scheme=colors)
            crp_logo.style_spines(visible=False)
            crp_logo.style_spines(spines=['left', 'bottom'], visible=True)
            crp_logo.style_xticks(rotation=90, fmt='%d', anchor=0)

            crp_logo.ax.set_title(cat.capitalize() + f" ({len(cat_sectors)})",
                                  fontweight="bold", fontsize="x-small")

            residue_num = results.loc[
                results[f"xcor_{sect}"]].sort_values(
                    f"IC{sect}",
                    ascending=False)["pdb_named_pos"].values
            crp_logo.ax.set_xticklabels(
                residue_num[:len(cat_sectors[0])],
                fontsize="xx-small")
            crp_logo.ax.set_yticks([])
    fig.savefig(f"figures/{dataset}/logo_{sect}.pdf")
