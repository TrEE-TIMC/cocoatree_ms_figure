import os
import pandas as pd
import argparse
from cocoatree.io import load_MSA, load_tree_ete3
from cocoatree.visualization import update_tree_ete3_and_return_style
from utils.postprocessing import annotate_results
from utils.colors_and_labels import halabi_cmap

parser = argparse.ArgumentParser()
parser.add_argument("dataset")
parser.add_argument("coevolution_metric")
parser.add_argument("correction")
parser.add_argument("tree")
args = parser.parse_args()

dataset = args.dataset
coevolution_metric = args.coevolution_metric
correction = args.correction
tree = args.tree

# Load metadata
if dataset == 'halabi':
    annot_file = '../data/Trypsin/Halabi/halabi_metadata.csv'
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
results = pd.read_csv(f"results/cocoatree_gt/{dataset}/cocoatree_{coevolution_metric}_{correction}.csv")
results = annotate_results(results)
results = results.loc[~results["pdb_pos"].isna()]
num_sectors = int((len(results.columns) - 11)/4)


for sect in range(1, num_sectors+1):
    # Load the tree each iteration to clean out tree_style
    tree_ete3 = load_tree_ete3(tree)
    # Load sector sequence as fasta file
    sector_file = f"results/cocoatree_gt/{dataset}/cocoatree_sector_{sect}_{coevolution_metric}_{correction}.fasta"
    sector = load_MSA(sector_file, 'fasta')
    sector_id = sector["sequence_ids"]
    sector_seq = sector["alignment"]

    tree_style = update_tree_ete3_and_return_style(
        tree_ete3,
        df_annot,
        sector_id,
        sector_seq,
        meta_data=['most_abundant_prot_specificity', 'vertebrate',
                   'most_abundant_class'],
        show_leaf_name=False,
        fig_title="",
        linewidth=4,
        metadata_colors=halabi_cmap,
        t_sector_seq=True,
        t_sector_heatmap=True,
        colormap='Blues'
    )

    os.makedirs(f"figures/{dataset}", exist_ok=True)
    tree_ete3.render(f"figures/{dataset}/tree_sector_{sect}_{coevolution_metric}_{correction}.png", tree_style=tree_style, dpi=300)
    tree_ete3.render(f"figures/{dataset}/tree_sector_{sect}_{coevolution_metric}_{correction}.pdf", tree_style=tree_style, dpi=300)
