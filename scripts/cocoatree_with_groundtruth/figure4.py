import pandas as pd
from cocoatree.io import load_MSA, load_tree_ete3
from cocoatree.visualization import update_tree_ete3_and_return_style
from utils.colors_and_labels import halabi_longer_cmap
import argparse

parser = argparse.ArgumentParser(prog='plot_figure_tree',
                                 description='This script allows to select \
                                     different datasets, metrics, and \
                                     corrections to plot Cocoatree\'s \
                                     figure of a phylogenetic tree \
                                     along with sequence information and \
                                     a heatmap of sequence identity')
parser.add_argument("tree",
                    help="path to the newick file")
parser.add_argument("metadata",
                    help="path to the csv of metadata associated to the tree")
parser.add_argument("xcor_fasta",
                    help="path to the fasta of the XCoR")
parser.add_argument('metadata_list', metavar='N', type=str, nargs='+',
                    help='a list of metadata to display')
parser.add_argument('title',
                    help="figure title")
parser.add_argument('cmap',
                    help="matplotlib colormap for the heatmap")
parser.add_argument("output",
                    help="path to the output")
args = parser.parse_args()

tree = args.tree
metadata = args.metadata
xcor_file = args.xcor_fasta
outname = args.output
metadata_list = args.metadata_list
title = args.title
cmap = args.cmap

# Load metadata file
df_annot = pd.read_csv(metadata)
# Load tree file
tree_ete3 = load_tree_ete3(tree)
# Load sector sequence as fasta file
xcor = load_MSA(xcor_file, 'fasta')
xcor_id = xcor["sequence_ids"]
xcor_seq = xcor["alignment"]


tree_style, _ = update_tree_ete3_and_return_style(
    tree_ete3,
    df_annot,
    xcor_id,
    xcor_seq,
    meta_data=metadata_list,
    show_leaf_name=False,
    fig_title=title,
    linewidth=4,
    metadata_colors=halabi_longer_cmap,
    t_xcor_seq=True,
    t_xcor_heatmap=True,
    colormap=cmap
)

tree_ete3.render(outname, tree_style=tree_style)
