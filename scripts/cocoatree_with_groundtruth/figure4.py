import pandas as pd
from cocoatree.io import load_MSA, load_tree_ete3
from cocoatree.visualization import update_tree_ete3_and_return_style
from utils.colors_and_labels import halabi_longer_cmap
import argparse

parser = argparse.ArgumentParser(prog='plot_figure_tree',
                                 description='This script allows to select different datasets, metrics, and corrections' \
                                 'to plot Cocoatree\'s figure of a phylogenetic tree along with sequence information and' \
                                 'a heatmap of sequence identity')
parser.add_argument("tree")
parser.add_argument("metadata")
parser.add_argument("sector_fasta")
parser.add_argument('metadata_list', metavar='N', type=str, nargs='+',
                    help='a list of metadata to display')
parser.add_argument('title')
parser.add_argument('cmap')
parser.add_argument("output")
args = parser.parse_args()

tree = args.tree
# tree = '/home/jullimar/Documents/Postdoc_TIMC/Trypsines/data/Halabi/4_IQTREE/halabi_subset_aln_kpsg.fasta.treefile'
metadata = args.metadata
# metadata = '/home/jullimar/Documents/Postdoc_TIMC/2023-margaux-cocoatree/data/Trypsin/Halabi/halabi_metadata.csv'
sector_file = args.sector_fasta
# sector_file = '/home/jullimar/Documents/Postdoc_TIMC/2023-margaux-cocoatree/scripts/results/cocoatree_gt/halabi/cocoatree_xcor_1_SCA_none.fasta'
# sector_file = '/home/jullimar/Documents/Postdoc_TIMC/halabi_xcor_1_SCA_none.fasta'
outname = args.output
# outname = '/home/jullimar/Documents/Postdoc_TIMC/test_plot_figure3C.png'
metadata_list = args.metadata_list
# metadata_list = ['Protein_type', 'Subphylum', 'Class']
title = args.title
cmap = args.cmap


#print(metadata_list)
#print(type(metadata_list))

# Load metadata file
df_annot = pd.read_csv(metadata)
# Load tree file
tree_ete3 = load_tree_ete3(tree)
# Load sector sequence as fasta file
sector = load_MSA(sector_file, 'fasta')
sector_id = sector["sequence_ids"]
sector_seq = sector["alignment"]
#subsector_seq = [res[:5] for res in sector_seq]

#subsector_seq = []
#for seq in sector_seq:
#    subseq = str(seq[1]+seq[2]+seq[4])
#    subsector_seq.append(subseq)


tree_style, _ = update_tree_ete3_and_return_style(
    tree_ete3,
    df_annot,
    sector_id,
    sector_seq,
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
