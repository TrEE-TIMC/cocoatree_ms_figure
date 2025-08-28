import pandas as pd
from cocoatree.io import load_MSA, load_tree_ete3
from cocoatree.visualization import update_tree_ete3_and_return_style
from utils.postprocessing import annotate_results
from utils.colors_and_labels import halabi_longer_cmap

dataset = 'halabi'
coevolution_metric = 'SCA'
correction = 'none'
tree = '/home/jullimar/Documents/Postdoc_TIMC/Trypsines/data/Halabi/5_ITOL/halabi_82_seqs_reordered_v2.txt'
annot_file = 'data/Trypsin/Halabi/halabi_metadata.csv'
df_annot = pd.read_csv(annot_file)
metadata_list = ['protein specificity', 'Subphylum', 'Class']

tree_ete3 = load_tree_ete3(tree)
# Load sector sequence as fasta file
sector_file = f"scripts/results/cocoatree_gt/{dataset}/cocoatree_sector_1_{coevolution_metric}_{correction}.fasta"
sector = load_MSA(sector_file, 'fasta')
sector_id = sector["sequence_ids"]
sector_seq = sector["alignment"]

tree_style = update_tree_ete3_and_return_style(
    tree_ete3,
    df_annot,
    sector_id,
    sector_seq,
    meta_data=metadata_list,
    show_leaf_name=False,
    fig_title="",
    linewidth=4,
    metadata_colors=halabi_longer_cmap,
    t_sector_seq=True,
    t_sector_heatmap=True,
    colormap='GnBl'
)

tree_ete3.render("", tree_style=tree_style)