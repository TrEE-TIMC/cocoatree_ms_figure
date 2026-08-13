# DCA Results

from cocoatree.datasets import load_S1A_serine_proteases
from cocoatree.io import load_MSA, load_tree_ete3
from cocoatree.visualization import update_tree_ete3_and_return_style
from utils.colors_and_labels import halabi_cmap

# Get original MSA
serine_dataset = load_S1A_serine_proteases(paper='halabi')
loaded_seqs = serine_dataset["alignment"]
loaded_seqs_id = serine_dataset["sequence_ids"]
n_loaded_pos, n_loaded_seqs = len(loaded_seqs[0]), len(loaded_seqs)
df_annot = serine_dataset["metadata"]

# Positions of DCA couples on original MSA
couples = [[507, 673], [248, 262], [218, 244], [216, 510], [261, 402],
           [281, 401], [660, 716], [250, 410], [279, 666], [263, 687]]

couple_seq = []
for seq in range(len(loaded_seqs)):
    init = ''
    for cpl in couples:
        if cpl != [263, 687]:
            init += loaded_seqs[seq][cpl[0]] + loaded_seqs[seq][cpl[1]] + '---'
        else:
            init += loaded_seqs[seq][cpl[0]] + loaded_seqs[seq][cpl[1]]
    couple_seq.append(init)

# Write a fasta file of the ten first DCA couples
with open('dca_couples.fasta', 'w') as outfasta:
    for seq in range(len(loaded_seqs)):
        outfasta.write('>' + str(loaded_seqs_id[seq]) + '\n')
        outfasta.write(couple_seq[seq] + '\n')

# Plot the DCA couples along the phylogenetic tree
fasta = load_MSA('dca_couples.fasta', 'fasta')
fasta_id = fasta['sequence_ids']
fasta_seq = fasta['alignment']

tree_file = '../../data/Trypsin/Halabi/halabi_cov_80_id_35_aln.fasta.treefile'
tree_ete3 = load_tree_ete3(tree_file)

tree_style, column_layout = update_tree_ete3_and_return_style(
    tree_ete3, df_annot,
    xcor_id=fasta_id,
    xcor_seq=fasta_seq,
    meta_data=('Protein_type', 'Subphylum', 'Class'),
    fig_title='Best 10 DCA couples',
    metadata_colors=halabi_cmap,
    linewidth=3,
    show_leaf_name=False,
    t_xcor_seq=True,
    t_xcor_heatmap=False,
    colormap='GnBu'
    )
tree_ete3.render('tree_DCA_couples.pdf', tree_style=tree_style)
