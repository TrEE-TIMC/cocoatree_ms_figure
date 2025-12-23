import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from plotmastery.utils_subfigure import add_letter_and_title

from cocoatree.datasets import load_S1A_serine_proteases
from cocoatree.io import load_MSA
from cocoatree.msa import filter_seq_id, compute_seq_identity

halabi = load_S1A_serine_proteases('halabi')
halabi_seq = halabi['alignment']
halabi_id = halabi['sequence_ids']

sector_file = '/home/jullimar/Documents/Postdoc_TIMC/2023-margaux-cocoatree/scripts/cocoatree_with_groundtruth/results/cocoatree_gt/halabi/cocoatree_xcor_1_SCA_none.fasta'
sector = load_MSA(sector_file, 'fasta')
sector_seq = sector['alignment']
sector_id = sector['sequence_ids']

metadata = '/home/jullimar/Documents/Postdoc_TIMC/2023-margaux-cocoatree/data/Trypsin/Halabi/halabi_metadata.csv'
df_annot = pd.read_csv(metadata)

trypsin = df_annot[df_annot['Protein_type'] == 'Trypsin'].Seq_ID
trypsin = list(map(str, trypsin))
chymotrypsin = df_annot[df_annot['Protein_type'] == 'Chymotrypsin'].Seq_ID
chymotrypsin = list(map(str, chymotrypsin))

_, trypsin_seq_id, trypsin_seq = filter_seq_id(halabi_seq, halabi_id, trypsin)
_, trypsin_xcor_id, trypsin_xcor_seq = filter_seq_id(sector_seq, sector_id,
                                                     trypsin)
_, chymotrypsin_seq_id, chymotrypsin_seq = filter_seq_id(halabi_seq, halabi_id,
                                                         chymotrypsin)
_, chymotrypsin_xcor_id, chymotrypsin_xcor_seq = filter_seq_id(sector_seq,
                                                               sector_id,
                                                               chymotrypsin)

idmat_trypsin_sequences = compute_seq_identity(trypsin_seq)
idmat_trypsin_xcor = compute_seq_identity(trypsin_xcor_seq)
idmat_chymotrypsin_sequences = compute_seq_identity(chymotrypsin_seq)
idmat_chymotrypsin_xcor = compute_seq_identity(chymotrypsin_xcor_seq)

# Scatter plot
fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(8, 4), tight_layout=True)
ax = axes[0]
ax.plot(idmat_trypsin_sequences, idmat_trypsin_xcor, marker='o', linewidth=0,
        markersize=4, zorder=50, color='#ff0000', alpha=0.2)
ax.plot([0, 1], [0, 1], transform=ax.transAxes, color='black', ls='dashed')
ax.tick_params(axis='both', which='both', labelsize='xx-small',
               bottom=True, top=False, labeltop=False, labelbottom=True,
               left=True, right=False, labelleft=True, labelright=False)
ax.spines[['right', 'top']].set_visible(False)
ax.set_xlim(0.5, 1)
ax.set_ylim(0, 1)
ax.set_xlabel('Full sequence', fontsize='x-small', fontweight='bold', labelpad=2)
ax.set_ylabel('XCoR sequence', fontsize='x-small', fontweight='bold', labelpad=2)
add_letter_and_title(ax, 'A.', 'Trypsin')

ax = axes[1]
ax.plot(idmat_chymotrypsin_sequences, idmat_chymotrypsin_xcor, marker='o',
        linewidth=0, markersize=4, zorder=50, color='#00a08a', alpha=0.2)
ax.plot([0, 1], [0, 1], transform=ax.transAxes, color='black', ls='dashed')
ax.tick_params(axis='both', which='both', labelsize='xx-small',
               bottom=True, top=False, labeltop=False, labelbottom=True,
               left=True, right=False, labelleft=True, labelright=False)
ax.spines[['right', 'top']].set_visible(False)
ax.set_xlim(0.5, 1)
ax.set_ylim(0, 1)
ax.set_xlabel('Full sequence', fontsize='x-small', fontweight='bold', labelpad=2)
ax.set_ylabel('XCoR sequence', fontsize='x-small', fontweight='bold', labelpad=2)
add_letter_and_title(ax, 'B.', 'Chymotrypsin')
fig.savefig('results/cocoatree_gt/halabi/trypsin_vs_chymotrypsin_pairwise_seqid_xcor_1_scatter.pdf')
fig.savefig('results/cocoatree_gt/halabi/trypsin_vs_chymotrypsin_pairwise_seqid_xcor_1_scatter.png')

# Keep only lower triangular identity matrix as array
tril_trypsin_xcor = idmat_trypsin_xcor[np.tril_indices(
    idmat_trypsin_xcor.shape[0])]
tril_trypsin_seq = idmat_trypsin_sequences[
    np.tril_indices(idmat_trypsin_sequences.shape[0])]
tril_chymotrypsin_xcor = idmat_chymotrypsin_xcor[np.tril_indices(
    idmat_chymotrypsin_xcor.shape[0])]
tril_chymotrypsin_seq = idmat_chymotrypsin_sequences[np.tril_indices(
    idmat_chymotrypsin_sequences.shape[0])]

# Distributions
fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(8, 4), tight_layout=True)
ax = axes[0]
bins = np.arange(0, 1, 0.05)
ax.hist(tril_trypsin_seq, bins=bins, alpha=0.7, label='Trypsins',
        color='#ff0000', density=True)
ax.hist(tril_chymotrypsin_seq, bins=bins, alpha=0.7, label='Chymotrypsins',
        color='#00a08a', density=True)
ax.legend(loc='upper left')
ax.tick_params(axis='both', which='both', labelsize='xx-small',
               bottom=True, top=False, labeltop=False, labelbottom=True,
               left=True, right=False, labelleft=True, labelright=False)
ax.set_xlabel('Pairwise sequence identity', fontsize='x-small',
              fontweight='bold', labelpad=2)
add_letter_and_title(ax, 'A.', 'Full sequences')

ax = axes[1]
ax.hist(tril_trypsin_xcor, bins=bins, alpha=0.7, label='Trypsins',
        color='#ff0000', density=True)
ax.hist(tril_chymotrypsin_xcor, bins=bins, alpha=0.7, label='Chymotrypsins',
        color='#00a08a', density=True)
ax.legend(loc='upper right')
ax.tick_params(axis='both', which='both', labelsize='xx-small',
               bottom=True, top=False, labeltop=False, labelbottom=True,
               left=True, right=False, labelleft=True, labelright=False)
ax.set_xlabel('Pairwise sequence identity', fontsize='x-small',
              fontweight='bold', labelpad=2)
add_letter_and_title(ax, 'B.', 'XCoR sequences')
fig.savefig('results/cocoatree_gt/halabi/trypsin_vs_chymotrypsin_pairwise_seqid_xcor_1_hist.pdf')
fig.savefig('results/cocoatree_gt/halabi/trypsin_vs_chymotrypsin_pairwise_seqid_xcor_1_hist.png')
fig.savefig('results/cocoatree_gt/halabi/trypsin_vs_chymotrypsin_pairwise_seqid_xcor_1_hist.svg')
