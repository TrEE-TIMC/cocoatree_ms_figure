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

sector_file = 'results/cocoatree_gt/halabi/cocoatree_xcor_1_SCA_none.fasta'
sector = load_MSA(sector_file, 'fasta')
sector_seq = sector['alignment']
sector_id = sector['sequence_ids']

metadata = '../../data/Trypsin/Halabi/halabi_metadata.csv'
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
np.median(idmat_trypsin_xcor)
# np.float64(0.6086956521739131)
idmat_chymotrypsin_sequences = compute_seq_identity(chymotrypsin_seq)
idmat_chymotrypsin_xcor = compute_seq_identity(chymotrypsin_xcor_seq)
np.median(idmat_chymotrypsin_xcor)
# np.float64(0.30434782608695654)

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
fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(8, 3), tight_layout=True)
ax = axes[0]
bins = np.arange(0, 1, 0.05)
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
add_letter_and_title(ax, 'A.', 'XCoR sequences')

ax = axes[1]
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
add_letter_and_title(ax, 'B.', 'Full sequences')
fig.savefig('figures/halabi/trypsin_vs_chymotrypsin_pw_seqid_xcor_1_hist.pdf')
fig.savefig('figures/halabi/trypsin_vs_chymotrypsin_pw_seqid_xcor_1_hist.svg')

# Same thing but with vertebrate vs invertebrate in Green XCoR
sector_file = 'results/cocoatree_gt/halabi/cocoatree_xcor_2_SCA_none.fasta'
sector = load_MSA(sector_file, 'fasta')
sector_seq = sector['alignment']
sector_id = sector['sequence_ids']

vertebrate = df_annot[df_annot['Subphylum'] == 'vertebrate'].Seq_ID
vertebrate = list(map(str, vertebrate))
invertebrate = df_annot[df_annot['Subphylum'] == 'invertebrate'].Seq_ID
invertebrate = list(map(str, invertebrate))

_, vertebrate_seq_id, vertebrate_seq = filter_seq_id(halabi_seq, halabi_id,
                                                     vertebrate)
_, vertebrate_xcor_id, vertebrate_xcor_seq = filter_seq_id(sector_seq,
                                                           sector_id,
                                                           vertebrate)
_, invertebrate_seq_id, invertebrate_seq = filter_seq_id(halabi_seq, halabi_id,
                                                         invertebrate)
_, invertebrate_xcor_id, invertebrate_xcor_seq = filter_seq_id(sector_seq,
                                                               sector_id,
                                                               invertebrate)

idmat_vertebrate_sequences = compute_seq_identity(trypsin_seq)
idmat_vertebrate_xcor = compute_seq_identity(trypsin_xcor_seq)
idmat_invertebrate_sequences = compute_seq_identity(invertebrate_seq)
idmat_invertebrate_xcor = compute_seq_identity(invertebrate_xcor_seq)

# Keep only lower triangular identity matrix as array
tril_vertebrate_xcor = idmat_vertebrate_xcor[np.tril_indices(
    idmat_vertebrate_xcor.shape[0])]
tril_vertebrate_seq = idmat_vertebrate_sequences[
    np.tril_indices(idmat_vertebrate_sequences.shape[0])]
tril_invertebrate_xcor = idmat_invertebrate_xcor[np.tril_indices(
    idmat_invertebrate_xcor.shape[0])]
tril_invertebrate_seq = idmat_invertebrate_sequences[np.tril_indices(
    idmat_invertebrate_sequences.shape[0])]

fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(8, 3), tight_layout=True)
ax = axes[0]
bins = np.arange(0, 1, 0.05)
ax.hist(tril_vertebrate_xcor, bins=bins, alpha=0.7, label='Vertebrates',
        color='#798e87', density=True)
ax.hist(tril_invertebrate_xcor, bins=bins, alpha=0.7, label='Not vertebrates',
        color='#c27d38', density=True)
ax.legend(loc='upper right')
ax.tick_params(axis='both', which='both', labelsize='xx-small',
               bottom=True, top=False, labeltop=False, labelbottom=True,
               left=True, right=False, labelleft=True, labelright=False)
ax.set_xlabel('Pairwise sequence identity', fontsize='x-small',
              fontweight='bold', labelpad=2)
add_letter_and_title(ax, 'A.', 'XCoR sequences')

ax = axes[1]
ax.hist(tril_vertebrate_seq, bins=bins, alpha=0.7, label='Vertebrates',
        color='#798e87', density=True)
ax.hist(tril_invertebrate_seq, bins=bins, alpha=0.7, label='Not vertebrates',
        color='#c27d38', density=True)
ax.legend(loc='upper left')
ax.tick_params(axis='both', which='both', labelsize='xx-small',
               bottom=True, top=False, labeltop=False, labelbottom=True,
               left=True, right=False, labelleft=True, labelright=False)
ax.set_xlabel('Pairwise sequence identity', fontsize='x-small',
              fontweight='bold', labelpad=2)
add_letter_and_title(ax, 'B.', 'Full-length sequences')
fig.savefig('figures/halabi/vertebrate_vs_invert_pw_seqid_xcor_2_hist.pdf')
fig.savefig('figures/halabi/vertebrate_vs_invert_pw_seqid_xcor_2_hist.svg')

# Same thing but with vertebrate vs invertebrate in Blue XCoR
sector_file = 'results/cocoatree_gt/halabi/cocoatree_xcor_3_SCA_none.fasta'
sector = load_MSA(sector_file, 'fasta')
sector_seq = sector['alignment']
sector_id = sector['sequence_ids']

vertebrate = df_annot[df_annot['Subphylum'] == 'vertebrate'].Seq_ID
vertebrate = list(map(str, vertebrate))
invertebrate = df_annot[df_annot['Subphylum'] == 'invertebrate'].Seq_ID
invertebrate = list(map(str, invertebrate))

_, vertebrate_seq_id, vertebrate_seq = filter_seq_id(halabi_seq, halabi_id,
                                                     vertebrate)
_, vertebrate_xcor_id, vertebrate_xcor_seq = filter_seq_id(sector_seq,
                                                           sector_id,
                                                           vertebrate)
_, invertebrate_seq_id, invertebrate_seq = filter_seq_id(halabi_seq, halabi_id,
                                                         invertebrate)
_, invertebrate_xcor_id, invertebrate_xcor_seq = filter_seq_id(sector_seq,
                                                               sector_id,
                                                               invertebrate)

idmat_vertebrate_sequences = compute_seq_identity(trypsin_seq)
idmat_vertebrate_xcor = compute_seq_identity(trypsin_xcor_seq)
idmat_invertebrate_sequences = compute_seq_identity(invertebrate_seq)
idmat_invertebrate_xcor = compute_seq_identity(invertebrate_xcor_seq)

# Keep only lower triangular identity matrix as array
tril_vertebrate_xcor = idmat_vertebrate_xcor[np.tril_indices(
    idmat_vertebrate_xcor.shape[0])]
tril_vertebrate_seq = idmat_vertebrate_sequences[
    np.tril_indices(idmat_vertebrate_sequences.shape[0])]
tril_invertebrate_xcor = idmat_invertebrate_xcor[np.tril_indices(
    idmat_invertebrate_xcor.shape[0])]
tril_invertebrate_seq = idmat_invertebrate_sequences[np.tril_indices(
    idmat_invertebrate_sequences.shape[0])]

fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(8, 4), tight_layout=True)
ax = axes[0]
bins = np.arange(0, 1, 0.05)
ax.hist(tril_vertebrate_xcor, bins=bins, alpha=0.7, label='Vertebrates',
        color='#798e87', density=True)
ax.hist(tril_invertebrate_xcor, bins=bins, alpha=0.7, label='Not vertebrates',
        color='#c27d38', density=True)
ax.legend(loc='upper right')
ax.tick_params(axis='both', which='both', labelsize='xx-small',
               bottom=True, top=False, labeltop=False, labelbottom=True,
               left=True, right=False, labelleft=True, labelright=False)
ax.set_xlabel('Pairwise sequence identity', fontsize='x-small',
              fontweight='bold', labelpad=2)
add_letter_and_title(ax, 'A.', 'XCoR sequences')

ax = axes[1]
ax.hist(tril_vertebrate_seq, bins=bins, alpha=0.7, label='Vertebrates',
        color='#798e87', density=True)
ax.hist(tril_invertebrate_seq, bins=bins, alpha=0.7, label='Not vertebrates',
        color='#c27d38', density=True)
ax.legend(loc='upper left')
ax.tick_params(axis='both', which='both', labelsize='xx-small',
               bottom=True, top=False, labeltop=False, labelbottom=True,
               left=True, right=False, labelleft=True, labelright=False)
ax.set_xlabel('Pairwise sequence identity', fontsize='x-small',
              fontweight='bold', labelpad=2)
add_letter_and_title(ax, 'B.', 'Full-length sequences')
fig.savefig('figures/halabi/vertebrate_vs_invert_pw_seqid_xcor_3_hist.pdf')
fig.savefig('figures/halabi/vertebrate_vs_invert_pw_seqid_xcor_3_hist.svg')
