import matplotlib.pyplot as plt
from cocoatree.datasets import load_S1A_serine_proteases
from cocoatree.msa import filter_sequences
from cocoatree.statistics.pairwise import compute_sca_matrix
from cocoatree.deconvolution import extract_principal_components, \
    extract_independent_components, extract_xcors_from_ICs, \
    remove_global_correlations
from plotmastery.utils_subfigure import add_letter_and_title
import numpy as np

serine_dataset = load_S1A_serine_proteases(paper='halabi')
loaded_seqs = serine_dataset["alignment"]
loaded_seqs_id = serine_dataset["sequence_ids"]
n_loaded_pos, n_loaded_seqs = len(loaded_seqs[0]), len(loaded_seqs)

seq_kept, seq_id_kept, pos_kept = filter_sequences(loaded_seqs, loaded_seqs_id)

SCA_matrix = compute_sca_matrix(seq_kept)

#################
# Plot SCA matrix
fig, axes = plt.subplots(figsize=(7.5, 2.7), nrows=1, ncols=3, squeeze=False,
                         tight_layout=True)
ax = axes[0, 0]
im = ax.imshow(SCA_matrix, vmin=0, vmax=1.4, cmap='inferno')

ax.set_xlabel('Residues', fontweight="bold", fontsize="small", labelpad=2)
ax.set_ylabel('Residues', fontweight="bold", fontsize="small", labelpad=2)
# fig.colorbar(im, shrink=0.7)
add_letter_and_title(axes[0, 0], "A.", "SCA matrix")

n_components = 3
principal_components = extract_principal_components(SCA_matrix)
idpt_components = extract_independent_components(SCA_matrix,
                                                 n_components=n_components)
xcors = extract_xcors_from_ICs(idpt_components, SCA_matrix)


xcor_sizes = [len(x) for x in xcors]
cumul_sizes = sum(xcor_sizes)
sorted_pos = [p for xcor in xcors for p in xcor]

###################################################
# Plot reduced SCA matrix sorted according to XCoRs
ax = axes[0, 1]
im = ax.imshow(SCA_matrix[np.ix_(sorted_pos, sorted_pos)],
               vmin=0, vmax=2,
               interpolation='none', aspect='equal',
               extent=[0, cumul_sizes, cumul_sizes, 0],
               cmap='inferno')
# cb = fig.colorbar(im)
# cb.set_label("coevolution level")

line_index = 0
n_xcors = len(xcors)
for i in range(n_xcors):
    ax.plot([line_index + xcor_sizes[i], line_index + xcor_sizes[i]],
            [0, cumul_sizes], 'w', linewidth=2)
    ax.plot([0, cumul_sizes],
            [line_index + xcor_sizes[i], line_index + xcor_sizes[i]],
            'w', linewidth=2)
    line_index += xcor_sizes[i]

ticks = []
for ix in range(len(xcors)):
    shift = np.sum([len(xcors[j]) for j in range(ix)])
    ticks.append(shift+len(xcors[ix])/2)

ax.set_xticks(ticks)
ax.set_xticklabels(['XCoR_%d' % ix for ix in range(1, len(xcors)+1)],
                   fontweight="bold", fontsize="small")
ax.set_yticks(ticks)
ax.set_yticklabels(['XCoR_%d' % ix for ix in range(1, len(xcors)+1)],
                   rotation=90, va='center',
                   fontweight="bold", fontsize="small")
add_letter_and_title(axes[0, 1], "B.", "XCoR SCA matrix")

#####################
# Remove global mode
SCA_matrix_ngm = remove_global_correlations(SCA_matrix)
ax = axes[0, 2]
im = ax.imshow(SCA_matrix_ngm[np.ix_(sorted_pos, sorted_pos)],
               vmin=0, vmax=2,
               interpolation='none', aspect='equal',
               extent=[0, cumul_sizes, cumul_sizes, 0],
               cmap='inferno')
# cb = fig.colorbar(im)
# cb.set_label("coevolution level")

line_index = 0
n_xcors = len(xcors)
for i in range(n_xcors):
    ax.plot([line_index + xcor_sizes[i], line_index + xcor_sizes[i]],
            [0, cumul_sizes], 'w', linewidth=2)
    ax.plot([0, cumul_sizes],
            [line_index + xcor_sizes[i], line_index + xcor_sizes[i]],
            'w', linewidth=2)
    line_index += xcor_sizes[i]

ticks = []
for ix in range(len(xcors)):
    shift = np.sum([len(xcors[j]) for j in range(ix)])
    ticks.append(shift+len(xcors[ix])/2)

ax.set_xticks(ticks)
ax.set_xticklabels(['XCoR_%d' % ix for ix in range(1, len(xcors)+1)],
                   fontweight="bold", fontsize="small")
ax.set_yticks(ticks)
ax.set_yticklabels(['XCoR_%d' % ix for ix in range(1, len(xcors)+1)],
                   rotation=90, va='center',
                   fontweight="bold", fontsize="small")
add_letter_and_title(axes[0, 2], "C.", "without global mode")

fig.savefig("figures/figure_2.pdf")
fig.savefig("figures/figure_2.png", dpi=300)
