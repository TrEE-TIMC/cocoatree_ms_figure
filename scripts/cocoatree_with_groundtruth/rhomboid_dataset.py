"""
Aligning the PDB with the alignment
"""


from cocoatree import datasets
from Bio import Align


def align_rhomboid_pdb():
    data = datasets.load_rhomboid_proteases()
    aligner = Align.PairwiseAligner()

    pdb_seq = data["pdb_sequence"][0]
    ref_seq = data["alignment"][0]

    msa = aligner.align(pdb_seq, ref_seq)
    return msa[0]
