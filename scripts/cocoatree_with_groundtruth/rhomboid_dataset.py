"""
Aligning the PDB with the alignment
"""


from Bio import Align


def align_pdb(data):
    aligner = Align.PairwiseAligner()

    pdb_seq = data["pdb_sequence"][0]
    ref_seq = data["alignment"][0]

    msa = aligner.align(pdb_seq, ref_seq)
    return msa[0]
