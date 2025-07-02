python
import numpy as np
import pandas as pd
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt

# #154c3b
colors = {
    "cocoatree": "#AB0000",
    "other": "#006c4b",
    "both": "#006ab2",
    "default": "#9cadbc",
    "sector_1": "#AB0000", # Red
    "sector_3": "#0023B2", #"#006ab2", # Green
    "sector_2": "#078400", #"#006c4b", # Blue
    "sector_4": "#845EC2",
    "sector_5": "crimson",
    "sector_6": "orange"}




residue_swap = pd.read_csv(
"/home/nelle/Projects/research/2025-margaux-cocoatree/scripts/cocoatree_with_groundtruth/results/cocoatree_gt/halabi/cocoatree_SCA_none.csv")
residue_swap = residue_swap.loc[~residue_swap["pdb_pos"].isna()]

orig_columns = [c for c in residue_swap.columns if c.startswith("orig_sector_")]
columns = [c for c in residue_swap.columns if c.startswith("sector_")]
columns.sort()

residue_lists = [
    residue_swap.loc[residue_swap[c], "pdb_named_pos"].to_numpy()
    for c in orig_columns]


for i, residue_pos in enumerate(residue_lists):
    color = colors[f"sector_{i+1}"]
    r, g, b = mcolors.to_rgb(color)

    for i in range(len(residue_pos)):
        residue = residue_pos[i]

        color_name = f"residue_{residue}"
        cmd.set_color(color_name, [r, g, b])
        cmd.color(color_name, f"resi {residue}")
        cmd.show("spheres", f"resi {residue}")
python end
