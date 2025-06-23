import numpy as np


# Now create some form of confusion matrix
def compute_IOU_metric(set1, set2, metric="intersection"):
    union = len(set(set1).union(set2))
    intersection = len(set(set1).intersection(set2))
    if metric == "IOU":
        return intersection / union
    else:
        return intersection



def get_best_ordered_sectors(res, dataset="halabi", type="SCA"):
    sector_cols = [c for c in res.columns if c.startswith("sector")]
    orig_sector_cols = [c for c in res.columns if c.startswith("orig_sector")]

    # Get original sectors in the original order
    orig_sectors = [np.where(res[c])[0] for c in orig_sector_cols]

    sectors = []
    # It's not the best strategy, but just get iteratively the "best
    # sector"

    all_scores = np.zeros((len(sector_cols), len(sector_cols)))
    for i, s in enumerate(orig_sectors):
        for j, s1 in enumerate(sector_cols):
            all_scores[i, j] = compute_IOU_metric(
                    s,
                    np.where(res[s1])[0])
    order = all_scores.argmax(axis=1)
    if len(np.unique(order)) != len(sector_cols):
        order = np.arange(len(sector_cols))
        if dataset == "halabi" and type == "NMI":
            order = [2, 1, 0]
        if dataset == "halabi" and type == "MI":
            order = [1, 0, 2]
        if dataset == "rhomboid" and type == "MI":
            order = [0, 2, 1]

    return order

