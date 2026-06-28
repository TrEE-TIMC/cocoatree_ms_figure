import numpy as np


# Now create some form of confusion matrix
def compute_IOU_metric(set1, set2, metric="intersection"):
    union = len(set(set1).union(set2))
    intersection = len(set(set1).intersection(set2))
    if metric == "IOU":
        return intersection / union
    else:
        return intersection


def compute_IOU_metric_all_vs_all(res):
    sector_cols = [c for c in res.columns if c.startswith("xcor")]
    sector_cols.sort()
    orig_sector_cols = [c for c in res.columns if c.startswith("orig_sector")]
    orig_sector_cols.sort()

    # Get original sectors in the original order
    orig_sectors = [np.where(res[c])[0] for c in orig_sector_cols]

    all_scores = np.zeros((len(sector_cols), len(sector_cols)))
    for i, s in enumerate(orig_sectors):
        for j, s1 in enumerate(sector_cols):
            all_scores[i, j] = compute_IOU_metric(
                    s,
                    np.where(res[s1])[0])
    return all_scores


def compute_all_vs_all(sectors_per_method, comp1=0, comp2=0):
    IOU_metric = np.zeros(
        (len(sectors_per_method), len(sectors_per_method)))
    for i, met1 in enumerate(sectors_per_method):
        for j, met2 in enumerate(sectors_per_method):
            IOU_metric[i, j] = compute_IOU_metric(met1[comp1], met2[comp2])
    return IOU_metric


def get_best_ordered_sectors(res, dataset="halabi", type="SCA",
                             correction=None):
    sector_cols = [c for c in res.columns if c.startswith("xcor")]

    all_scores = compute_IOU_metric_all_vs_all(res)

    order = all_scores.argmax(axis=1)

    if len(np.unique(order)) != len(sector_cols):
        order = np.arange(len(sector_cols))

    if dataset == "halabi" and type == "NMI":
        order = [2, 1, 0]
    elif dataset == "halabi" and type == "MI" and correction is None:
        order = [1, 0, 2]
    elif dataset == "halabi" and type == "MI" and correction == "APC":
        order = [0, 1, 2]

    elif dataset == "rhomboid" and type == "MI" and correction is None:
        order = [1, 2, 0]
    elif dataset == "rhomboid" and type == "MI" and correction == "APC":
        order = [1, 0, 2]
    elif dataset == "DHFR" and type == "MI" and correction == "APC":
        order = [3, 1, 2, 0]
    elif dataset == "DHFR" and type == "NMI":
        order = [0, 2, 3, 1]
    elif dataset == "DHFR" and type == "MI" and correction is None:
        order = [2, 3, 0, 1]

    return order


