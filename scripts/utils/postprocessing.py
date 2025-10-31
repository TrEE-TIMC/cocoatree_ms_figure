import warnings


def annotate_results(results):
    warnings.simplefilter("ignore")
    results = results.loc[~results["filtered_msa_pos"].isna()]

    sca_xcors = [c for c in results.columns if c.startswith("xcor")]
    orig_sectors = [c for c in results.columns if c.startswith("orig_")]

    results.loc[:, "is_cocoatree_sector"] = results[sca_xcors].sum(
        axis=1) > 0

    results.loc[:, "is_orig_sector"] = results[orig_sectors].sum(
        axis=1) > 0
    results.loc[:, "is_both"] = (
        results["is_cocoatree_sector"] & results["is_orig_sector"])
    results.loc[:, "is_only_cocoatree"] = (
        results["is_cocoatree_sector"] & ~results["is_orig_sector"])
    results.loc[:, "is_only_orig"] = (
        ~results["is_cocoatree_sector"] & results["is_orig_sector"])
    return results
