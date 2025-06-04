

def annotate_results(results):
    results = results.loc[~results["filtered_msa_pos"].isna()]

    sca_sectors = [c for c in results.columns if c.startswith("sector")]
    orig_sectors = [c for c in results.columns if c.startswith("orig_")]

    results["is_cocoatree_sector"] = results[sca_sectors].sum(
        axis=1) > 0
    results["is_orig_sector"] = results[orig_sectors].sum(
        axis=1) > 0
    results["is_both"] = (
        results["is_cocoatree_sector"] & results["is_orig_sector"])
    results["is_only_cocoatree"] = (
        results["is_cocoatree_sector"] & ~results["is_orig_sector"])
    results["is_only_orig"] = (
        ~results["is_cocoatree_sector"] & results["is_orig_sector"])
    return results
