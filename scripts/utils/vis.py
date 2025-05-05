import numpy as np
from .colors import colors


def plot_scatter_sectors(ax, results, ica, icb, annotate=True):
    """
    Plot sectors as a scatter plot

    Parameters
    ----------

    ax : matplotlib.Axes object

    ica : name of the first axis

    icb : name of the second axis
    """
    mask = ~(results["is_cocoatree_sector"] | results["is_orig_sector"])
    ax.plot(
        results.loc[mask, ica],
        results.loc[mask, icb],
        marker="o",
        linewidth=0,
        color=colors["default"],)

    mask = results["is_only_cocoatree"]
    ax.plot(
        results.loc[mask, ica],
        results.loc[mask, icb],
        marker="o",
        linewidth=0,
        color=colors["cocoatree"],
        label="Only cocoatree",
        )

    mask = results["is_only_orig"]
    ax.plot(
        results.loc[mask, ica],
        results.loc[mask, icb],
        marker="o",
        linewidth=0,
        label="Only orig.",
        color=colors["other"])

    mask = results["is_both"]
    ax.plot(
        results.loc[mask, ica],
        results.loc[mask, icb],
        marker="o",
        linewidth=0,
        label="Both",
        color=colors["both"])

    if annotate:
        mask = results["is_cocoatree_sector"] | results["is_orig_sector"]
        for i, (x, y, text) in results.loc[mask,
                                           [ica,
                                            icb,
                                            "pdb_named_pos"]].iterrows():
            if isinstance(text, float):
                if np.isnan(text):
                    text = "."
                else:
                    text = "%d" % text
            ax.text(x, y+.01, text,
                    fontsize=8,
                    ha='center', va='center')

    ax.spines["top"].set_linewidth(0)
    ax.spines["right"].set_linewidth(0)

    ax.set_xlabel(ica, fontweight="bold")
    ax.set_ylabel(icb, fontweight="bold")
