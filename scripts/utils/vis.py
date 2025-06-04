import numpy as np
from .colors import colors
from matplotlib.lines import Line2D
import matplotlib.pyplot as plt
from .colors import sectors_cm


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
        marker=".",
        linewidth=0,
        markersize=4,
        color=colors["default"],)

    for sector in ["default", "sector_1", "sector_2", "sector_3"]:
        if sector != "default":
            marker_mask = results[sector]
        else:
            marker_mask = ~results["is_cocoatree_sector"]
        color = colors[sector]
        mask = results["is_only_cocoatree"]
        ax.plot(
            results.loc[mask & marker_mask, ica],
            results.loc[mask & marker_mask, icb],
            marker="o",
            markeredgecolor=color,
            linewidth=0,
            markersize=4,
            markerfacecolor="none",
            )

        mask = results["is_only_orig"]
        ax.plot(
            results.loc[mask & marker_mask, ica],
            results.loc[mask & marker_mask, icb],
            marker="x",
            linewidth=0,
            markersize=4,
            color=color,
            )

        mask = results["is_both"]
        ax.plot(
            results.loc[mask & marker_mask, ica],
            results.loc[mask & marker_mask, icb],
            marker="o",
            linewidth=0,
            markersize=4,
            color=color)

        ax.tick_params(
            axis='both', which='both', labelsize='x-small',
            bottom=True, top=False, labeltop=False, labelbottom=True,
            left=True, right=False, labelleft=True, labelright=False)
        ax.xaxis.set_major_locator(plt.MaxNLocator(3))
        ax.yaxis.set_major_locator(plt.MaxNLocator(3))

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

    ax.set_xlabel(ica, fontsize="small",
                  fontweight="bold",
                  labelpad=2)
    ax.set_ylabel(icb, fontsize="small", fontweight="bold", labelpad=2)


def create_legend():
    sectors = [
        Line2D([0], [0], linewidth=0, marker='o', color=colors[f"sector_{i}"],
               label=f"sector {i}",
               markersize=4) for i in range(1, 4)]
    methods = [
        Line2D([0], [0], marker='o', color="0",
               linewidth=0,
               label="both",
               markersize=4),
        Line2D([0], [0], marker='o', color="none",
               markeredgecolor="0",
               linewidth=0,
               label="cocoatree",
               markersize=4),
        Line2D([0], [0], marker='x', color="0",
               label="orig.", linewidth=0,
               markersize=4),
        Line2D([0], [0], marker='.', color="0",
               label="none", linewidth=0,
               markersize=4),

               ]
    return {"sectors": sectors, "methods": methods}


def plot_sectors(ax, results,
                 columns=["sector_1", "sector_2", "sector_3"]):
    """
    Plot sectors

    Parameters
    ----------
    """
    for i, col in enumerate(columns):
        sec = results[col].values[np.newaxis, :].astype(float)
        sec[sec == 0] = np.nan
        ax.matshow(sec, aspect="auto", cmap=sectors_cm[f"sector_{i+1}"],
                   vmin=0)

    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_yticks([0])
    ax.set_yticklabels(["cocoatree"], fontweight="bold")

    ax.tick_params(axis='both', which='both', labelsize='x-small',
                   bottom=False, top=False, labeltop=False, labelbottom=False,
                   left=False, right=False, labelleft=False, labelright=True)
