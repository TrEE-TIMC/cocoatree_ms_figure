import numpy as np
from matplotlib.lines import Line2D
import matplotlib.pyplot as plt
from PIL import Image

from .colors_and_labels import colors
from .colors_and_labels import sectors_cm, labels


def plot_scatter_sectors(ax, results, ica, icb, annotate=True,
                         add_labels=True):
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
        zorder=1,
        color=colors["default"],)

    sectors_columns = [c for c in results.columns if c.startswith("sector")]

    for sector in sectors_columns + ["default"]:
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
            zorder=45,
            markerfacecolor="none",
            )

        mask = results["is_only_orig"]
        ax.plot(
            results.loc[mask & marker_mask, ica],
            results.loc[mask & marker_mask, icb],
            marker="x",
            linewidth=0,
            markersize=4,
            zorder=40,
            color=color,
            )

        mask = results["is_both"]
        ax.plot(
            results.loc[mask & marker_mask, ica],
            results.loc[mask & marker_mask, icb],
            marker="o",
            linewidth=0,
            markersize=4,
            zorder=50,
            color=color)

        ax.tick_params(
            axis='both', which='both', labelsize='xx-small',
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

    if add_labels:
        ax.set_xlabel(ica, fontsize="x-small",
                      fontweight="bold",
                      labelpad=2)
        ax.set_ylabel(icb, fontsize="x-small", fontweight="bold", labelpad=2)


def create_legend():
    sector_colors = ["Red", "Green", "Blue", "Purple", "Crimson"]
    sectors = [
        Line2D([0], [0], linewidth=0, marker='o', color=colors[f"sector_{i}"],
               label=f"{sector_colors[i-1]} sector",
               markersize=4) for i in range(1, 5)]
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
                 columns=["sector_1", "sector_2", "sector_3"], title=""):
    """
    Plot sectors

    Parameters
    ----------
    """
    for col in columns:
        i = int(col.split("_")[-1])
        sec = results[col].values[np.newaxis, :].astype(float)
        sec[sec == 0] = np.nan
        ax.matshow(sec, aspect="auto", cmap=sectors_cm[f"sector_{i}"],
                   vmin=0, vmax=1.2)

    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_yticks([0])
    ax.set_yticklabels([title], fontweight="bold")

    ax.tick_params(axis='both', which='both', labelsize='x-small',
                   bottom=False, top=False, labeltop=False, labelbottom=False,
                   left=False, right=False, labelleft=False, labelright=True)


def plot_coev_mat_sectors(fig, ax, results, coev_mat):
    """
    Plot coevolution matrix of sectors

    Parameters
    ----------
    ax : matplotlib.Axes object

    results :

    coev_mat :
    """

    # Get the number of sectors
    # /!\ depends on the structure of the results file
    sector_columns = [
        col for col in results.columns if col.startswith("sector_")]
    sector_columns.sort()
    num_sectors = len(sector_columns)
    sectors_list = []
    for sect in sector_columns:
        sect_id = sect.split("_")[-1]
        weights = results.loc[results[sect], f"IC{sect_id}"]
        weights = weights.sort_values(ascending=False)
        sect_pos = results.loc[weights.index, 'filtered_msa_pos']
        sectors_list.append(sect_pos.astype(int).values)

    sector_sizes = [len(sec) for sec in sectors_list]
    cumul_sizes = sum(sector_sizes)
    sorted_pos = np.concat(sectors_list)

    submatrix = coev_mat.loc[sorted_pos, sorted_pos].values
    submatrix[np.diag_indices_from(submatrix)] = np.nan

    # Get extent for the colorbars
    vmax = min(2, np.nanmax(submatrix))
    vmin = min(0, np.nanmin(submatrix))

    # Plot coevolution matrix
    im = ax.imshow(submatrix,
                   vmin=vmin, vmax=vmax,
                   interpolation='none', aspect='equal',
                   origin="lower",
                   extent=[0, cumul_sizes, 0, cumul_sizes], cmap='RdBu_r')
    cb = fig.colorbar(im)
    cb.ax.tick_params(labelsize="x-small")
    cb.ax.locator = plt.MaxNLocator(nbins=3)
    cb.update_ticks()
    cb.set_label("Coevolution metric", fontweight="bold")

    line_index = 0
    label_index = []
    for i in range(num_sectors):
        ax.axvline(line_index + sector_sizes[i],
                   color='0', linewidth=1)
        ax.axhline(line_index + sector_sizes[i],
                   color='0', linewidth=1)

        label_index += [np.sum(line_index) + sector_sizes[i] / 2]
        line_index += sector_sizes[i]

    ax.tick_params(axis='both', which='both', labelsize='small',
                   labelright=False, labelleft=True, right=False, left=True,
                   labeltop=False, labelbottom=True, bottom=True, top=False)

    ax.set_xticks(label_index)
    ax.set_yticks(label_index)
    ax.set_xticklabels(labels[:len(label_index)], fontweight="bold")
    ax.set_yticklabels(labels[:len(label_index)], fontweight="bold")

    return cb


def plot_image(ax, path, trim=100):
    """
    Plot an image
    """
    image = np.asarray(
        Image.open(path))

    image = image[trim:-trim, trim:-trim]

    [ax.spines[s].set_linewidth(0) for s in ax.spines]
    ax.set_xticks([])
    ax.set_yticks([])

    ax.imshow(image)
