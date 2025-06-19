from matplotlib.colors import LinearSegmentedColormap

# #4ab38d
# #154c3b
lightcolors = {
    "cocoatree": "#AB0000",
    "other": "#006c4b",
    "both": "#006ab2",
    "default": "#9cadbc",
    "sector_2": "#AB0000",
    "sector_3": "#006ab2",
    "sector_1": "#006c4b",
    "sector_4": "#845EC2",
    "sector_5": "yellow",
    "sector_6": "orange"}


sectors_cm = {
    "sector_1": LinearSegmentedColormap.from_list(
        "sector_1", ['white', lightcolors["sector_1"]]),
    "sector_2": LinearSegmentedColormap.from_list(
        "sector_2", ['white', lightcolors["sector_2"]]),
    "sector_3": LinearSegmentedColormap.from_list(
        "sector_3", ['white', lightcolors["sector_3"]]),
    "sector_4": LinearSegmentedColormap.from_list(
        "sector_4", ['white', lightcolors["sector_4"]]),

    }

darkcolors = {
    "cocoatree": "#31775d",
    "other": "#796392",
    "both": "#2a719a",
    "default": "#9cadbc"
}

colors = lightcolors
