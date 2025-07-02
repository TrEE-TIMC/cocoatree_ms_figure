from matplotlib.colors import LinearSegmentedColormap

# #4ab38d
# #154c3b
lightcolors = {
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


sectors_cm = {
    "sector_1": LinearSegmentedColormap.from_list(
        "sector_1", ['white', lightcolors["sector_1"]]),
    "sector_2": LinearSegmentedColormap.from_list(
        "sector_2", ['white', lightcolors["sector_2"]]),
    "sector_3": LinearSegmentedColormap.from_list(
        "sector_3", ['white', lightcolors["sector_3"]]),
    "sector_4": LinearSegmentedColormap.from_list(
        "sector_4", ['white', lightcolors["sector_4"]]),
    "sector_5": LinearSegmentedColormap.from_list(
        "sector_5", ['white', lightcolors["sector_5"]]),
    "sector_6": LinearSegmentedColormap.from_list(
        "sector_6", ['white', lightcolors["sector_6"]]),
    "others": LinearSegmentedColormap.from_list(
        "other", ["white", "dimgray"]),

    }

darkcolors = {
    "cocoatree": "#31775d",
    "other": "#796392",
    "both": "#2a719a",
    "default": "#9cadbc"
}

colors = lightcolors
