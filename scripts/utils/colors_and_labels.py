from matplotlib.colors import LinearSegmentedColormap

# #4ab38d
# #154c3b
lightcolors = {
    "cocoatree": "#AB0000",
    "other": "#006c4b",
    "both": "#006ab2",
    "default": "#9cadbc",
    "sector_1": "#AB0000",  # Red
    "sector_3": "#0023B2",  # "#006ab2", # Green
    "sector_2": "#078400",  # "#006c4b", # Blue
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

labels = ["Red", "Green", "Blue", "Purple", "Crimson", "Orange"]


halabi_cmap_from_paper = {
    'vertebrate': 'black',
    'not vertebrate': 'grey',
    'fungi': '#ffe4b5',
    'bacteria': '#afdde9',
    'chymotrypsin': 'darkblue',
    'trypsin': 'magenta',
    'tryptase': 'yellow',
    'kallikrein': 'darkorange',
    'granzyme': 'lime'
    }

halabi_cmap = {
    'vertebrate': '#798e87',
    'not vertebrate': '#c27d38',
    'fungi': '#ccc591',
    'bacteria': '#29211f',
    'chymotrypsin': '#ff0000',
    'trypsin': '#00a08a',
    'tryptase': '#f2ad00',
    'kallikrein': '#f98400',
    'granzyme': '#5bbcd6',
    'Mammalia': '#ff9898',
    'Actinopterygii': '#d9636c',
    'Amphibia': '#a91e45',
    'Malacostraca': '#691238',
    'other': 'lightgrey',
    'Insecta': '#251714'
}
