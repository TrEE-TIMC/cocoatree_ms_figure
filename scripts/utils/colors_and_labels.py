from matplotlib.colors import LinearSegmentedColormap

# #4ab38d
# #154c3b
lightcolors = {
    "cocoatree": "#AB0000",
    "other": "#006c4b",
    "both": "#006ab2",
    "default": "#9cadbc",
    "xcor_1": "#AB0000",  # Red
    "xcor_3": "#0023B2",  # "#006ab2", # Green
    "xcor_2": "#078400",  # "#006c4b", # Blue
    "xcor_4": "#845EC2",
    "xcor_5": "crimson",
    "xcor_6": "orange"}


sectors_cm = {
    "sector_1": LinearSegmentedColormap.from_list(
        "sector_1", ['white', lightcolors["xcor_1"]]),
    "sector_2": LinearSegmentedColormap.from_list(
        "sector_2", ['white', lightcolors["xcor_2"]]),
    "sector_3": LinearSegmentedColormap.from_list(
        "sector_3", ['white', lightcolors["xcor_3"]]),
    "sector_4": LinearSegmentedColormap.from_list(
        "sector_4", ['white', lightcolors["xcor_4"]]),
    "sector_5": LinearSegmentedColormap.from_list(
        "sector_5", ['white', lightcolors["xcor_5"]]),
    "sector_6": LinearSegmentedColormap.from_list(
        "sector_6", ['white', lightcolors["xcor_6"]]),
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

IC_labels = ["IC_Red", "IC_Green", "IC_Blue", "IC_Purple", "IC_Crimson", "IC_Orange"]


halabi_cmap_from_paper = {
    'vertebrate': 'black',
    'invertebrate': 'grey',
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
    'invertebrate': '#c27d38',
    'fungi': '#ccc591',
    'bacteria': '#29211f',
    'Chymotrypsin': '#00a08a',
    'Trypsin': '#ff0000',
    'Tryptase': '#f2ad00',
    'Kallikrein': '#f98400',
    'Granzyme': '#5bbcd6',
    'Mammalia': '#ff9898',
    'Actinopterygii': '#d9636c',
    'Amphibia': '#a91e45',
    'Malacostraca': '#691238',
    'other': 'lightgrey',
    'Insecta': '#251714'
}

halabi_longer_cmap = {
    'vertebrate': '#798e87',
    'invertebrate': '#c27d38',
    'fungi': '#ccc591',
    'bacteria': '#29211f',
    'Chymotrypsin': '#00a08a',
    'Trypsin': '#ff0000',
    'Tryptase': '#f2ad00',
    'Kallikrein': '#f98400',
    'Granzyme': '#5bbcd6',
    'Mammalia': '#C969A1',
    'Actinopterygii': '#CE4441',
    'Amphibia': '#EE8577',
    'Malacostraca': '#EB7926',
    'other': 'lightgrey',
    'Insecta': '#FFBB44',
    'Actinobacteria': '#859B6C',
    'Arachnida': '#62929A',
    'Oligochaeta': '#004F63'
}
