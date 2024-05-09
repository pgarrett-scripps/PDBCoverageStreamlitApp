import os

AAS = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V',
       'W', 'X', 'Y', 'Z']
AA_TO_THREE_LETTER_CODE = {
    "A": "Ala",
    "R": "Arg",
    "N": "Asn",
    "D": "Asp",
    "C": "Cys",
    "E": "Glu",
    "Q": "Gln",
    "G": "Gly",
    "H": "His",
    "I": "Ile",
    "L": "Leu",
    "K": "Lys",
    "M": "Met",
    "F": "Phe",
    "P": "Pro",
    "S": "Ser",
    "T": "Thr",
    "W": "Trp",
    "Y": "Tyr",
    "V": "Val",
    "B": "Asx",
    "Z": "Glx",
    "X": "Xaa",
    "U": "Sec",
    "O": "Pyl",
    "J": "Xle",
}

# Create a list of available color maps
color_maps = ["coolwarm" ,"rainbow", "viridis", "plasma", "inferno", "magma", "cividis", "twilight", "twilight_shifted",
              "turbo", "nipy_spectral", "gist_ncar", "gist_rainbow", "hsv", "flag", "prism", "ocean", "gist_stern",
              "gnuplot", "gnuplot2", "CMRmap", "cubehelix", "brg", "gist_earth", "terrain", "gist_heat", "hot",
              "afmhot", "copper", "pink", "spring", "autumn", "cool", "winter", "bwr", "seismic", "bone", "cividis",
              "twilight", "twilight_shifted", "turbo", "nipy_spectral", "gist_ncar", "gist_rainbow", "hsv", "flag",
              "prism", "ocean", "gist_stern", "gnuplot", "gnuplot2", "CMRmap", "cubehelix", "brg", "gist_earth",
              "terrain", "gist_heat", "hot", "afmhot", "copper", "pink", "spring", "autumn", "cool", "winter", "bwr",
              "seismic", "bone"]


DEFAULT_PROTEIN_ID = 'P60174'
DEFAULT_COVERAGE_ARRAY = [0, 0, 0, 0, 0, 2, 2, 2, 2, 2, 2, 2, 2, 2, 0, 0, 0, 0, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
                          2, 2, 4, 4,
                          4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 0, 0, 0, 0, 0, 0, 3, 3, 3, 3, 3, 3, 3,
                          3, 3, 3, 1,
                          1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 5,
                          5, 5, 5, 5,
                          5, 5, 5, 5, 5, 5, 5, 5, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 0, 0, 0, 0, 0,
                          0, 0, 0, 0,
                          0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
                          4, 4, 4, 4,
                          3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 0, 0, 0, 0, 0, 0, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
                          3, 3, 3, 3,
                          3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
                          3, 3, 3, 3,
                          3, 3, 3, 0]
DEFAULT_PEPTIDES = ['KFFVGGNWK', 'KFFVGGNWK', 'KQSLGELIGTLNAAK', 'QSLGELIGTLNAAK', 'VPADTEVVCAPPTAYIDFAR',
                    'VPADTEVVCAPPTAYIDFAR', 'VPADTEVVCAPPTAYIDFAR', 'VPADTEVVCAPPTAYIDFAR', 'IAVAAQNCYK',
                    'IAVAAQNCYK', 'IAVAAQNCYK', 'VTNGAFTGEISPGMIK', 'DCGATWVVLGHSER', 'DCGATWVVLGHSER',
                    'RHVFGESDELIGQK', 'HVFGESDELIGQK', 'HVFGESDELIGQK', 'HVFGESDELIGQK', 'HVFGESDELIGQK',
                    'VAHALAEGLGVIACIGEK', 'VAHALAEGLGVIACIGEK', 'VAHALAEGLGVIACIGEK', 'VVLAYEPVWAIGTGK',
                    'VVLAYEPVWAIGTGK', 'VVLAYEPVWAIGTGK', 'VVLAYEPVWAIGTGK', 'TATPQQAQEVHEK', 'TATPQQAQEVHEK',
                    'TATPQQAQEVHEK', 'SNVSDAVAQSTR', 'SNVSDAVAQSTR', 'SNVSDAVAQSTR', 'IIYGGSVTGATCK', 'IIYGGSVTGATCK',
                    'IIYGGSVTGATCK', 'ELASQPDVDGFLVGGASLKPEFVDIINAK', 'ELASQPDVDGFLVGGASLKPEFVDIINAK',
                    'ELASQPDVDGFLVGGASLKPEFVDIINAK']


def get_env_str(var_name, default):
    return os.getenv(var_name, default)


PDB_APP_URL = get_env_str('PDB_APP_URL', 'http://localhost:8502/')
