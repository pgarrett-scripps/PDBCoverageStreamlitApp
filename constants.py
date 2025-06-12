import os

# Create a list of available color maps
COLOR_MAPS = list({"coolwarm" ,"rainbow", "viridis", "plasma", "inferno", "magma", "cividis", "twilight", "twilight_shifted",
              "turbo", "nipy_spectral", "gist_ncar", "gist_rainbow", "hsv", "flag", "prism", "ocean", "gist_stern",
              "gnuplot", "gnuplot2", "CMRmap", "cubehelix", "brg", "gist_earth", "terrain", "gist_heat", "hot",
              "afmhot", "copper", "pink", "spring", "autumn", "cool", "winter", "bwr", "seismic", "bone", "cividis",
              "twilight", "twilight_shifted", "turbo", "nipy_spectral", "gist_ncar", "gist_rainbow", "hsv", "flag",
              "prism", "ocean", "gist_stern", "gnuplot", "gnuplot2", "CMRmap", "cubehelix", "brg", "gist_earth",
              "terrain", "gist_heat", "hot", "afmhot", "copper", "pink", "spring", "autumn", "cool", "winter", "bwr",
              "seismic", "bone"})
DEFAULT_COLOR_MAP = 'coolwarm'


DEFAULT_PROTEIN_ID = 'P60174'
_DEFAULT_PEPTIDES = ['KFFVGGNWK/2', 'KFFVGGNWK/3', 'KQSLGELIGTLNAAK', 'QSLGELIGTLNAAK', 'VPADTEVVCAPPTAYIDFAR',
                    'VPADTEVVC[+20]APPTAYIDFAR', '[Acetyl]-VPADTEVVCAPPTAYIDFAR', 'VPADTEVVCAPPTAYIDFAR', 'IAVAAQNCYK',
                    'IAVAAQNCYK', 'IAVAAQNCYK', 'VTNGAFTGEISPGMIK', 'DCGATWVVLGHSER', 'DCGATWVVLGHSER',
                    'RHVFGESDELIGQK', 'HVFGESDELIGQK', 'HVFGESDELIGQK', 'HVFGESDELIGQK', 'HVFGESDELIGQK',
                    'VAHALAEGLGVIACIGEK', 'VAHALAEGLGVIACIGEK', 'VAHALAEGLGVIACIGEK', 'VVLAYEPVWAIGTGK',
                    'VVLAYEPVWAIGTGK', 'VVLAYEPVWAIGTGK', 'VVLAYEPVWAIGTGK', 'TATPQQAQEVHEK', 'TATPQQAQEVHEK',
                    'TATPQQAQEVHEK', 'SNVSDAVAQSTR', 'SNVSDAVAQSTR', 'SNVSDAVAQSTR', 'IIYGGSVTGATCK', 'IIYGGSVTGATCK',
                    'IIYGGSVTGATCK', 'ELASQPDVDGFLVGGASLKPEFVDIINAK', 'ELASQPDVDGFLVGGASLKPEFVDIINAK',
                    'ELASQPDVDGFLVGGASLKPEFVDIINAK']

DEFAULT_PEPTIDES = '\n'.join(peptide for peptide in _DEFAULT_PEPTIDES)

def get_env_str(var_name, default):
    return os.getenv(var_name, default)


PDB_APP_URL = get_env_str('PDB_APP_URL', 'https://pdb-cov.streamlit.app/')


DEFAULT_PROTEIN_SEQUENCE = 'MAPSRKFFVGGNWKMNGRKQSLGELIGTLNAAKVPADTEVVCAPPTAYIDFARQKLDPKIAVAAQNCYKVTNGAFTGEISPGMIKDCGATWVVLGHSERRHVFGESDELIGQKVAHALAEGLGVIACIGEKLDEREAGITEKVVFEQTKVIADNVKDWSKVVLAYEPVWAIGTGKTATPQQAQEVHEKLRGWLKSNVSDAVAQSTRIIYGGSVTGATCKELASQPDVDGFLVGGASLKPEFVDIINAKQ'