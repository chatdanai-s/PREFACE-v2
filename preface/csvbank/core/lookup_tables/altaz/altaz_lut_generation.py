import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from preface import CSV_core_folder

# Creates lookup tables (LUT) containing local AltAz positions of the Sun and Moon
# First create a table with 5-min precision, then interpolate for 1-min precision.