import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from preface import wipe_intermediate_csvs

wipe_intermediate_csvs()
