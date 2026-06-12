import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from preface import telescope_list, csvbank_location, get_available_filters_list, open_scope_csv

print(telescope_list)
print(csvbank_location)

instrument = 'VLT FORS2 (200kHz)'
available_filters = get_available_filters_list(instrument)
print(available_filters)

open_scope_csv()
