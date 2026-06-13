from .configs import TelescopeConfigurations, OutputConfigurations, MoonlightNoiseConfigurations, MultiprocessingConfigurations
from . import P1_ModCheck
from . import P1_ImpactMerger
from . import P1_ExoplanetseuImpactMerger
from . import P1_WorkingTEPSetBuilder
from . import P1_RankMaker
from . import P1_Cutter
from . import P1_ViabilitySplitter
# from . import P2_MultiprocessingWrapper
# from . import P2_PostCleaner

import pandas as pd
import numpy as np
import os
import subprocess
import sys
from pathlib import Path

# Locate CSV files
PACKAGE_ROOT = Path(__file__).resolve().parent
CSV_core_folder = PACKAGE_ROOT / "csvbank" / "core"
CSV_intermediate_folder = PACKAGE_ROOT / "csvbank" / "intermediate"

# Scope dataframe load and useful variables
scope_df = pd.read_csv(CSV_core_folder / 'Scope.csv')
telescope_list = scope_df['Telescope'].tolist()
csvbank_location = PACKAGE_ROOT / "csvbank"

# Useful utility functions
def wipe_intermediate_csvs():
    for path in CSV_intermediate_folder.rglob("*"):
        if path.is_file() or path.is_symlink():
            path.unlink()
    print('All intermediate CSV files in PREFACE successfully wiped from system.')


def get_available_filters_list(instrument: str):
    # Find the row corresponding to the instrument
    row = scope_df.loc[scope_df["Telescope"] == instrument]
    if row.empty:
        raise ValueError(f"Instrument '{instrument}' not found.")

    row = row.iloc[0]

    # Find filters that have both mzp_ and msky_ columns defined as numbers
    mzp_filters = {
        col[len("mzp_"):]
        for col in scope_df.columns
        if col.startswith("mzp_") and pd.notna(row[col])
    }
    msky_filters = {
        col[len("msky_"):]
        for col in scope_df.columns
        if col.startswith("msky_") and pd.notna(row[col])
    }
    
    # Return filters that exist in both
    available_filters = sorted(list(mzp_filters & msky_filters))
    return available_filters


def open_scope_csv():
    csv_path = CSV_core_folder / "Scope.csv"
    
    if sys.platform.startswith("win"):
        os.startfile(csv_path)
    elif sys.platform == "darwin":  # macOS
        subprocess.run(["open", str(csv_path)], check=True)
    else:  # Linux
        subprocess.run(["xdg-open", str(csv_path)], check=True)


##### Main preface function #####
def run_preface(TelescopeConfigurations: TelescopeConfigurations,
                OutputConfigurations: OutputConfigurations,
                MoonlightnoiseConfigurations: MoonlightNoiseConfigurations,
                MultiprocessingConfigurations: MultiprocessingConfigurations):
    
    # Unpack variables (Some unused variables left for code readability)
    instrument, filter_name, run_mode, toggle_sky_noise, toggle_defocus = TelescopeConfigurations.unpack
    obs_start, obs_end, output_folder, metric_mode, viable_cumulative_cut, toggle_graph_outputs, event_weight_graph_threshold = OutputConfigurations.unpack
    toggle_moonlight_noise, scattering_aod, absorption_aod, asymmetry_factor, moonlight_amplification_factor = MoonlightnoiseConfigurations.unpack
    toggle_multiprocessing, total_cores, cores_to_leave_out, cores_used = MultiprocessingConfigurations.unpack

    scope_idx = np.where(scope_df['Telescope'] == instrument)[0][0]

    # Filter availability check if moonlight noise enabled
    if toggle_moonlight_noise == True:
        available_filters = ['U','B','V','R','I','u','g','r','i','z']
        if filter_name not in available_filters:
            print("[InputCheck] Invalid filter_name for moonlight noise modeling detected. Is it UBVRI or ugriz?")
            sys.exit('[InputCheck] Invalid input(s) found -- PREFACE terminated.')


    # Run pipeline (Phase 1)
    print('\nRunning Phase One of PREFACE...\n')

    P1_ModCheck.Check(CSV_core_folder)
    P1_ImpactMerger.ExoOrgImpacts(CSV_core_folder)
    P1_ExoplanetseuImpactMerger.ExoeuImpacts(CSV_core_folder)
    P1_WorkingTEPSetBuilder.WorkBuilder(CSV_core_folder)

    P1_RankMaker.RankMaker(CSV_core_folder, CSV_intermediate_folder,
                           scope_df, scope_idx, *TelescopeConfigurations.unpack)

    RMin = P1_Cutter.RankMaker(CSV_core_folder, CSV_intermediate_folder,
                               scope_df, *TelescopeConfigurations.unpack, metric_mode, viable_cumulative_cut)
    
    P1_ViabilitySplitter.Splitter(CSV_core_folder, CSV_intermediate_folder, output_folder,
                                  *TelescopeConfigurations.unpack, metric_mode, viable_cumulative_cut, RMin)


    # Run pipeline (Phase 2)
    print('\nRunning Phase Two of PREFACE...\n')

    # P2_MultiprocessingWrapper.P2Wrap(CSV_core_folder, CSV_intermediate_folder, output_folder,
    #                                  scope_df, scope_idx, *TelescopeConfigurations.unpack, metric_mode, viable_cumulative_cut,
    #                                  *MoonlightnoiseConfigurations.unpack,
    #                                  toggle_graph_outputs, event_weight_graph_threshold,
    #                                  obs_start, obs_end,
    #                                  toggle_multiprocessing, cores_used)
    
    # P2_PostCleaner.Cleaner(CSV_core_folder, CSV_intermediate_folder, output_folder,
    #                        scope_df, scope_idx, *TelescopeConfigurations.unpack, metric_mode, viable_cumulative_cut,
    #                        obs_start, obs_end,
    #                        toggle_multiprocessing, cores_used)
    
    print('PREFACE run complete.')
