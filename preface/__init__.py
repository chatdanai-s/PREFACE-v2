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
import sys
from pathlib import Path

# Locate CSV files
PACKAGE_ROOT = Path(__file__).resolve().parent
CSV_core_folder = PACKAGE_ROOT / "csvbank" / "core"
CSV_intermediate_folder = PACKAGE_ROOT / "csvbank" / "intermediate"

# Scope dataframe load
scope_df = pd.read_csv(CSV_core_folder / 'Scope.csv')
telescope_list = scope_df['Telescope'].tolist()

# Useful CSV Bank functions
def locate_csvbanks():
    print(PACKAGE_ROOT)

def clear_intermediate_csvs():
    for path in CSV_intermediate_folder.rglob("*"):
        if path.is_file() or path.is_symlink():
            path.unlink()
    print('All intermediate CSV files in pipeline removed from system.')

# TBA functions (More pressing things are due first)
# def list_available_filters(instrument): 
#     pass


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
    print('\nRunning Phase 1 of PREFACE...\n')

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


    # # Run pipeline (Phase 2)
    # print('\nRunning Phase 2 of preface...\n')

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
