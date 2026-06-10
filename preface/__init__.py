from .configs import TelescopeConfigurations, OutputConfigurations, MoonlightNoiseConfigurations, MultiprocessingConfigurations
import P1_ModCheck, P1_ImpactMerger, P1_ExoplanetseuImpactMerger
import P1_WorkingTEPSetBuilder, P1_RankMaker, P1_Cutter, P1_ViabilitySplitter
import P2_MultiprocessingWrapper, P2_PostCleaner

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


# def list_available_filters(instrument): # TBA (More pressing things are due first)
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
    P1_ModCheck.Check(CSV_core_folder)
    P1_ImpactMerger.ExoOrgImpacts(CSV_core_folder)
    P1_ExoplanetseuImpactMerger.ExoeuImpacts(CSV_core_folder)

    # Above is finished, below is unfinished.

    P1_WorkingTEPSetBuilder.WorkBuilder(CSV_core_folder)

    # P1_RankMaker.RankMaker(CSV_core_folder, CSV_intermediate_folder,
    #                        scope_df, scope_idx, *TelescopeConfigurations.unpack)
    
    # min_rank = P1_Cutter.RankMaker(CSV_core_folder, CSV_intermediate_folder, output_folder,
    #                            scope_df, *TelescopeConfigurations.unpack, metric_mode, viable_cumulative_cut)
    
    # P1_ViabilitySplitter.Splitter(CSV_core_folder, CSV_intermediate_folder, 
    #                               *TelescopeConfigurations.unpack, metric_mode, viable_cumulative_cut, min_rank)


    # # Run pipeline (Phase 2)
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
