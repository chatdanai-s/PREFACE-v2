from .classes import TelescopeConfigurations, OutputConfigurations, MoonlightNoiseConfigurations, MultiprocessingConfigurations
import P0_InputCheck, P0_ModCheck, P0_ImpactMerger, P0_ExoplanetseuImpactMerger
import P1_WorkingTEPSetBuilder, P1_RankMaker, P1_Cutter, P1_ViabilitySplitter
import P2_MultiprocessingWrapper, P2_PostCleaner

import pandas as pd
import numpy as np

CSV_core_folder = None
CSV_intermediate_folder = None

scope_df = pd.read_csv(rf'{CSV_core_folder}/Scope.csv')
telescope_list = scope_df['Telescope'].tolist()

def list_available_filters(instrument):
    pass

def run_preface(TelescopeConfigurations,
                OutputConfigurations,
                MoonlightnoiseConfigurations,
                MultiprocessingConfigurations):
    
    # Unpack variables

    # Run pipeline

    pass