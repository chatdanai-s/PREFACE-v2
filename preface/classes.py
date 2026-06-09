import os
import sys
from datetime import timedelta, datetime
import numpy as np
import pandas as pd
import joblib

# Checks if variable is a number
def isNumber(var):
    return isinstance(var, (int, float)) and not np.isnan(var)


# Telescope configs
class TelescopeConfigurations:
    def __init__(
        self,
        instrument: str,
        filter_name: str,
        run_mode: str,
        toggle_sky_noise=True,
        toggle_defocus=True,
    ):  
        # Store configuration keys

        # Input check
        InputErrorFlag = False

# Output configs
class OutputConfigurations:
    def __init__(
        self,
        observation_start: datetime,
        observation_end: datetime,
        output_path: str,
        metric_mode="Rank",
        viable_cumulative_cut=0.97,
        toggle_graph_outputs=True,
        event_weight_graph_threshold=0.5,
    ):
        # Store configuration keys
        
        # Input check
        InputErrorFlag = False


# Moonlight noise configs (Default: False)
class MoonlightNoiseConfigurations:
    def __init__(
        self,
        toggle_moonlight_noise=False,
        scattering_aod="Default",
        absorption_aod="Default",
        asymmetry_factor="Default",
        moonlight_amplification_factor=5,
    ):
        # Store configuration keys
        
        # Input check
        InputErrorFlag = False

# Multiprocessing configs (Default: True, use all cores except one)
class MultiprocessingConfigurations:
    def __init__(
        self,
        toggle_multiprocessing=True,
        cores_to_leave_out=1,
    ):
        # Store configuration keys
        
        # Input check
        InputErrorFlag = False

