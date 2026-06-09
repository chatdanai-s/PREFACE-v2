import os
import sys
from datetime import timedelta, datetime
import numpy as np
import pandas as pd
import joblib

from preface import scope_df


# Checks if variable is a number or boolean
def isNumber(var):
    return isinstance(var, (int, float)) and not np.isnan(var)
def isBoolean(var):
    return isinstance(var, bool)

# Checks InputErrorFlag, termminates if False
def checkInputErrorFlag(InputErrorFlag):
    if InputErrorFlag == True:
        sys.exit('[InputCheck] Invalid input(s) found -- PREFACE terminated.')


# Telescope configs
class TelescopeConfigurations:
    def __init__(
        self,
        instrument: str,
        filter_name: str,
        run_mode: str,
        toggle_sky_noise=True,
        toggle_defocus=True
    ):
        
        # Check instrument and filter_name
        InputErrorFlag = False
        try:
            scope_idx = np.where(scope_df['Telescope'] == instrument)[0][0]

            mzp = scope_df[f'mzp_{filter_name}'].iloc[scope_idx]
            msky = scope_df[f'msky_{filter_name}'].iloc[scope_idx]
            if not (isNumber(mzp) and isNumber(msky)):
                InputErrorFlag = True
                print(f"[InputCheck] filter_name {filter_name} is not available for this telescope!")
        
        except IndexError:
            scope_idx = None
            InputErrorFlag = True
            print(f"[InputCheck] {instrument} is not a valid telescope name! "
                  "See preface.telescope_list for list of available telescopes.")
        except KeyError:
            InputErrorFlag = True
            print(f"[InputCheck] {filter_name} is not a valid filter_name!")

        # Check run_mode
        if run_mode not in ['Half_Well', 'Spectral_Half_Well', 'IR_Half_Well']:
            InputErrorFlag = True
            print(f"[InputCheck] {run_mode} is not a valid run_mode! "
                  "Is it 'Half_Well', 'Spectral_Half_Well', or 'IR_Half_Well'?")

        if scope_idx is not None:
            Lambda_Cent = scope_df['Lambda_Cent'].iloc[scope_idx]
            Dispersion = scope_df['Dispersion'].iloc[scope_idx]
            if (run_mode in ['Spectral_Half_Well', 'IR_Half_Well']) and not (isNumber(Lambda_Cent) and isNumber(Dispersion)):
                InputErrorFlag = True
                print(f"[InputCheck] run_mode {run_mode} is not available for this telescope! "
                      "Set run_mode = 'Half_Well'.")
        elif scope_idx == None:
            if run_mode in ['Spectral_Half_Well', 'IR_Half_Well']:
                print("[InputCheck] run_mode validity cannot be checked due to invalid telescope input.")

        # Check toggle sky noise
        if not isBoolean(toggle_sky_noise):
            InputErrorFlag = True
            print(f"[InputCheck] toggle_sky_noise {toggle_sky_noise} is not a boolean!")

        # Check toggle defocus
        if not isBoolean(toggle_defocus):
            InputErrorFlag = True
            print(f"[InputCheck] toggle_defocus {toggle_defocus} is not a boolean!")

        if scope_idx is not None:
            theta_DF = scope_df['Theta_DF'].iloc[scope_idx]
            if toggle_defocus and not isNumber(theta_DF):
                InputErrorFlag = True
                print("[InputCheck] Defocus is not available for this telescope! Set toggle_defocus to False.")
        elif scope_idx == None:
            if toggle_defocus:
                print("[InputCheck] Defocus validity cannot be checked due to invalid telescope index.")

        # Terminate if issues arise
        checkInputErrorFlag(InputErrorFlag)

        # Store valid telescope configuration keys
        self.instrument = instrument
        self.filter_name = filter_name
        self.run_mode = run_mode
        self.toggle_sky_noise = toggle_sky_noise
        self.toggle_defocus = toggle_defocus


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
        
        # Store valid output keys
        self.observation_start = observation_start
        self.observation_end = observation_end
        self.output_path = output_path
        self.metric_mode = metric_mode
        self.viable_cumulative_cut = viable_cumulative_cut
        self.toggle_graph_outputs = toggle_graph_outputs
        self.event_weight_graph_threshold = event_weight_graph_threshold


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
        # Store valid moonlight noise keys
        self.toggle_moonlight_noise = toggle_moonlight_noise
        self.scattering_aod = scattering_aod
        self.absorption_aod = absorption_aod
        self.asymmetry_factor = asymmetry_factor
        self.moonlight_amplification_factor = moonlight_amplification_factor


# Multiprocessing configs (Default: True, use all cores except one)
class MultiprocessingConfigurations:
    def __init__(
        self,
        toggle_multiprocessing=True,
        cores_to_leave_out=1,
    ):
        # Check multiprocessing configurations
        InputErrorFlag = False

        if not isBoolean(toggle_multiprocessing):
            InputErrorFlag = True
            print(f"[InputCheck] toggle_multiprocessing {toggle_multiprocessing} is not a boolean!")

        if not (
            isNumber(cores_to_leave_out)
            and int(cores_to_leave_out) == cores_to_leave_out
            and cores_to_leave_out >= 0
        ):
            InputErrorFlag = True
            print(f"[InputCheck] cores_to_leave_out {cores_to_leave_out} is not a non-negative integer!")

        cpu_count = joblib.cpu_count()
        ncores = cpu_count - int(cores_to_leave_out)

        if (ncores <= 0) or (ncores > cpu_count):
            InputErrorFlag = True
            print(f"[InputCheck] Invalid core configuration for multiprocessing. "
                  f"({ncores}/{cpu_count} cores to be used)")
        elif ncores == cpu_count:
            print("[InputCheck] WARNING: All CPU cores will be used. Your machine may become unresponsive.")
            if input(' '*13 + "Continue anyways? [Y/N]: ").upper() != "Y":
                sys.exit("[InputCheck] PREFACE terminated.")

        # Terminate if issues arise
        checkInputErrorFlag(InputErrorFlag)

        # Store valid multiprocessing keys
        self.toggle_multiprocessing = toggle_multiprocessing
        self.cores_to_leave_out = cores_to_leave_out
        self.total_cores = cpu_count
        self.cores_used = ncores

