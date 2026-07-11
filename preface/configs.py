import os
import sys
from datetime import timedelta, datetime
import numpy as np
import pandas as pd
import joblib
from pathlib import Path


# Scope df (Cannot import from __init__.py due to circular import)
PACKAGE_ROOT = Path(__file__).resolve().parent
CSV_core_folder = PACKAGE_ROOT / "csvbank" / "core"
scope_df = pd.read_csv(CSV_core_folder / 'Scope.csv')

# Checks if variable is a number
def isNumber(var):
    return isinstance(var, (int, float)) and not np.isnan(var)
# Checks if variable is a boolean
def isBoolean(var):
    return isinstance(var, bool)

# Checks InputErrorFlag (bool), terminates if False
def checkInputErrorFlag(InputErrorFlag):
    if InputErrorFlag == True:
        sys.exit('[InputCheck] Invalid input(s) found -- PREFACE terminated.')


# Telescope configs
class TelescopeConfigurations:
    """
    Configuration container for telescope and instrument settings.
    """

    def __init__(
        self,
        instrument: str,
        filter_name: str,
        run_mode: str,
        toggle_sky_noise=True,
        toggle_defocus=False
    ):
        """
        Initialize and validate telescope-specific configuration parameters
        before storing them for use throughout the PREFACE pipeline.

        Parameters
        ----------
        instrument : str
            Name of the observing instrument or telescope.
        filter_name : str
            Photometric filter to be used for calculations.
        run_mode : str
            Observation mode. Must be one of ``"Half_Well"``,
            ``"Spectral_Half_Well"``, or ``"IR_Half_Well"``.
        toggle_sky_noise : bool, default=True
            Whether to include sky background noise in calculations.
        toggle_defocus : bool, default=False
            Whether to apply telescope defocus modelling when available.

        Returns
        -------
        None

        Raises
        ------
        SystemExit
            If any supplied configuration parameter is invalid.
        """

        # Check instrument and filter_name
        InputErrorFlag = False
        scope_idx = None
        
        try:
            scope_idx = np.where(scope_df['Telescope'] == instrument)[0][0]

            mzp = scope_df[f'mzp_{filter_name}'].iloc[scope_idx]
            msky = scope_df[f'msky_{filter_name}'].iloc[scope_idx]
            if not (isNumber(mzp) and isNumber(msky)):
                InputErrorFlag = True
                print(f"[InputCheck] filter_name {filter_name} is not available for this telescope! "
                       "See preface.get_available_filters_list(instrument) for list of available filters for telescope")
        
        except IndexError:
            InputErrorFlag = True
            print(f"[InputCheck] {instrument} is not a valid telescope name! "
                  "See preface.telescope_list for list of available telescopes.")
            
        except KeyError:
            InputErrorFlag = True
            print(f"[InputCheck] {filter_name} is not a valid filter_name! "
                  "See preface.get_available_filters_list(instrument) for list of available filters for telescope")


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
            if (toggle_defocus == True) and not isNumber(theta_DF):
                InputErrorFlag = True
                print("[InputCheck] Defocus is not available for this telescope! Set toggle_defocus to False.")
        elif scope_idx == None:
            if toggle_defocus:
                print("[InputCheck] Defocus validity cannot be checked due to invalid telescope input.")


        # Terminate if issues arise
        checkInputErrorFlag(InputErrorFlag)


        # Store valid telescope configuration keys
        self.instrument = instrument
        self.filter_name = filter_name
        self.run_mode = run_mode
        self.toggle_sky_noise = toggle_sky_noise
        self.toggle_defocus = toggle_defocus

    @property
    def unpack(self):
        """
        Return the stored telescope configuration as an ordered list.

        Returns
        -------
        list
            Values in the order:
            ``[instrument, filter_name, run_mode,
            toggle_sky_noise, toggle_defocus]``.
        """

        return [
            self.instrument,
            self.filter_name,
            self.run_mode,
            self.toggle_sky_noise,
            self.toggle_defocus
        ]


# Output configs
class OutputConfigurations:
    """
    Configuration container controlling PREFACE output generation
    and exoplanet ranking settings.
    """

    def __init__(
        self,
        observation_start: datetime,
        observation_end: datetime,
        output_folder: str | Path,
        metric_mode="Rank",
        viable_cumulative_cut=0.97,
        toggle_graph_outputs=True,
        event_weight_graph_threshold=0.75,
    ):
        """
        Initialize and validate output-related configuration parameters
        before storing them for use throughout the PREFACE pipeline.

        Parameters
        ----------
        observation_start : datetime.datetime
            Beginning of the observation interval.
        observation_end : datetime.datetime
            End of the observation interval.
        output_folder : str or pathlib.Path
            Existing directory in which output files will be written.
        metric_mode : str, default="Rank"
            Ranking metric used for event prioritization.
            Must be one of ``"Rank"``, ``"Habitable_Rank"``,
            ``"Multi_Transit_Rank"``, or ``"Multi_Transit_Habitable_Rank"``.
        viable_cumulative_cut : float, default=0.97
            Cumulative viability threshold between 0 and 1.
        toggle_graph_outputs : bool, default=True
            Whether diagnostic plots for each transit event should be generated.
        event_weight_graph_threshold : float, default=0.75
            Minimum event weight required for graph generation.
            Must be between 0 and 1.

        Returns
        -------
        None

        Raises
        ------
        SystemExit
            If any supplied configuration parameter is invalid.
        """

        # Check observation times
        InputErrorFlag = False

        if observation_end <= observation_start:
            InputErrorFlag = True
            print('[InputCheck] Invalid Observation times. observation_end is earlier than observation_start!')
        elif observation_end - observation_start < timedelta(hours=1):
            InputErrorFlag = True
            print('[InputCheck] Invalid Observation times. Observation period must be at least one hour!')


        # Check if output folder exists
        if (not isinstance(output_folder, (str, Path))) or (not os.path.isdir(output_folder)):
            InputErrorFlag = True
            print(f"[InputCheck] {output_folder} is not a valid output_folder! "
                  "Does the directory exist?")


        # Check metric_mode
        if metric_mode not in ['Rank', 'Habitable_Rank', 'Multi_Transit_Rank', 'Multi_Transit_Habitable_Rank']:
            InputErrorFlag = True
            print(f"[InputCheck] {metric_mode} is not a valid metric_mode! "
                  "Is it 'Rank', 'Habitable_Rank', 'Multi_Transit_Rank', or 'Multi_Transit_Habitable_Rank'?") 


        # Check viable_cumulative_cut
        if (not isNumber(viable_cumulative_cut)) or (not (0 < viable_cumulative_cut <= 1)):
            InputErrorFlag = True
            print(f"[InputCheck] {viable_cumulative_cut} is not a valid viable_cumulative_cut! "
                  "Is it a number between 0 and 1?") 


        # Check toggle_graph_outputs
        if not isBoolean(toggle_graph_outputs):
            InputErrorFlag = True
            print(f"[InputCheck] toggle_multiprocessing {toggle_graph_outputs} is not a boolean!")

        
        # Check event_weight_graph_threshold
        if (not isNumber(event_weight_graph_threshold)) or (not (0 < event_weight_graph_threshold <= 1)):
            InputErrorFlag = True
            print(f"[InputCheck] {event_weight_graph_threshold} is not a valid event_weight_graph_threshold! "
                  "Is it a number between 0 and 1 (except 0)?")


        # Terminate if issues arise
        checkInputErrorFlag(InputErrorFlag)


        # Store valid output keys
        self.observation_start = observation_start
        self.observation_end = observation_end
        self.output_folder = Path(output_folder)
        self.metric_mode = metric_mode
        self.viable_cumulative_cut = viable_cumulative_cut
        self.toggle_graph_outputs = toggle_graph_outputs
        self.event_weight_graph_threshold = event_weight_graph_threshold

    @property
    def unpack(self):
        """
        Return the stored output configuration as an ordered list.

        Returns
        -------
        list
            Values in the order:
            ``[observation_start, observation_end, output_folder,
            metric_mode, viable_cumulative_cut,
            toggle_graph_outputs, event_weight_graph_threshold]``.
        """

        return [
            self.observation_start,
            self.observation_end,
            self.output_folder,
            self.metric_mode,
            self.viable_cumulative_cut,
            self.toggle_graph_outputs,
            self.event_weight_graph_threshold,
        ]


# Moonlight noise configs (Default: False)
class MoonlightNoiseConfigurations:
    """
    Configuration container for moonlight background noise modelling.
    """

    def __init__(
        self,
        toggle_moonlight_noise=False,
        scattering_aod=0.2,
        absorption_aod=0.3,
        asymmetry_factor=0.6,
        moonlight_amplification_factor=10,
    ):
        """
        Initialize and validate moonlight noise model parameters
        before storing them for use throughout the PREFACE pipeline.

        Parameters
        ----------
        toggle_moonlight_noise : bool, default=False
            Whether moonlight noise modelling is enabled.
        scattering_aod : float, default=0.2
            Atmospheric scattering aerosol optical depth.
            Must be at least 0.
        absorption_aod : float, default=0.3
            Atmospheric absorption aerosol optical depth.
            Must be at least 0.
        asymmetry_factor : float, default=0.6
            Scattering phase-function asymmetry parameter.
            Must be between -1 and +1.
        moonlight_amplification_factor : float, default=10
            Empirical scaling factor applied to the moonlight model.
            Defined as the effective decrease in lunar magnitude used to
            amplify modeled moonlight brightness and its impact on target SNR.

        Returns
        -------
        None

        Raises
        ------
        SystemExit
            If any supplied configuration parameter is invalid.
        """
        
        # Check toggle_moonlight_noise
        InputErrorFlag = False

        if not isBoolean(toggle_moonlight_noise):
            InputErrorFlag = True
            print(f"[InputCheck] toggle_moonlight_noise {toggle_moonlight_noise} is not a boolean!")

        # If False, then do not bother
        if toggle_moonlight_noise == True:
            # Check scattering_aod
            if (not isNumber(scattering_aod) or scattering_aod < 0):
                InputErrorFlag = True
                print(f"[InputCheck] '{scattering_aod}' is not a valid scattering_aod! "
                      "Is it a number >= 0?")

            # Check absorption_aod
            if (not isNumber(absorption_aod) or absorption_aod < 0):
                InputErrorFlag = True
                print(f"[InputCheck] '{absorption_aod}' is not a valid scattering_aod! "
                      "Is it a number >= 0?")

            # Check asymmetry_factor
            if (not isNumber(asymmetry_factor) or not (-1 < asymmetry_factor < 1)):
                InputErrorFlag = True
                print(f"[InputCheck] '{asymmetry_factor}' is not a valid asymmetry_factor! "
                      "Is it a number between -1 and +1? (0.5-0.8 recommended)")

            # Check moonlight_amplification_factor
            if not isNumber(moonlight_amplification_factor):
                InputErrorFlag = True
                print(f"[InputCheck] '{moonlight_amplification_factor}' is not a valid moonlight_amplification_factor! "
                      "Is it a number? (default=10)")


        # Terminate if issues arise
        checkInputErrorFlag(InputErrorFlag)


        # Store valid moonlight noise keys
        self.toggle_moonlight_noise = toggle_moonlight_noise
        self.scattering_aod = scattering_aod
        self.absorption_aod = absorption_aod
        self.asymmetry_factor = asymmetry_factor
        self.moonlight_amplification_factor = moonlight_amplification_factor
      
    @property
    def unpack(self):
        """
        Return the stored moonlight noise configuration as an ordered list.

        Returns
        -------
        list
            Values in the order:
            ``[toggle_moonlight_noise, scattering_aod,
            absorption_aod, asymmetry_factor,
            moonlight_amplification_factor]``.
        """

        return [
            self.toggle_moonlight_noise,
            self.scattering_aod,
            self.absorption_aod,
            self.asymmetry_factor,
            self.moonlight_amplification_factor,
        ]


# Multiprocessing configs (Default: True, use all cores except one)
class MultiprocessingConfigurations:
    """
    Configuration container controlling multiprocessing behaviour
    and CPU resource allocation.
    """

    def __init__(
        self,
        toggle_multiprocessing=True,
        cores_to_leave_out=1,
    ):
        """
        Initialize and validate multiprocessing configuration parameters
        before storing them for use throughout the PREFACE pipeline.

        Parameters
        ----------
        toggle_multiprocessing : bool, default=True
            Whether multiprocessing should be enabled.
        cores_to_leave_out : int, default=1
            Number of logical CPU cores reserved from computation.

        Returns
        -------
        None

        Raises
        ------
        SystemExit
            If the multiprocessing configuration is invalid or the user
            declines confirmation when all CPU cores would be utilized.
        """

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
        self.total_cores = cpu_count
        self.cores_to_leave_out = cores_to_leave_out
        self.cores_used = ncores

    @property
    def unpack(self):
        """
        Return the stored multiprocessing configuration as an ordered list.

        Returns
        -------
        list
            Values in the order:
            ``[toggle_multiprocessing, total_cores,
            cores_to_leave_out, cores_used]``.
        """

        return [
            self.toggle_multiprocessing,
            self.total_cores,
            self.cores_to_leave_out,
            self.cores_used
        ]
