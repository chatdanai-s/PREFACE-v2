import sys
import datetime as dt
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from preface import run_preface
from preface.configs import TelescopeConfigurations, OutputConfigurations, MoonlightNoiseConfigurations, MultiprocessingConfigurations


ObsStart = dt.datetime(2025,10,1)
ObsEnd = dt.datetime(2026,1,1)
outputFolder = rf'C:\Users\WBS\Documents\PREFACE_test_output'  # Change as needed in local device

TelescopeConfigs = TelescopeConfigurations(
    instrument='TNT ULTRASPEC',
    filter_name='r',
    run_mode='Half_Well',
    toggle_sky_noise=True,
    toggle_defocus=False
)
OutputConfigs = OutputConfigurations(
    observation_start=ObsStart,
    observation_end=ObsEnd,
    output_folder=outputFolder,
    metric_mode='Rank',
    viable_cumulative_cut=0.90,
    toggle_graph_outputs=True,
    event_weight_graph_threshold=0.75
)
MoonlightConfigs = MoonlightNoiseConfigurations(
    toggle_moonlight_noise=True,
    scattering_aod=0.2,
    absorption_aod=0.3,
    asymmetry_factor=0.6,
    moonlight_amplification_factor=10
)
MultiprocessingConfigs = MultiprocessingConfigurations(
    toggle_multiprocessing=True,
    cores_to_leave_out=2
)

run_preface(TelescopeConfigurations=TelescopeConfigs,
            OutputConfigurations=OutputConfigs,
            MoonlightNoiseConfigurations=MoonlightConfigs,
            MultiprocessingConfigurations=MultiprocessingConfigs)

