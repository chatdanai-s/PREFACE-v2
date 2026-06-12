import sys
import datetime as dt
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from preface import run_preface
from preface.configs import TelescopeConfigurations, OutputConfigurations, MoonlightNoiseConfigurations, MultiprocessingConfigurations


ObsStart = dt.datetime(2025,10,1,0,0,0)
ObsEnd = dt.datetime(2026,5,31,0,0,0)
ThisFolder = Path(__file__).resolve().parent

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
    output_folder=ThisFolder,
    metric_mode='Rank',
    viable_cumulative_cut=0.97,
    toggle_graph_outputs=True,
    event_weight_graph_threshold=0.5
)
MoonlightConfigs = MoonlightNoiseConfigurations(
    toggle_moonlight_noise=False,
    scattering_aod="Default",
    absorption_aod="Default",
    asymmetry_factor="Default",
    moonlight_amplification_factor=5
)
MultiprocessingConfigs = MultiprocessingConfigurations(
    toggle_multiprocessing=True,
    cores_to_leave_out=1
)

run_preface(TelescopeConfigurations=TelescopeConfigs,
            OutputConfigurations=OutputConfigs,
            MoonlightnoiseConfigurations=MoonlightConfigs,
            MultiprocessingConfigurations=MultiprocessingConfigs)
