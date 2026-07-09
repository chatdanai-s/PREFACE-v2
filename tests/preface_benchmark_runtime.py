import sys
import datetime as dt
from pathlib import Path
import time
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from preface import run_preface, wipe_intermediate_csvs
from preface.configs import TelescopeConfigurations, OutputConfigurations, MoonlightNoiseConfigurations, MultiprocessingConfigurations

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------

OUTPUT_FOLDER = Path(r"C:\Users\WBS\Documents\PREFACE_benchmarks\time_benchmark")
OUTPUT_FOLDER.mkdir(parents=True, exist_ok=True)

ObsStart = dt.datetime(2027,1,1)

TOTAL_CPU_CORES = 12
CORES_TO_TEST = [10]           # Actual cores used
MONTHS_TO_TEST = [1, 2, 3, 4]

results = []

# -----------------------------------------------------------------------------
# Benchmark
# -----------------------------------------------------------------------------

for used_cores in CORES_TO_TEST:

    cores_to_leave_out = TOTAL_CPU_CORES - used_cores

    for months in MONTHS_TO_TEST:
        ObsEnd = ObsStart + pd.DateOffset(months=months)

        print("=" * 70)
        print(f"Running benchmark: {used_cores} core(s), {months} month(s)")
        print("=" * 70)

        TelescopeConfigs = TelescopeConfigurations(
            instrument="TNT ULTRASPEC",
            filter_name="r",
            run_mode="Half_Well",
            toggle_sky_noise=True,
            toggle_defocus=False,
        )

        OutputConfigs = OutputConfigurations(
            observation_start=ObsStart,
            observation_end=ObsEnd,
            output_folder=str(OUTPUT_FOLDER),
            metric_mode="Rank",
            viable_cumulative_cut=0.97,
            toggle_graph_outputs=True,
            event_weight_graph_threshold=0.75,
        )

        MoonlightConfigs = MoonlightNoiseConfigurations(
            toggle_moonlight_noise=True,
            scattering_aod=0.2,
            absorption_aod=0.3,
            asymmetry_factor=0.6,
            moonlight_amplification_factor=10,
        )

        MultiprocessingConfigs = MultiprocessingConfigurations(
            toggle_multiprocessing=True,
            cores_to_leave_out=cores_to_leave_out,
        )

        start = time.perf_counter()

        run_preface(
            TelescopeConfigurations=TelescopeConfigs,
            OutputConfigurations=OutputConfigs,
            MoonlightNoiseConfigurations=MoonlightConfigs,
            MultiprocessingConfigurations=MultiprocessingConfigs,
        )

        elapsed = time.perf_counter() - start

        wipe_intermediate_csvs()

        print(f"Finished in {elapsed/60:.2f} min\n")

        results.append(
            {
                "cores_used": used_cores,
                "cores_left_out": cores_to_leave_out,
                "months": months,
                "observation_start": ObsStart,
                "observation_end": ObsEnd,
                "runtime_seconds": elapsed,
                "runtime_minutes": elapsed / 60,
            }
        )

        # Save after every run in case benchmarking is interrupted
        pd.DataFrame(results).to_csv(
            OUTPUT_FOLDER / "runtime_benchmark.csv",
            index=False,
        )

print("\nAll benchmarks completed.")