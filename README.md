# PREFACE

**P**rioritization and **R**anking of **E**xoplanets **F**or **A**stronomical **C**haracterization and **E**xploration (**PREFACE**) is a Python package for selection of promising exoplanet transmission spectroscopy observations based on their expected scientific return and observational feasibility.

PREFACE combines telescope characteristics, target properties, observing windows, detector effects, airmass, and optional moonlight modelling to produce numerically ranked observing opportunities for transit follow-up campaigns.

---

## Installation

<!-- Install the latest stable release from PyPI: -->
You will be able to install the latest stable release from PyPI:

```bash
pip install preface-spearnet
```

---

## Usage

Using `preface` consists of four steps:

1. Configure the observing instrument with `TelescopeConfigurations`.
2. Define the observing window and output options with `OutputConfigurations`.
3. Optionally configure moonlight modelling and multiprocessing with `MoonlightNoiseConfigurations` and `MultiprocessingConfigurations`.
4. Execute the complete pipeline with `run_preface()`.

Input validation is performed automatically before pipeline execution.

### Example

```python
import datetime as dt

from preface import run_preface
from preface.configs import (
    TelescopeConfigurations,
    OutputConfigurations,
    MoonlightNoiseConfigurations,
    MultiprocessingConfigurations,
)

ObsStart = dt.datetime(2025, 10, 1)
ObsEnd = dt.datetime(2026, 5, 31)
OutputFolder = r"C:\PREFACE_Output"

TelescopeConfigs = TelescopeConfigurations(
    instrument="TNT ULTRASPEC",
    filter_name="r",
    run_mode="Half_Well",
    toggle_sky_noise=True,
    toggle_defocus=False
)

OutputConfigs = OutputConfigurations(
    observation_start=ObsStart,
    observation_end=ObsEnd,
    output_folder=OutputFolder,
    metric_mode="Rank",
    viable_cumulative_cut=0.97
)

MoonlightConfigs = MoonlightNoiseConfigurations(
    toggle_moonlight_noise=True,
    scattering_aod=0.2,
    absorption_aod=0.3,
    asymmetry_factor=0.6,
    moonlight_amplification_factor=5
)

MultiprocessingConfigs = MultiprocessingConfigurations(
    toggle_multiprocessing=True,
    cores_to_leave_out=2
)

run_preface(
    TelescopeConfigurations=TelescopeConfigs,
    OutputConfigurations=OutputConfigs,
    MoonlightnoiseConfigurations=MoonlightConfigs,
    MultiprocessingConfigurations=MultiprocessingConfigs
)
```

---

## Configuration Classes

### `TelescopeConfigurations`

Defines the observing setup.

| Parameter          | Type   | Description                                       | Valid Values / Restrictions                                 |
| ------------------ | ------ | ------------------------------------------------- | ----------------------------------------------------------- |
| `instrument`       | `str`  | Telescope or instrument name.                     | Must correspond to a supported instrument.                  |
| `filter_name`      | `str`  | Photometric filter to be used.                    | Must be valid for the selected instrument.                  |
| `run_mode`         | `str`  | Detector operating mode.                          | `"Half_Well"`, `"Spectral_Half_Well"`, or `"IR_Half_Well"`. |
| `toggle_sky_noise` | `bool` | Include sky background noise in pipeline.         | `True` or `False`.                                          |
| `toggle_defocus`   | `bool` | Include telescope defocus effects when supported. | `True` or `False`.                                          |

### `OutputConfigurations`

Defines observation timing and ranking behaviour.

| Parameter                      | Type       | Description                                | Valid Values / Restrictions      |
| ------------------------------ | ---------- | ------------------------------------------ | -------------------------------- |
| `observation_start`            | `datetime` | Start of observing period.                 | Must precede `observation_end` by at least one hour.  |
| `observation_end`              | `datetime` | End of observing period.                   | Must follow `observation_start` by at least one hour. |
| `output_folder`                | `str`      | Existing directory for generated outputs.  | Directory must already exist. |
| `metric_mode`                  | `str`      | Event ranking metric.                      | `"Rank"`, `"Habitable_Rank"`, `"Multi_Transit_Rank"`, or `"Multi_Transit_Habitable_Rank"`. |
| `viable_cumulative_cut`        | `float`    | Fraction of cumulative viability retained. | Between `0` and `1` (`0.97` recommended). |
| `toggle_graph_outputs`         | `bool`     | Produce diagnostic plots per transit.      | `True` or `False`. |
| `event_weight_graph_threshold` | `float`    | Minimum event weight for graph generation. | Non-negative (`0.50` to `0.99` recommended). |

### `MoonlightNoiseConfigurations`

Controls optional moonlight background modelling.

| Parameter                        | Type    | Description                                                 | Valid Values / Restrictions       |
| -------------------------------- | ------- | ----------------------------------------------------------- | --------------------------------- |
| `toggle_moonlight_noise`         | `bool`  | Enable moonlight noise calculations.                        | `True` or `False`.                |
| `scattering_aod`                 | `float` | Atmospheric scattering aerosol optical depth.               | Non-negative.                     |
| `absorption_aod`                 | `float` | Atmospheric absorption aerosol optical depth.               | Non-negative.                     |
| `asymmetry_factor`               | `float` | Scattering asymmetry parameter.                             | Between `-1` and `1`.             |
| `moonlight_amplification_factor` | `float` | Empirical scaling controlling modeled moonlight influence.  | Positive value ≥ `5` recommended. |

### `MultiprocessingConfigurations`

Controls CPU utilization during multiprocessing in Phase Two.

| Parameter                | Type   | Description                                            | Valid Values / Restrictions |
| ------------------------ | ------ | ------------------------------------------------------ | --------------------------- |
| `toggle_multiprocessing` | `bool` | Enable multiprocessing.                                | `True` or `False`.          |
| `cores_to_leave_out`     | `int`  | Number of logical CPU cores reserved from computation. | Integer ≥ `0`.              |


### Configuration Notes
* Some run modes and defocus calculations are instrument-specific. The program will stop and report accordingly if an unsupported combination is detected.
* Moonlight modelling currently supports only UBVRI and ugriz photometric filters, on top of telescope-valid filters.
* By default, multiprocessing reserves one logical CPU core to maintain system responsiveness. If no cores are reserved (0 cores left), the program prompts for confirmation before continuing, as using all logical cores may cause the system to become unresponsive.

---

## Built-in Utilities

The package also contains several convenience utilities:

```python
import preface

# List supported telescopes
preface.telescope_list

# Retrieve supported filters for an instrument
preface.get_available_filters_list("TNT ULTRASPEC")

# Open bundled telescope reference table
preface.open_scope_csv()

# Remove temporary intermediate files to clear up space
preface.wipe_intermediate_csvs()
```

---

## Outputs

PREFACE writes all generated files to the specified `output_folder`. The output directory is organized according to two processing stages:

```text
output_folder/
├── phase_1/
│   ├── nonviable_target_list/
│   └── viable_target_list/
└── phase_2/
    ├── cumulative_observability_scores/
    ├── full_ranked_event_list/
    ├── graphs/
    └── individual_planets/
```

### Phase One (`phase_1`)

Selects the most promising planetary candidates for transmission spectroscopy from catalogues, by ranking them proportionally to expected data signal-to-noise ratio.

| Directory                | Contents                                                                                                   |
| ------------------------ | ---------------------------------------------------------------------------------------------------------- |
| `nonviable_target_list/` | Non-viable targets with metric values below the input cumulative ranking cutoff `viable_cumulative_cut`.   |
| `viable_target_list/`    | Promising targets retained after applying the cumulative ranking cutoff and carried forward for Phase Two. |


### Phase Two (`phase_2`)

Determines the best exoplanet observation dates from predicted transits observed at selected telescope location.

| Directory                          | Contents                                                      |
| ---------------------------------- | ------------------------------------------------------------- |
| `cumulative_observability_scores/` | Observability score accumulated from all transits per planet. |
| `full_ranked_event_list/`          | Complete ranked list of individual observable transit events. |
| `graphs/`                          | Diagnostic plots for each transit (when enabled).             |
| `individual_planets/`              | Per-planet computed transits and each of their scores.        |

---
<!-- 
## Authors

**Name 1**
Affilation 1

**Name 2**
Affilation 2

**Name 3**
Affilation 3

---

## Citation

If you use **PREFACE** in academic work, please cite the associated publications:

```bibtex
@software{NameYear_PREFACE,
  author    = {Authors},
  title     = {Title},
  year      = {Year},
  publisher = {Publisher}
}
```

```bibtex
@phdthesis{Morgan2020_SPEARNET,
  author       = {Morgan, Jake S.},
  title        = {SPEARNET: A Pilot Exoplanet Transmission Spectroscopy Survey},
  school       = {The University of Manchester},
  year         = {2020},
  type         = {PhD thesis},
  date         = {2020-11-18},
  month        = nov,
  url          = {https://research.manchester.ac.uk/en/studentTheses/spearnet-a-pilot-exoplanet-transmission-spectroscopy-survey}
}
```
-->
