# PREFACE

**P**rioritization and **R**anking of **E**xoplanets **F**or **A**stronomical **C**haracterization and **E**xploration (**PREFACE**) is a Python package for selection of promising exoplanet transmission spectroscopy observations based on their expected scientific return and observational feasibility.

---

## To-do List
As this project is work in progress, here is the current to-do list:
- [ ] Finish documentation
- [ ] Fix a bug where the moonlight noise metric calculation sometimes return NaN
- [ ] (Maybe) Determine the default moonlight amplification factor that more properly and sensibly punishes full moon nights
- [ ] (Maybe) Reduce multiprocessing overhead
- [ ] (Maybe) Month- and location-dependent aerosol scattering parameters via end-to-end AERONET data retrieval

---

## Installation

Install the latest stable release from PyPI:
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

Full documentation (configuration reference, PREFACE workflow and output descriptions, and API) is available at **[preface.readthedocs.io](https://preface.readthedocs.io)**.

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
    moonlight_amplification_factor=5,
    toggle_graph_outputs=True,
    event_weight_graph_threshold=0.75
)

MultiprocessingConfigs = MultiprocessingConfigurations(
    toggle_multiprocessing=True,
    cores_to_leave_out=2
)

run_preface(
    TelescopeConfigurations=TelescopeConfigs,
    OutputConfigurations=OutputConfigs,
    MoonlightNoiseConfigurations=MoonlightConfigs,
    MultiprocessingConfigurations=MultiprocessingConfigs
)
```

---
## Authors

**Jake Staberg Morgan** (Original author)

**Chatdanai Sawangwong** (Current maintainer)\
email: chatdanai.saw@gmail.com

**Supachai Awiphan**\
email: supachai@narit.or.th

**Orarik Tasuya**\
email: orarik@narit.or.th

**Napaporn A-thano**\
email: napaporn@narit.or.th

<!-- ---

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
``` -->
