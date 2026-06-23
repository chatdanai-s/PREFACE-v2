Quickstart
==========

Using PREFACE involves four steps:

1. Configure the observing instrument with :class:`~preface.configs.TelescopeConfigurations`.
2. Define the observing window and output options with :class:`~preface.configs.OutputConfigurations`.
3. Optionally configure moonlight modelling and multiprocessing with
   :class:`~preface.configs.MoonlightNoiseConfigurations` and
   :class:`~preface.configs.MultiprocessingConfigurations`.
4. Execute the complete pipeline with :func:`~preface.run_preface`.

Input validation is performed automatically before pipeline execution.

Example
-------

The following example runs PREFACE over an observing window from October 2025
to May 2026, using the TNT ULTRASPEC instrument with moonlight modelling and
multiprocessing enabled:

.. code-block:: python

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
        toggle_graph_outputs=True,
        event_weight_graph_threshold=0.75
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
        MoonlightNoiseConfigurations=MoonlightConfigs,
        MultiprocessingConfigurations=MultiprocessingConfigs
    )

