Configuration
=============

PREFACE uses four configuration classes to control pipeline behaviour. Each is
passed as a keyword argument to :func:`~preface.run_preface`.

.. contents:: Configuration Classes
   :local:
   :depth: 1

----

TelescopeConfigurations
-----------------------

Defines the observing instrument setup.

.. code-block:: python

    from preface.configs import TelescopeConfigurations

    TelescopeConfigs = TelescopeConfigurations(
        instrument="TNT ULTRASPEC",
        filter_name="r",
        run_mode="Half_Well",
        toggle_sky_noise=True,
        toggle_defocus=False
    )

.. list-table::
   :header-rows: 1
   :widths: 25 10 35 20 10

   * - Parameter
     - Type
     - Description
     - Valid Values
     - Default
   * - ``instrument``
     - ``str``
     - Telescope or instrument name.
     - Must correspond to a supported instrument. See ``preface.telescope_list``.
     - *required*
   * - ``filter_name``
     - ``str``
     - Photometric filter to be used.
     - Must be valid for the selected instrument. See :func:`~preface.get_available_filters_list`.
     - *required*
   * - ``run_mode``
     - ``str``
     - Detector operating mode.
     - ``"Half_Well"``, ``"Spectral_Half_Well"``, or ``"IR_Half_Well"``.
     - *required*
   * - ``toggle_sky_noise``
     - ``bool``
     - Include sky background noise in pipeline.
     - ``True`` or ``False``.
     - ``True``
   * - ``toggle_defocus``
     - ``bool``
     - Include telescope defocus effects when supported.
     - ``True`` or ``False``.
     - ``False``

.. note::
   Some run modes and defocus calculations are instrument-specific. The pipeline
   will stop and report accordingly if an unsupported combination is detected.

----

OutputConfigurations
--------------------

Defines the observation timing window, output directory, and ranking behaviour.

.. code-block:: python

    from preface.configs import OutputConfigurations
    import datetime as dt

    OutputConfigs = OutputConfigurations(
        observation_start=dt.datetime(2025, 10, 1),
        observation_end=dt.datetime(2026, 5, 31),
        output_folder=r"C:\PREFACE_Output",
        metric_mode="Rank",
        viable_cumulative_cut=0.90,
        toggle_graph_outputs=True,
        event_weight_graph_threshold=0.75
    )

.. list-table::
   :header-rows: 1
   :widths: 30 12 35 13 10

   * - Parameter
     - Type
     - Description
     - Valid Values
     - Default
   * - ``observation_start``
     - ``datetime``
     - Start of observing period.
     - Must precede ``observation_end`` by at least one hour.
     - *required*
   * - ``observation_end``
     - ``datetime``
     - End of observing period.
     - Must follow ``observation_start`` by at least one hour.
     - *required*
   * - ``output_folder``
     - ``str``
     - Existing directory for generated outputs.
     - Directory must already exist.
     - *required*
   * - ``metric_mode``
     - ``str``
     - Event ranking metric.
     - ``"Rank"``, ``"Habitable_Rank"``, ``"Multi_Transit_Rank"``, or ``"Multi_Transit_Habitable_Rank"``.
     - *required*
   * - ``viable_cumulative_cut``
     - ``float``
     - Fraction of cumulative viability retained.
     - Between ``0`` and ``1`` (``0.90`` to ``0.97`` recommended).
     - ``0.90``
   * - ``toggle_graph_outputs``
     - ``bool``
     - Produce diagnostic plots per transit.
     - ``True`` or ``False``.
     - ``True``
   * - ``event_weight_graph_threshold``
     - ``float``
     - Minimum event weight for graph generation.
     - Non-negative (``0.50`` to ``0.99`` recommended).
     - ``0.75``

----

MoonlightNoiseConfigurations
-----------------------------

Controls optional moonlight background modelling.

.. code-block:: python

    from preface.configs import MoonlightNoiseConfigurations

    MoonlightConfigs = MoonlightNoiseConfigurations(
        toggle_moonlight_noise=True,
        scattering_aod=0.2,
        absorption_aod=0.3,
        asymmetry_factor=0.6,
        moonlight_amplification_factor=5
    )

.. list-table::
   :header-rows: 1
   :widths: 35 10 35 10 10

   * - Parameter
     - Type
     - Description
     - Valid Values
     - Default
   * - ``toggle_moonlight_noise``
     - ``bool``
     - Enable moonlight noise calculations.
     - ``True`` or ``False``.
     - ``False``
   * - ``scattering_aod``
     - ``float``
     - Atmospheric scattering aerosol optical depth.
     - Non-negative.
     - ``0.2``
   * - ``absorption_aod``
     - ``float``
     - Atmospheric absorption aerosol optical depth.
     - Non-negative.
     - ``0.3``
   * - ``asymmetry_factor``
     - ``float``
     - Scattering asymmetry parameter.
     - Between ``-1`` and ``1``.
     - ``0.6``
   * - ``moonlight_amplification_factor``
     - ``float``
     - Empirical scaling controlling modelled moonlight influence.
     - Positive value ≥ ``5`` recommended.
     - ``5``

.. note::
   Moonlight modelling currently supports only UBVRI and ugriz photometric
   filters, in addition to telescope-valid filters.

----

MultiprocessingConfigurations
------------------------------

Controls CPU utilization during Phase Two multiprocessing.

.. code-block:: python

    from preface.configs import MultiprocessingConfigurations

    MultiprocessingConfigs = MultiprocessingConfigurations(
        toggle_multiprocessing=True,
        cores_to_leave_out=2
    )

.. list-table::
   :header-rows: 1
   :widths: 30 10 40 10 10

   * - Parameter
     - Type
     - Description
     - Valid Values
     - Default
   * - ``toggle_multiprocessing``
     - ``bool``
     - Enable multiprocessing.
     - ``True`` or ``False``.
     - ``True``
   * - ``cores_to_leave_out``
     - ``int``
     - Number of logical CPU cores reserved from computation.
     - Integer ≥ ``0``.
     - ``1``

.. warning::
   By default, one logical CPU core is reserved to maintain system
   responsiveness. Setting ``cores_to_leave_out=0`` uses all available cores
   and may cause the system to become unresponsive. PREFACE will prompt
   for confirmation before continuing in this case.
