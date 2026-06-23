Outputs
=======

PREFACE writes all generated files to the directory specified by
``OutputConfigurations.output_folder``. Outputs are organized across two
processing phases:

.. code-block:: text

    output_folder/
    ├── phase_1/
    │   ├── nonviable_target_list/
    │   └── viable_target_list/
    └── phase_2/
        ├── cumulative_observability_scores/
        ├── full_ranked_event_list/
        ├── graphs/
        └── individual_planets/

----

Phase One
---------

Phase One selects the most promising planetary candidates for transmission
spectroscopy from catalogues, ranking them proportionally to expected data
signal-to-noise ratio.

.. list-table::
   :header-rows: 1
   :widths: 35 65

   * - Directory
     - Contents
   * - ``nonviable_target_list/``
     - Targets with metric values below the ``viable_cumulative_cut`` threshold.
   * - ``viable_target_list/``
     - Targets retained after applying the cumulative ranking cutoff, carried
       forward to Phase Two.

----

Phase Two
---------

Phase Two determines the best exoplanet observation dates from predicted transits
at the selected telescope location.

.. list-table::
   :header-rows: 1
   :widths: 40 60

   * - Directory
     - Contents
   * - ``cumulative_observability_scores/``
     - Observability score accumulated from all transits per planet.
   * - ``full_ranked_event_list/``
     - Complete ranked list of individual observable transit events.
   * - ``graphs/``
     - Diagnostic plots for each transit (when ``toggle_graph_outputs=True``).
   * - ``individual_planets/``
     - Per-planet computed transits and their individual scores.

