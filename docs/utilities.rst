Miscellaneous Utilities
==================

PREFACE ships with several convenience utilities for inspecting supported
instruments and managing intermediate files.

Listing Supported Telescopes
-----------------------------

.. code-block:: python

    import preface

    print(preface.telescope_list)

Returns a list of all instrument names recognised by PREFACE.

Retrieving Available Filters
-----------------------------

.. code-block:: python

    import preface

    preface.get_available_filters_list("TNT ULTRASPEC")

Returns the photometric filters supported by the specified instrument.

Opening the Telescope Reference Table
--------------------------------------

.. code-block:: python

    import preface

    preface.open_scope_csv()

Opens the bundled telescope reference CSV in your default application, providing
a full overview of instrument parameters used internally by PREFACE.

Clearing Intermediate Files
----------------------------

.. code-block:: python

    import preface

    preface.wipe_intermediate_csvs()

Removes temporary intermediate CSV files generated during pipeline execution,
freeing up disk space after a run is complete.

