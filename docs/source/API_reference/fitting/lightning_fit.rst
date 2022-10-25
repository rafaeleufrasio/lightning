LIGHTNING_FIT
=============

Name
----
LIGHTNING_FIT

Purpose
-------
Generates the models and priors specified in the Lightning configuration
structure, and passes them and the SED data to the specified fitting
algorithm.

Calling Sequence
----------------
::

    lightning_fit, input_dir, sed_data, config

Inputs
------
``input_dir`` : string scalar
    The path to the file containing the input SED data.
``sed_data`` : structure
    A structure containing the SED luminosities and uncertainties, filter
    labels, distances, redshifts, and optional X-ray data. (See
    ``lightning_input.pro`` for details and contents.)
``config`` : structure
    A Lightning configuration structure. (See
    ``lightning_configure_defaults.pro`` for details and contents.)

Modification History
--------------------
- 2022/05/13: Created (Keith Doore)
- 2022/07/11: Moved models generation to separate function (Keith Doore)
- 2022/07/27: Added non-parametric SFH bin cutting (Keith Doore)
- 2022/08/09: Added ``galactic_nh`` as input for xray emission from sed_data (Keith Doore)
- 2022/08/17: Fixed issue where ``config.PSI`` was being updated and returned to
  ``lightning.pro`` if ``config.MAXCPU = 1`` (Keith Doore)
- 2022/09/15: Removed different calls to ``lightning_models.pro`` for X-ray emission and
  no X-ray emission, since ``lightning_input.pro`` now handles including X-ray emission if
  requested in configuration. (Keith Doore)
- 2022/10/25: Renamed SPS to SSP (Keith Doore)

