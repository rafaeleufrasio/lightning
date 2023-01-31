LIGHTNING_CONFIGURE_INTERACTIVE
===============================

Name
----
LIGHTNING_CONFIGURE_INTERACTIVE

Purpose
-------
Generates the Lightning configuration structure via interactive
terminal prompts. Inputs are checked for errors to make sure all
inputs are correctly formatted and can be understood by Lightning.
Additionally, already existing configuration structures can be edited
interactively if it is optionally input. Further details can be
found at :ref:`configure-setting-label`.

Calling Sequence
----------------
::

    config = lightning_configure_interactive([config_edit = ])

Optional Input
--------------
``config_edit`` : structure
    A Lightning configuration structure that is to be
    interactively edited.

Output
------
``config`` : structure
    A Lightning configuration structure. (See
    ``lightning_configure_defaults.pro`` for details and contents.)

Modification History
--------------------
- 2022/04/28: Created (Keith Doore)
- 2022/05/05: Allowed for interactively editing structure (Keith Doore)
- 2022/05/16: Added ability to choose no stellar model (Keith Doore)
- 2022/05/17: Added ability to choose cosmology (Keith Doore)
- 2022/05/20: Removed loguniform and lognormal prior options (Keith Doore)
- 2022/08/30: Placed cosmology parameters at end of Core section for online documentation purposes (Keith Doore)
- 2022/09/01: Replaced ``XRAY_STAT`` with ``XRAY_UNC`` (Erik B. Monson)
- 2022/09/14: Updates to allow fitting with X-ray fluxes (Erik B. Monson)
- 2022/10/24: Added option to choose stranded walker deviation value for affine MCMC (Keith Doore)
- 2022/10/25: Renamed SPS to SSP (Keith Doore)
- 2022/12/13: Prevented ``XRAY_UNC`` from begin set if ``XRAY_UNIT='FLUX'`` (Keith Doore)
- 2023/01/31: Added ``OUTPUT_FILENAME`` option to allow for setting of post-processed filename (Keith Doore)

