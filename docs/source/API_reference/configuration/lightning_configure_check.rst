LIGHTNING_CONFIGURE_CHECK
=========================

Name
----
LIGHTNING_CONFIGURE_CHECK

Purpose
-------
Checks the Lightning configuration structure to make sure
all inputs are correctly formatted and can be understood by
Lightning. All unused parameters are removed from the
structure.

Calling Sequence
----------------
::

    config_checked = lightning_configure_check(config)

Input
-----
``config`` : structure
    A Lightning configuration structure. (See
    ``lightning_configure_defaults.pro`` for details and contents.)

Output
------
``config_checked`` : structure
    The error checked Lightning configuration structure with
    unused tags removed.

Modification History
--------------------
- 2022/04/26: Created (Keith Doore)
- 2022/05/03: Updated method for removing key from hash (Keith Doore)
- 2022/05/17: Added ability to check chosen cosmology (Keith Doore)
- 2022/05/20: Removed loguniform and lognormal prior options (Keith Doore)
- 2022/06/16: Added check for number of parallel chains if using affine MCMC (Keith Doore)
- 2022/06/29: Added conversion of ``nebular_extinction`` to ``no_nebular_extinction`` (Keith Doore)
- 2022/06/29: Added conversion of ``emission_lines`` to ``no_emission_lines`` (Keith Doore)
- 2022/07/08: Added check to make sure ``energy_balance`` is ``0`` if only using a dust model (Keith Doore)
- 2022/08/16: Added MPFIT keywords (Keith Doore)
- 2022/08/18: Added checks to ensure MCMC post-processing values are reasonable for number of trials (Keith Doore)
- 2022/09/01: Replaced ``XRAY_STAT`` with ``XRAY_UNC`` (Erik B. Monson)
- 2022/09/14: Updates to allow fitting with X-ray fluxes (Erik B. Monson)
- 2022/09/14: Fixed issue with MCMC post-processing size check for affine MCMC (Keith Doore)
- 2022/10/24: Added option to choose stranded walker deviation value for affine MCMC (Keith Doore)
- 2022/10/25: Renamed SPS to SSP (Keith Doore)
- 2022/12/13: Prevented ``XRAY_UNC`` from begin checked if ``XRAY_UNIT='FLUX'`` (Keith Doore)

