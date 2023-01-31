LIGHTNING_POSTPROCESSING
========================

Name
----
LIGHTNING_POSTPROCESSING

Purpose
-------
Post-processes the fitting results into more readily usable outputs.
The final post-processed output is saved to a table in the first
extension of a FITS file.

Calling Sequence
----------------
::

    lightning_postprocessing, input_dir, config, sed_id

Inputs
------
``input_dir`` : string scalar
    The path to the file containing the input SED data.
``config`` : structure
    A Lightning configuration structure. (See
    ``lightning_configure_defaults.pro`` for details and contents.)
``sed_id`` : string array(Nsed)
    The ID of each SED.

Output
------
A data table saved to the first extension of a FITS file with a file name (and path) of
``<input_dir>/lightning_output/postprocessed_data_<date_in_UTC>.fits.gz``. The table contains the
post-processed data for all SEDs.
See :ref:`postprocessing-label` for full details and possible contents.

Modification History
--------------------
- 2022/07/12: Created (Keith Doore)
- 2022/07/29: Added check to remove occasional stranded walker in affine-MCMC from final distribution (Keith Doore)
- 2022/08/01: Fixed chain thinning bug (Keith Doore)
- 2022/08/01: Changed to element-wise flattening from segment-wise flattening for affine MCMC (Keith Doore)
- 2022/08/07: Fixed high resolution model indexing bug (Keith Doore)
- 2022/08/09: Fixed incorrect initialization size of high resolution xray models (Keith Doore)
- 2022/08/09: Fixed issue where ``config.HIGH_RES_MODEL_FRACTION=0`` was getting no models vs best fit (Keith Doore)
- 2022/08/11: Updated high resolution xray model function name (Keith Doore)
- 2022/08/17: Updated to include MPFIT outputs (Keith Doore)
- 2022/08/17: Included ``PARAMETER_NAMES`` tag to determine certain corresponding outputs (Keith Doore)
- 2022/08/18: Added progress printing (Keith Doore)
- 2022/08/18: Fixed burn-in and thinning issues with MCMC chain if ``Ntrials`` is small (Keith Doore)
- 2022/08/23: Added date to the end of output file name (Keith Doore)
- 2022/08/24: Updated to allow for different xray bandpasses, which were not checked at input (Keith Doore)
- 2022/08/24: Updated to make fixed parameters only be scalar vs array of same value (Keith Doore)
- 2022/09/01: Added handling for user-supplied X-ray count uncertainties (Erik B. Monson)
- 2022/09/01: Added ``chi2`` output in addition to ``lnprob`` (Erik B. Monson)
- 2022/09/01: Added X-ray emission to PPC (Erik B. Monson)
- 2022/09/14: Updates to allow fitting with X-ray fluxes (Erik B. Monson)
- 2022/09/21: Few bug fixes with X-ray fluxes (Keith Doore)
- 2022/09/22: Updated how we identify stranded walkers in affine MCMC (Keith Doore)
- 2022/09/23: Made ``autocorr_flag`` unique to affine MCMC (Keith Doore)
- 2022/09/26: Added ``DOF`` to MPFIT output (Keith Doore)
- 2022/10/24: Updated stranded walker search to use configuration input value (Keith Doore)
- 2022/10/25: Renamed SPS to SSP (Keith Doore)
- 2023/01/16: Fixed issue if ``lnprob`` of MPFIT is same for multiple solvers (Keith Doore)
- 2023/01/23: Included all solvers parameter values is using MPFIT for user convergence testing (Keith Doore)
- 2023/01/31: Updated to use ``config.OUTPUT_FILENAME`` to set the filename (Keith Doore)
- 2023/01/31: Replaced ``STEPS_MSTAR_COEFF`` output with ``MSTAR`` and ``STEPS_MSTAR`` (Keith Doore)

