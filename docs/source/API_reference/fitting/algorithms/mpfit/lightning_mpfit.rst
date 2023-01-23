LIGHTNING_MPFIT
===============

Name
----
LIGHTNING_MPFIT

Purpose
-------
Fits an observed SED with the Lightning models using the IDL MPFIT package,
which utilizes the Levenberg-Marquardt gradient decent algorithm. Initial
parameter positions are randomly selected from the specified prior distribution.

Calling Sequence
----------------
::

    lightning_mpfit, input_dir, sed_data, config_nopriors, models, priors

Inputs
------
``input_dir`` : string scalar
    The path to the file containing the input SED data.
``sed_data`` : structure
    A structure containing the SED luminosities and uncertainties, filter
    labels, distances, redshifts, and optional X-ray data. (See 
    ``lightning_input.pro`` for details and contents.)
``config_nopriors`` : structure
    A Lightning configuration structure edited to remove the prior 
    substructures for each parameter. (See ``lightning_configure_defaults.pro``
    for details and contents.)
``models`` : structure
    A structure containing each model structure (stellar, dust, AGN, 
    X-ray) as a substructure. (See ``lightning_models.pro`` for details
    and contents.)
``priors`` : structure
     A structure containing the prior hyper-parameters. (See
     ``generate_prior_struct.pro`` for details and contents.)

Output
------
An IDL save file saved to ``<input_dir>/lightning_output/output_sav_files/`` named 
``lightning_output_<galaxy_id>.sav'``, containing the resulting MPFIT best fit
parameters and log probability for each set of parameter starting points, and the
convergence metrics of the fits.

References
----------
- `Markwardt, C. B. 2009, ASP Conf. Ser. 411, 251 <https://ui.adsabs.harvard.edu/abs/2009ASPC..411..251M/abstract>`_
- `Moré, J. 1978, Numerical Analysis, 630, 105 <https://ui.adsabs.harvard.edu/abs/1978LNM...630..105M/abstract>`_

Modification History
--------------------
- 2022/01/01: Created (Rafael Eufrasio)
- 2022/08/15: Major update to include new implementation (e.g., prior, config, etc.) (Keith Doore)
- 2022/08/18: Added progress printing (Keith Doore)
- 2023/01/23: Fixed bug where if ``status <= 0`` then the corresponding ``lnprob_mpfit = 0``.
  Now have ``lnprob_mpfit = NaN`` if ``status <= 0`` (Keith Doore)

