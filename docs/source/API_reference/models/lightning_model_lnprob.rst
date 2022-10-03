LIGHTNING_MODEL_LNPROB
======================

Name
----
LIGHTNING_MODEL_LNPROB

Purpose
-------
Calculates the natural log of the model(s) probability
(i.e. ln(prior * likelihood)) given the data for a set
of parameters.

Calling Sequence
----------------
::

    lnprob = lightning_model_lnprob(sed_data, parameters, config_nopriors, models, priors [, $
                                    /negative, /error_check, UVIR_chi2=UVIR_chi2, $
                                    xray_chi2=xray_chi2, lnlike=lnlike, lnprior=lnprior, $
                                    Lnu_mod=Lnu_mod, xray_counts_mod=xray_counts_mod, $
                                    Lnu_xray_mod=Lnu_xray_mod])

Inputs
------
``sed_data`` : structure
    A structure containing the SED luminosities and uncertainties, filter
    labels, distances, redshifts, and optional X-ray data. (See
    ``lightning_input.pro`` for details and contents.)
``parameters`` : int, float, or double array(Nparam, Nmodels)
    Parameters of the model(s). The actual parameters contained in this
    array depend on the chosen model(s) during configuration.
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

Optional Inputs
---------------
``negative`` : flag
    If set, the function returns the negative of the log probability
    (for minimization purposes).
``error_check`` : flag
    If set, all inputs are checked for errors. Otherwise, all inputs are
    assumed to be of correct format.

Output
------
``lnprob`` : double array(Nmodels)
    The log probability of the model(s) for the specified parameters.

Optional Outputs
----------------
``UVIR_chi2`` : double array(Nmodels)
    The :math:`\chi^2` of the UV-to-IR portion of the model(s).
``xray_chi2`` : double array(Nmodels)
    The X-ray statistic of the model(s).
``lnlike`` : double array(Nmodels)
    The total likelihood log probability of the model(s).
``lnprior`` : double array(Nmodels)
    The total prior log probability of each model(s).
``Lnu_mod`` : double array(Nfilters, Nmodels)
    The model Lnu for each filter of each model generated from the
    specified parameters :math:`[L_\odot\ {\rm Hz}^{-1}]`.
``xray_counts_mod`` : double array(Nxray, Nmodels)
    The model X-ray counts for each band of each model generated from
    the specified parameters :math:`[{\rm counts}]`.
``Lnu_xray_mod`` : double array(Nxray, Nmodels)
    The X-ray model Lnu for each band of each model generated from the
    specified parameters :math:`[L_\odot\ {\rm Hz}^{-1}]`.

Modification History
--------------------
- 2022/02/01: Created (Erik Monson)
- 2022/06/07: Major update to include new implementation (e.g., prior, config, etc.) (Keith Doore)
- 2022/07/05: Renamed config to ``config_nopriors`` to reflect removed parameters (Keith Doore)
- 2022/07/25: Fixed bug where in bounds indexing for ``AGN_mass`` and ``AGN_logmdot`` in qsosed ``L2500`` computation was missing (Keith Doore)
- 2022/09/01: Added handling for user-supplied X-ray count uncertainties (Erik B. Monson)
- 2022/09/14: Updates to allow fitting with X-ray fluxes (Erik B. Monson)

