LIGHTNING_MODEL_LNU_HIGHRES
===========================

Name
----
LIGHTNING_MODEL_LNU_HIGHRES

Purpose
-------
Generates the observed-frame luminosity spectrum of a given lightning model
for a given set (or sets) of parameters.

Calling Sequence
----------------
::

    Lnu_mod_highres = lightning_model_lnu_highres(parameters, parameter_name, models [, sps = , $
                                                  sfh = , dust_model = , agn_model = , /energy_balance, $
                                                  /error_check, /get_wave, Lnu_stellar_highres=Lnu_stellar_highres, $
                                                  wave=wave, Lnu_unred_stellar_highres=Lnu_unred_stellar_highres, $
                                                  Lnu_AGN_highres=Lnu_AGN_highres, Lnu_unred_AGN_highres=Lnu_unred_AGN_highres, $
                                                  Lnu_dust_highres=Lnu_dust_highres, LTIR=LTIR, $
                                                  Lbol_AGN_model=Lbol_AGN_model, L2500=L2500, _extra = ])

Inputs
------
``parameters`` : int, float, or double array(Nparam, Nmodels)
    Parameters of the model(s). The actual parameters contained in this
    array depend on the chosen model(s) during configuration.
``parameter_name`` : string array(Nparam)
    The names associated with the parameters of the model(s) given in the
    same order as ``parameters``.
``models`` : structure
    A structure containing each model structure (stellar, dust, AGN,
    X-ray) as a substructure. (See ``lightning_models.pro`` for details
    and contents.)

Optional Inputs
---------------
``sps`` : string scalar
    The stellar population synthesis (SPS) models to use for the stellar
    population. Current options  are: ``'PEGASE'`` or ``'NONE'``.
    (Default = ``'PEGASE'``)
``sfh`` : string scalar
    The type of SFH to assume if the model is to include a stellar
    population. Current options are: ``'NON-PARAMETRIC'``.
    (Default = ``'NON-PARAMETRIC'``)
``dust_model`` : string scalar
    The dust emission model to use. Current options are: ``'DL07'`` or ``'NONE'``.
    (Default = ``'DL07'``)
``agn_model`` : string scalar
    The UV-to-IR AGN emission model to use. Current options are: ``'SKIRTOR'`` or ``'NONE'``.
    (Default = ``'NONE'``)
``energy_balance`` : flag
    If set, the total integrated IR luminosity (normalization) of the dust emission
    is tied to the total absorbed stellar (and, if set, AGN) emission.
``error_check`` : flag
    If set, all inputs are checked for errors. Otherwise, all inputs are
    assumed to be of correct format.
``wave`` : int, float, or double array(Nwave)
    The common wavelength grid of the spectra :math:`[\mu \rm m]`.
``get_wave`` : flag
    If set, the function returns a value of ``-1`` after generating the default wavelength grid.
    This is useful when one wants to get the default wavelength grid without computing
    any high resolution spectra.
``L2500`` : double array(Nmodels)
    The rest-frame 2500 Angstrom monochromatic luminosity shifted to the observed
    frame :math:`[L_\odot\ {\rm Hz}^{-1}]`. (Used as input if using a QSOSED X-ray
    AGN model.)
``_extra`` : structure
    Additional optional inputs that are passed to ``binned_stellar_spectrum.pro`` and
    ``skirtor_spectrum.pro``.

Output
------
``Lnu_mod_highres`` : double array(Nwave, Nmodels)
    The model luminosity spectrum of each model generated on a common wavelength
    grid from the specified parameters :math:`[L_\odot\ {\rm Hz}^{-1}]`.

Optional Outputs
----------------
``wave`` : double array(Nwave)
    The common wavelength grid of the spectra :math:`[\mu \rm m]`.
``Lnu_stellar_highres`` : double array(Nwave, Nmodels)
    Stellar emission component of ``Lnu_mod_highres`` :math:`[L_\odot\ {\rm Hz}^{-1}]`.
``Lnu_unred_stellar_highres`` : double array(Nfilters, Nmodels)
    The unattenuated stellar emission of each model generated on a common wavelength
    grid from the specified parameters :math:`[L_\odot\ {\rm Hz}^{-1}]`.
``Lnu_AGN_highres`` : double array(Nwave, Nmodels)
    AGN emission component of ``Lnu_mod_highres`` :math:`[L_\odot\ {\rm Hz}^{-1}]`.
``Lnu_unred_AGN_highres`` : double array(Nwave, Nmodels)
    The unattenuated AGN emission of each model generated on a common wavelength
    grid from the specified parameters :math:`[L_\odot\ {\rm Hz}^{-1}]`.
``Lnu_dust_highres`` : double array(Nwave, Nmodels)
    Dust emission component of ``Lnu_mod_highres`` :math:`[L_\odot\ {\rm Hz}^{-1}]`.
``LTIR`` : double array(Nmodels)
    The total integrated IR luminosity (i.e., the bolometric luminosity)
    of the dust model for each set of model parameters :math:`[L_\odot]`.
``Lbol_AGN_model`` : double array(Nmodels)
    Bolometric luminosity of the AGN model :math:`[L_\odot]`.
``L2500`` : double array(Nmodels)
    The rest-frame 2500 Angstrom monochromatic luminosity shifted to the observed
    frame for the current AGN model or qsosed X-ray AGN model if input
    :math:`[L_\odot\ {\rm Hz}^{-1}]`.

Modification History
--------------------
- 2022/04/14: Created (Erik B. Monson)
- 2022/07/06: Major update to match ``lightning_model_lnu`` (Keith Doore)
- 2022/07/22: Added optional output of unreddened ``Lnu_stellar_highres`` (Keith Doore)
- 2022/07/22: Added optional output of the common wavelength grid (Keith Doore)
- 2022/07/22: Added optional output of unreddened ``Lnu_AGN_highres`` (Keith Doore)
- 2022/07/27: Renamed unreddened Lnus to prevent ambiguous keywords (Keith Doore)
- 2022/09/08: Allow user-specified arbitary wavelength grid (Erik B. Monson)
- 2022/10/25: Renamed SPS to SSP (Keith Doore)
- 2022/12/13: Fixed bug in ``wave`` error check (Keith Doore)

