LIGHTNING_MODEL_LNUXRAY
=======================

Name
----
lightning_model_lnuxray

Purpose
-------
Generates the observed-frame luminosity density of a given lightning
X-ray model for a given set (or sets) of parameters.

Calling Sequence
----------------
::

    Lnu_mod_xray = lightning_model_lnuxray(parameters, parameter_name, models [, Lbol_AGN_model = ,$
                                           L2500 = , xray_agn_model = , agn_model = , /error_check, $
                                           Lnu_xray_stellar=Lnu_xray_stellar, $
                                           Lnu_xray_unabs_stellar=Lnu_xray_unabs_stellar, $
                                           Lnu_xray_AGN=Lnu_xray_AGN, $
                                           Lnu_xray_unabs_AGN=Lnu_xray_unabs_AGN])

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
``Lbol_AGN_model`` : double array(Nmodels)
    Bolometric luminosity of the AGN model :math:`[L_\odot]`. (Required if using a power
    law X-ray AGN emission model.)
``L2500`` : double array(Nmodels)
    The rest-frame 2500 Angstrom monochromatic luminosity shifted to the observed
    frame :math:`[L_\odot\ {\rm Hz}^{-1}]`. (Required if using a power law X-ray AGN
    emission model.)
``agn_model`` : string scalar
    The UV-to-IR AGN emission model to use. Current options are: ``'SKIRTOR'`` or ``'NONE'``.
    (Default = ``'NONE'``)
``xray_agn_model`` : string scalar
    The X-ray AGN emission model to use. Current options are: ``'PLAW'``, ``'QSOSED'``, ``'NONE'``.
    (Default = ``'QSOSED'``)
``error_check`` : flag
    If set, all inputs are checked for errors. Otherwise, all inputs are
    assumed to be of correct format.

Output
------
``Lnu_mod_xray`` : double array(Nxray, Nmodels)
    The total X-ray SED, after absorption (intrinsic and Galactic), convolved with
    the bandpasses defined in ``models.XRAY_MODELS`` :math:`[L_\odot\ {\rm Hz}^{-1}]`.

Optional Outputs
----------------
``Lnu_xray_stellar`` : double array(Nwave, Nmodels)
    Stellar emission component of ``Lnu_mod_xray`` :math:`[L_\odot\ {\rm Hz}^{-1}]`.
``Lnu_xray_unabs_stellar`` : double array(Nwave, Nmodels)
    Intrinsic stellar X-ray SED for each set of parameters :math:`[L_\odot\ {\rm Hz}^{-1}]`.
``Lnu_xray_AGN`` : double array(Nwave, Nmodels)
    AGN emission component of ``Lnu_mod_xray`` :math:`[L_\odot\ {\rm Hz}^{-1}]`.
``Lnu_xray_unabs_AGN`` : double array(Nwave, Nmodels)
    Intrinsic AGN X-ray SED for each set of parameters :math:`[L_\odot\ {\rm Hz}^{-1}]`.

Modification History
--------------------
- 2022/09/08: Created (Erik B. Monson)
- 2022/09/15: Updated documentation (Keith Doore)
- 2022/11/02: Galactic NH is now in units of 1e20 cm-2 (Erik B. Monson)

