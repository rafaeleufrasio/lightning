LIGHTNING_MODEL_LNU
===================

Name
----
LIGHTNING_MODEL_LNU

Purpose
-------
Generates the observed-frame SED luminosities of a given lightning model
for a given set (or sets) of parameters.

Calling Sequence
----------------
::

    Lnu_mod = lightning_model_lnu(parameters, parameter_name, models [, sps = , sfh = , $
                                  dust_model = , agn_model = , /energy_balance, /error_check, $
                                  Lnu_stellar=Lnu_stellar, Lnu_unred_stellar=Lnu_unred_stellar,$
                                  Lnu_unred_AGN=Lnu_unred_AGN, Lnu_AGN=Lnu_AGN, Lnu_dust=Lnu_dust, $
                                  LTIR=LTIR, Lbol_AGN_model=Lbol_AGN_model, L2500=L2500, $
                                  _extra = ])

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
``L2500`` : double array(Nmodels)
    The rest-frame 2500 Angstrom monochromatic luminosity shifted to the observed 
    frame :math:`[L_\odot\ {\rm Hz}^{-1}]`. (Used as input if using a QSOSED X-ray AGN model.)
``_extra`` : structure
    Additional optional inputs that are passed to ``binned_stellar_sed.pro`` and
    ``skirtor_sed.pro``.

Output
------
``Lnu_mod`` : double array(Nfilters, Nmodels)
    The model Lnu for each band of each model generated from the
    specified parameters :math:`[L_\odot\ {\rm Hz}^{-1}]`.

Optional Outputs
----------------
``Lnu_stellar`` : double array(Nfilters, Nmodels)
    Stellar emission component of ``Lnu_mod`` :math:`[L_\odot\ {\rm Hz}^{-1}]`.
``Lnu_unred_stellar`` : double array(Nfilters, Nmodels)
    The unattenuated stellar emission for each band of each model
    :math:`[L_\odot\ {\rm Hz}^{-1}]`.
``Lnu_AGN`` : double array(Nfilters, Nmodels)
    AGN emission component of ``Lnu_mod`` :math:`[L_\odot\ {\rm Hz}^{-1}]`.
``Lnu_unred_AGN`` : double array(Nfilters, Nmodels)
    The unattenuated AGN emission for each band of each model :math:`[L_\odot\ {\rm Hz}^{-1}]`.
``Lnu_dust`` : double array(Nfilters, Nmodels)
    Dust emission component of ``Lnu_mod`` :math:`[L_\odot\ {\rm Hz}^{-1}]`.
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
- 2022/02/01: Created (Erik Monson)
- 2022/06/08: Major update to include new implementation (e.g., prior, config, etc.) (Keith Doore)
- 2022/06/27: Updated logic for ``LTIR`` to include a dust emission model (Keith Doore)
- 2022/06/30: Rearranged component Lnu generation to allow for more straightforward computations (Keith Doore)
- 2022/07/05: Removed ``config`` and replaced ``config`` tag calls with inputs (Keith Doore)
- 2022/07/05: Added ``_extra`` to hold extra optional inputs for stellar and agn models not directly used in this function (Keith Doore)
- 2022/07/05: Transposed ``parameters`` array to eliminate need to reform it after indexing. (Keith Doore)
- 2022/07/22: Added optional output of unreddened ``Lnu_stellar`` (Keith Doore)
- 2022/07/22: Added optional output of unreddened ``Lnu_AGN`` (Keith Doore)
- 2022/07/27: Renamed unreddened Lnus to prevent ambiguous keywords (Keith Doore)
- 2022/10/25: Renamed SPS to SSP (Keith Doore)

