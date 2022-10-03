LIGHTNING_MODEL_COUNTS
======================

Name
----
LIGHTNING_MODEL_COUNTS

Purpose
-------
Generates the X-ray detector counts (folded through the instrumental 
response) of a lightning model under a specified set(s) of parameters.

Calling Sequence
----------------
::

    counts = lightning_model_counts(parameters, parameter_name, models [, Lbol_AGN_model = , $
                                    L2500 = , agn_model = , xray_agn_model = , $
                                    /error_check, counts_xrb=counts_xrb, counts_agn=counts_agn])

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
    frame :math:`[L_\odot\ {\rm Hz}^{-1}]`. (Required if using a power law X-ray
    AGN emission model.)
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
``counts`` : double array(Nxray, Nmodels)
    The total instrumental counts produced by the X-ray model under the set
    observing conditions :math:`[{\rm counts}]`.

Optional Outputs
----------------
``counts_xrb`` : double array(Nxray, Nmodels)
    The instrumental counts produced by the X-ray binary model :math:`[{\rm counts}]`.
``counts_agn`` : double array(Nxray, Nmodels)
    The instrumental counts produced by the X-ray AGN model :math:`[{\rm counts}]`.

Modification History
--------------------
- 2022/06/07: Removed AGN covering factor from model (Erik B. Monson)
- 2022/06/20: Major update to include new implementation (e.g., prior, config, etc.) (Keith Doore)
- 2022/06/20: Replaced ``!cv`` with ``!lightning_cgs`` (Keith Doore)
- 2022/06/20: Updated documentation (Keith Doore)
- 2022/06/20: Added error handling (Keith Doore)
- 2022/06/20: Added ``error_check`` keyword to do error handling (Keith Doore)
- 2022/07/06: Removed ``config`` and replaced ``config`` tag calls with inputs (Keith Doore)
- 2022/07/06: Transposed parameters array to eliminate need to reform it after indexing. (Keith Doore)

