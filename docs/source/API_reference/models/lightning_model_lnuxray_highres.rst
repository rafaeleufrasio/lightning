LIGHTNING_MODEL_LNUXRAY_HIGHRES
===============================

Name
----
LIGHTNING_MODEL_LNUXRAY_HIGHRES

Purpose
-------
Generates the observed-frame luminosity density of a given lightning
X-ray model for a given set (or sets) of parameters.

Calling Sequence
----------------
::

    Lnu_xray_highres = lightning_model_lnuxray_highres(parameters, parameter_name, models [, Lbol_AGN_model = ,$
                                                  L2500 = , xray_agn_model = , agn_model = , /error_check, $
                                                  /get_wave, wave=wave, Lnu_xray_stellar_highres=Lnu_xray_stellar_highres, $
                                                  Lnu_xray_unabs_stellar_highres=Lnu_xray_unabs_stellar_highres, $
                                                  Lnu_xray_AGN_highres=Lnu_xray_AGN_highres, $
                                                  Lnu_xray_unabs_AGN_highres=Lnu_xray_unabs_AGN_highres, $
                                                  Lnu_xray_mod=Lnu_xray_mod])

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
``wave`` : int, float, or double array(Nwave)
    The common wavelength grid of the spectra :math:`[\mu \rm m]`.
``get_wave`` : flag
    If set, the function returns a value of ``-1`` after generating the default wavelength grid.
    This is useful when one wants to get the default wavelength grid without computing
    any high resolution spectra.

Output
------
``Lnu_xray_highres`` : double array(Nwave, Nmodels)
    The total, high resolution X-ray spectral model, after absorption (intrinsic and
    Galactic), calculated on a common wavelength grid for each set of parameters
    :math:`[L_\odot\ {\rm Hz}^{-1}]`.

Optional Outputs
----------------
``wave`` : double array(Nwave)
    The common wavelength grid of the spectra :math:`[\mu \rm m]`.
``Lnu_xray_stellar_highres`` : double array(Nwave, Nmodels)
    Stellar emission component of ``Lnu_xray_highres`` :math:`[L_\odot\ {\rm Hz}^{-1}]`.
``Lnu_xray_unabs_stellar_highres`` : double array(Nwave, Nmodels)
    Intrinsic high resolution stellar X-ray spectral model calculated on a
    common wavelength grid for each set of parameters :math:`[L_\odot\ {\rm Hz}^{-1}]`.
``Lnu_xray_AGN_highres`` : double array(Nwave, Nmodels)
    AGN emission component of ``Lnu_xray_highres`` :math:`[L_\odot\ {\rm Hz}^{-1}]`.
``Lnu_xray_unabs_AGN_highres`` : double array(Nwave, Nmodels)
    Intrinsic high resolution AGN X-ray spectral model calculated on a common
    wavelength grid for each set of parameters :math:`[L_\odot\ {\rm Hz}^{-1}]`.
``Lnu_xray_mod`` : double array(Nbands)
    The mean luminosity calculated from ``Lnu_xray_highres`` using each of the X-ray
    bandpasses :math:`[L_\odot\ {\rm Hz}^{-1}]`.

Modification History
--------------------
- 2022/07/25: Created (Erik B. Monson)
- 2022/07/25: Renamed variables to match common naming scheme (Keith Doore)
- 2022/07/25: Added error handling with ``error_check`` keyword (Keith Doore)
- 2022/07/27: Added fixed observed-frame wavelength grid (Keith Doore)
- 2022/08/11: Renamed function to match function naming scheme (Keith Doore)
- 2022/09/08: Allow user-specified arbitary wavelength grid (Erik B. Monson)
- 2022/11/02: Galactic NH is now in units of 1e20 cm-2 (Erik B. Monson)
- 2022/12/13: Fixed bug in ``wave`` error check (Keith Doore)

