BINNED_STELLAR_SPECTRUM
=======================

Name
----
BINNED_STELLAR_SPECTRUM

Purpose
-------
Generates the attenuated stellar model spectra for the non-parametric SFH.
The output model spectra depend on the age bins and attenuation.

Calling Sequence
----------------
::

    Lnu_spec_stellar = binned_stellar_spectrum(stellar_models, psi [, atten_curve = , $
                                    tauV_diff = , delta = , tauV_BC = , tauB_f = , F_clump = , $
                                    cosi = , b_to_d = , rold0_ages = , atten_models = , /uv_bump, $
                                    /error_check, steps_Lnu_spec=steps_Lnu_spec, $
                                    Lnu_spec_unred_stellar=Lnu_spec_unred_stellar,$
                                    Lbol_abs_stellar=Lbol_abs_stellar, steps_Lbol_abs=steps_Lbol_abs])

Inputs
------
``stellar_models`` : structure
    A structure containing the spectra, SEDs, and stellar parameters for the 
    non-parametric stellar model. (See ``binned_stellar_models.pro`` for
    details and contents.)
``psi`` : int, float, or double array(Nsteps, Nmodels)
    The non-parametric SFH coefficients :math:`[M_\odot\ {\rm yr}^{-1}]`.

Optional Inputs
---------------
``atten_curve`` : string scalar
    The name of the attenuation curve to apply to the stellar models. Current
    options are ``'CALZETTI00'``, ``'CALZETTI_MOD'``, and ``'DOORE21'``.
    (Default = ``'CALZETTI00'``)
``tauV_diff`` : int, float, or double array(Nmodels)
    The V-band optical depth of the diffuse dust for the Calzetti attenuation.
    (Default = ``1.0``)
``delta`` : int, float, or double array(Nmodels)
    The power law value to change the attenuation curve slope for the Calzetti
    attenuation. (Default = ``0.d0``)
``tauV_BC`` : int, float, or double array(Nmodels)
    The V-band optical depth of the birth cloud for the Calzetti attenuation.
    (Default = ``0.0``)
``tauB_f`` : int, float, or double array(Nmodels)
    The face-on optical depth in the B-band for the Doore attenuation.
    (Default = ``1.0``)
``F_clump`` : int, float, or double array(Nmodels)
    The clumpiness factor F for the Doore attenuation. (Default = ``0.0``)
``cosi`` : int, float, or double array(Nmodels)
    The inclination of the galactic disk in terms of cos(i) for the Doore
    attenuation. (Default = ``1.d0``)
``b_to_d`` : int, float, or double array(Nmodels)
    The bulge-to-disk ratio for the Doore attenuation. (Default = ``0.0``)
``rold0_ages`` : int, float, or double array(Nsteps)
    The binary parameter ``rold0``, designating each SFH age bin as part of
    the young or old population when using the Doore attenuation. A value
    of ``0`` for the corresponding age bin considers it to be part of the young
    population, and a value of ``1`` considers it to be part of the old
    populations. (Default = Ages < 500 Myr are ``0`` else ``1``)
``atten_models`` : structure
    A structure containing the preloaded files for the Doore attenuation.
``uv_bump`` : flag
    If set, then a 2175 Angstrom UV bump feature will be added to the 
    attenuation curve.
``error_check`` : flag
    If set, all inputs are checked for errors. Otherwise, all inputs are
    assumed to be of correct format.

Output
------
``Lnu_spec_stellar`` : double array(Nwave, Nmodels)
    The luminosity spectrum for each set of stellar model parameters at
    each wavelength :math:`[L_\odot\ {\rm Hz}^{-1}]`.

Optional Outputs
----------------
``steps_Lnu_spec`` : double array(Nwave, Nsteps, Nmodels)
    The luminosity spectrum for each SFH age bin and set of stellar model 
    parameters at each wavelength :math:`[L_\odot\ {\rm Hz}^{-1}]`.
``Lnu_spec_unred_stellar`` : double array(Nwave, Nmodels)
    The unattenuated luminosity spectrum for each set of stellar model
    parameters at each wavelength :math:`[L_\odot\ {\rm Hz}^{-1}]`.
``Lbol_abs_stellar`` : double array(Nmodels)
    The bolometric luminosity of the attenuated stellar emission for each
    set of stellar model parameters :math:`[L_\odot]`.
``steps_Lbol_abs`` : double array(Nsteps, Nmodels)
    The bolometric luminosity of the attenuated stellar emission for
    each SFH age bin and set of stellar model parameters :math:`[L_\odot]`.

Note
----
Defaults will only be set if the optional ``error_check`` input is set.

Modification History
--------------------
- 2022/04/13: Created (Keith Doore)
- 2022/06/30: Updated to match format of ``binned_stellar_sed.pro`` (Keith Doore)
- 2022/07/07: Change name of ``sfh_coeff`` to ``psi`` (Keith Doore)
- 2022/07/22: Added optional output of unreddened Lnu spectrum (Keith Doore)
- 2022/11/15: Fixed bug with variable redshift and ``rold0_ages`` (Keith Doore)

