BINNED_STELLAR_SED
==================

Name
----
BINNED_STELLAR_SED

Purpose
-------
Generates the attenuated stellar model SEDs for the non-parametric SFH.
The output model SEDs depend on the filters, age bins, and attenuation.

Calling Sequence
----------------
::

    mean_Lnu_stellar = binned_stellar_sed(stellar_models, psi [, atten_curve = , $
                          tauV_diff = , delta = , tauV_BC = , tauB_f = , F_clump = , $
                          cosi = , b_to_d = , rold0_ages = ,  atten_models = , $
                          /uv_bump, /error_check, steps_mean_Lnu=steps_mean_Lnu, $
                          mean_Lnu_unred_stellar=mean_Lnu_unred_stellar, $
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
``mean_Lnu_stellar`` : double array(Nfilters, Nmodels)
    The mean luminosity of each filter and set of stellar model parameters 
    :math:`[L_\odot\ {\rm Hz}^{-1}]`.

Optional Outputs
----------------
``steps_mean_Lnu`` : double array(Nfilters, Nsteps, Nmodels)
    The mean luminosity of each filter for each SFH age bin and set of stellar model 
    parameters :math:`[L_\odot\ {\rm Hz}^{-1}]`.
``mean_Lnu_unred_stellar`` : double array(Nfilters, Nmodels)
    The unattenuated mean luminosity of each filter and set of stellar model
    parameters :math:`[L_\odot\ {\rm Hz}^{-1}]`.
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
- 2016/05/01: Created (Rafael T. Eufrasio)
- 2020/04/27: Replaced if statements using ``n_elements`` on keywords to use ``keyword_set`` (Keith Doore)
- 2020/04/27: Fixed issue with ``doore21_Lbol_abs_table`` if a data point was at the final point in the table (Keith Doore)
- 2020/05/06: Added ability to use pure Calzetti curve (Keith Doore)
- 2020/05/06: Added ``_REF_EXTRA`` for keyword inheritance to cut down on list of keywords (Keith Doore)
- 2021/04/13: Modified Tuffs attenuation to have ``rdisk`` as intrinsic property and ``rbulge`` as B/D (Keith Doore)
- 2021/04/13: Added ``rold0_ages``, is ``1`` or ``0`` for the corresponding age bins (Keith Doore)
- 2022/03/16: Added proper error handling (Keith Doore)
- 2022/03/16: Renamed variables to standard format (Keith Doore)
- 2022/04/13: Made ``sfh_coeff`` an optional input (Keith Doore)
- 2022/04/13: Made output SED the sum of each ``sfh_coeff`` normalized step component, which is now an optional output (Keith Doore))
- 2022/04/13: Separated ``Lbol_abs`` from output SED and made it an optional output (Keith Doore)
- 2022/04/13: Made each parameter a unique optional input (Keith Doore)
- 2022/04/13: Made ``steps`` a required input (Keith Doore)
- 2022/04/13: Made the Doore+2021 ``steps_Lbol_abs`` lookup table a separate function (Keith Doore)
- 2022/06/09: Renamed some variables to match naming scheme (Keith Doore)
- 2022/06/09: Removed some keywords to replace with ``config`` (Keith Doore)
- 2022/06/09: Added ``error_check`` keyword to do error handling and added error handling (Keith Doore)
- 2022/06/30: Removed ``config`` and replaced with the three accessed configuration values (Keith Doore)
- 2022/06/30: Updated documentation (Keith Doore)
- 2022/06/30: Replaced ``!cv`` with ``!lightning_cgs`` (Keith Doore)
- 2022/07/07: Change name of ``sfh_coeff`` to ``psi`` (Keith Doore)
- 2022/07/07: Change name of ``Lbol_abs`` to ``Lbol_abs_stellar`` (Keith Doore)
- 2022/07/22: Added optional output of unreddened mean Lnu (Keith Doore)
- 2022/11/15: Fixed bug with variable redshift and ``rold0_ages`` (Keith Doore)

