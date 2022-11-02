SKIRTOR_SED
===========

Name
----
SKIRTOR_SED

Purpose
-------
Generates the AGN model SED for a given set of AGN model parameters and
filters. The AGN model is based on the SKIRTOR models (Stalevski et al. 2016),
which include the parameters of :math:`\tau_{9.7}` and :math:`i_{\rm AGN}`.

Calling Sequence
----------------
::

    mean_Lnu_AGN = skirtor_sed(skirtor [, tau97 = , i_agn = , tauV_diff = , delta = , $
                               atten_curve = , /uv_bump, /error_check, /xray_isotropic, $
                               mean_Lnu_unred_AGN=mean_Lnu_unred_AGN, Lbol_abs_AGN=Lbol_abs_AGN, $
                               Lbol_AGN=Lbol_AGN, L2500=L2500])

Input
-----
``skirtor`` : structure
    A structure containing the spectra, SEDs, and AGN parameters for the
    SKIRTOR model. (See ``skirtor_models.pro`` for details and contents.)

Optional Inputs
---------------
``tau97`` : int, float, or double array(Nmodels)
    Edge-on optical depth of the AGN torus at 9.7 um. (Default = ``3.0``)
``i_agn`` : int, float, or double array(Nmodels)
    Inclination from the polar axis to the line of sight :math:`[{\rm degrees}]`.
    (Default = ``0.0``)
``tauV_diff`` : int, float, or double array(Nmodel)
    The V-band optical depth of the diffuse dust. If this input is set,
    the model will be attenuated.
``delta`` : int, float, or double array(Nmodels)
    The power law value to change the attenuation curve slope.
``atten_curve`` : string scalar
    The name of the attenuation curve to apply to the AGN models. Current
    options are ``'CALZETTI00'`` and ``'CALZETTI_MOD'``.
``uv_bump`` : flag
    If set, then a 2175 Angstrom UV bump feature will be added to the
    attenuation curve.
``error_check`` : flag
    If set, all inputs are checked for errors. Otherwise, all inputs are
    assumed to be of correct format.
``xray_isotropic`` : flag
    Deprecated (?). If set, the output ``L2500`` is assumed to be isotropic from the
    face-on model ``L2500``. Otherwise, it is assumed to be a function of viewing angle.

Output
------
``mean_Lnu_AGN`` : double array(Nfilters, Nmodels)
    The mean luminosity of each filter and set of AGN model parameters
    :math:`[L_\odot\ {\rm Hz}^{-1}]`.

Optional Outputs
----------------
``mean_Lnu_unred_AGN`` : double array(Nfilters, Nmodels)
    The unattenuated mean luminosity of each filter and set of AGN model
    parameters :math:`[L_\odot\ {\rm Hz}^{-1}]`.
``Lbol_abs_AGN``:  double array(Nmodels)
    The bolometric luminosity of the attenuated AGN emission for each
    set of model parameters, integrated for inclination variation :math:`[L_\odot]`.
``Lbol_AGN`` : double array(Nmodels)
    The bolometric luminosity of the AGN model for each set of model parameters,
    integrated for inclination variation :math:`[L_\odot]`.
``L2500`` : double array(Nmodels)
    The rest-frame 2500 Angstrom monochromatic luminosity of the intrinsic accretion
    disk spectrum shifted to the observed frame :math:`[L_\odot\ {\rm Hz}^{-1}]`.

Note
----
Defaults will only be set if the optional ``error_check`` input is set.

Reference
---------
`Stalevski, M., Ricci, C., Ueda, Y., et al. 2016, MNRAS, 458, 2288 <https://ui.adsabs.harvard.edu/abs/2016MNRAS.458.2288S/abstract>`_

Modification History
--------------------
- 2021/03/03: Created (Erik B. Monson)
- 2021/03/19: Added optional output for rest-frame ``L2500`` (Erik B. Monson)
- 2021/09/14: Interpolation now proceeds in log(Lnu) (Erik B. Monson)
- 2022/03/15: Documentation update (Erik B. Monson)
- 2022/06/09: Changed parameter names to standardized format and updated documentation (Keith Doore)
- 2022/06/09: Added error handling (Keith Doore)
- 2022/06/09: Added ``error_check`` keyword to do error handling (Keith Doore)
- 2022/06/15: Replaced ``!cv`` with ``!lightning_cgs`` (Keith Doore)
- 2022/06/30: Updated documentation and renamed variables (Keith Doore)
- 2022/06/30: Removed ``config`` and replaced with the two accessed configuration values (Keith Doore)
- 2022/06/30: Fixed potential bug of ``atten_curve`` not being set (Keith Doore)
- 2022/07/22: Added optional output of unreddened mean Lnu (Keith Doore)

