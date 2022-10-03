DL07_SED
========

Name
----
DL07_SED

Purpose
-------
Generates the dust model SED for a given set of dust model parameters and
filters. The dust model is based on the Draine & Li (2007) dust model, which
include the parameters of :math:`\alpha`, :math:`U_{\rm min}`, :math:`U_{\rm max}`,
:math:`\gamma`, and :math:`q_{\rm PAH}`.

Calling Sequence
----------------
::

    mean_Lnu_dust = dl07_sed(dl07 [, alpha = , umin = , umax = , gam = , qPAH = , $
                             /error_check, LTIR=LTIR, pow_LTIR=pow_LTIR, $
                             del_LTIR=del_LTIR, pow_mean_Lnu=pow_mean_Lnu, $
                             del_mean_Lnu=del_mean_Lnu])

Input
-----
``dl07`` : structure
    A structure containing the spectra, SEDs, and dust parameters for 
    the DL07 model. (See ``dl07_models.pro`` for details and contents.)

Optional Inputs
---------------
``alpha`` : int, float, or double array(Nmodels)
    The exponent of the power-law distribution of heating starlight
    intensities between ``umin`` and ``umax``. (Default = ``2.0``)
``umin`` : int, float, or double array(Nmodels)
    The minimum radiation field intensity of the diffuse ISM radiation
    field heating the dust. (Default = ``1.0``)
``umax`` : int, float, or double array(Nmodels)
    The maximum radiation field intensity of the power-law distribution
    of heating starlight intensities. (Default = ``3e5``)
``gam`` : int, float, or double array(Nmodels)
    The fraction of the dust mass exposed to the power-law distribution of
    starlight intensities.  (Default = ``0.1``)
``qPAH`` : float or double array(Nmodels)
    The PAH index.  (Default = ``0.025``)
``error_check`` : flag
    If set, all inputs are checked for errors. Otherwise, all inputs are
    assumed to be of correct format.

Output
------
``mean_Lnu_dust`` : double array(Nfilters, Nmodels)
    The mean luminosity of each filter and set of dust model parameters
    :math:`[L_\odot\ {\rm Hz}^{-1}]`.

Optional Outputs
----------------
``LTIR`` : double array(Nmodels)
    The total integrated IR luminosity (i.e., the bolometric luminosity)
    of the dust model for each set of model parameters :math:`[L_\odot]`.
``pow_LTIR`` : double array(Nmodels)
    The power-law component of ``LTIR`` :math:`[L_\odot]`.
``del_LTIR`` : double array(Nmodels)
    The delta-function component of ``LTIR`` :math:`[L_\odot]`.
``pow_mean_Lnu`` : double array(Nfilters, Nmodels)
    The power-law component of ``mean_Lnu_dust`` :math:`[L_\odot\ {\rm Hz}^{-1}]`.
``del_mean_Lnu`` : double array(Nfilters, Nmodels)
    The delta-function component of ``mean_Lnu_dust`` :math:`[L_\odot\ {\rm Hz}^{-1}]`.

Note
----
Defaults will only be set if the optional ``error_check`` input is set.

Reference
---------
`Draine, B. T., & Li, A. 2007, ApJ, 657, 810 <https://ui.adsabs.harvard.edu/abs/2007ApJ...657..810D/abstract>`_

Modification History
--------------------
- 2020/03/01: Created (Rafael T. Eufrasio)
- 2020/09/30: Added ``gamma`` and ``del`` component (Rafael T. Eufrasio and Keith Doore)
- 2022/03/15: Added proper error handling (Keith Doore)
- 2022/03/18: Updated documentation (Keith Doore)
- 2022/04/12: Allowed for inputs to have degenerate dimensions (Keith Doore)
- 2022/04/12: Allowed for inputs to be scalars (Keith Doore)
- 2022/04/12: Allowed integer inputs (Keith Doore)
- 2022/04/12: Allowed for array of only one optional input to be given (Keith Doore)

