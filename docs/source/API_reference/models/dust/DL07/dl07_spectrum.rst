DL07_SPECTRUM
=============

Name
----
DL07_SPECTRUM

Purpose
-------
Generates the dust model spectra for a given set of dust model parameters.
The dust model is based on the Draine & Li (2007) dust model, which include
the parameters of :math:`\alpha`, :math:`U_{\rm min}`, :math:`U_{\rm max}`,
:math:`\gamma`, and :math:`q_{\rm PAH}`.

Calling Sequence
----------------
::

    Lnu_spec_dust = dl07_spectrum(dl07 [, alpha = , umin = , umax = , gam = , qPAH = , $
                                  /error_check, LTIR=LTIR, pow_LTIR=pow_LTIR, $
                                  del_LTIR=del_LTIR, pow_Lnu_spec=pow_Lnu_spec, $
                                  del_Lnu_spec=del_Lnu_spec])

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
``Lnu_spec_dust`` : double array(Nwave, Nmodels)
    The luminosity spectrum for each set of dust model parameters at
    each wavelength :math:`[L_\odot\ {\rm Hz}^{-1}]`.

Optional Outputs
----------------
``LTIR`` : double array(Nmodels)
    The total integrated IR luminosity (i.e., the bolometric luminosity)
    of the dust model for each set of model parameters :math:`[L_\odot]`.
``pow_LTIR`` : double array(Nmodels)
    The power-law component of ``LTIR`` :math:`[L_\odot]`
``del_LTIR`` : double array(Nmodels)
    The delta-function component of ``LTIR`` :math:`[L_\odot]`.
``pow_Lnu_spec`` : double array(Nwave, Nmodels)
    The power-law component of ``Lnu_spec_dust`` :math:`[L_\odot\ {\rm Hz}^{-1}]`.
``del_Lnu_spec`` : double array(Nwave, Nmodels)
    The delta-function component of ``Lnu_spec_dust`` :math:`[L_\odot\ {\rm Hz}^{-1}]`.

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
- 2022/06/30: Updated documentation and renamed variables (Keith Doore)

