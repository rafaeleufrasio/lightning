XRB_XAGN_MODELS
===============

Name
----
XRB_XAGN_MODELS

Purpose
-------
Generates the spectra, counts, and X-ray parameters for given a set of
X-ray bandpasses and ARF for an X-ray binary population (XRB). The
spectra and counts can include or not include Galactic absorption
and X-ray AGN emission.

Calling Sequence
----------------
::

    xray_models = xrb_xagn_models(xray_bandpass [, xray_exposure = , arf_E_lo = , arf_E_hi = , $
                                  arf_specresp = , redshift = , lumin_dist = , $
                                  galactic_nH = , xray_abs_model = , $
                                  xray_agn_model = , xray_unit = , /error_check])

Input
-----
``xray_bandpass`` : int, float or double array(2, Nxray)
    The bandpasses for the X-ray spectrum. The first column should be the lower
    energy bound, and the second should be the upper energy bound :math:`[{\rm keV}]`.

Optional Inputs
---------------
``xray_exposure`` : int, float or double array(Nxray)
    The exposure time of the observations, one per band :math:`[{\rm s}]`.
    Required to generate model count-rate spectrum.
``arf_E_lo`` : float or double array(Nchannels)
    Lower energy bounds of each channel in the ARF :math:`[{\rm keV}]`.
    Required to generate model count-rate spectrum.
``arf_E_hi`` : float or double array(Nchannels)
    Upper energy bounds of each channel in the ARF :math:`[{\rm keV}]`.
    Required to generate model count-rate spectrum.
``arf_specresp`` : float or double array(Nchannels)
    The spectral response of the ARF at each channel :math:`[{\rm cm}^2]`.
    Required to generate model count-rate spectrum.
``redshift`` : int, float, or double scalar
    The redshift of the model. (Default = ``0.0``)
``lumin_dist`` : int, float, double scalar
    The luminosity distance of the model :math:`[{\rm Mpc}]`. (Default = ``10``)
``galactic_nH`` : int, float, or double scalar
    Galactic, i.e. Milky Way, neutral Hydrogen column density along the line of sight
    :math:`[10^{20}\ {\rm cm}^{-2}]`. (Default = ``0``)
``xray_abs_model`` : string scalar
    The name of the X-ray absorption model to apply to the X-ray emission. Current options
    are: ``'TBABS-WILM'``, ``'ATTEN'``, and ``'NONE'``.
    (Default = ``'TBABS-WILM'``)
``xray_agn_model`` : string scalar
    The X-ray AGN emission model to use. Current options are: ``'PLAW'``, ``'QSOSED'``, ``'NONE'``.
    (Default = ``'QSOSED'``)
``xray_unit`` : string scalar
    The type of X-ray data to use in fitting. Current options are ``'COUNTS'`` and ``'FLUX'``.
    (Default = ``'COUNTS'``)
``error_check`` : flag
    If set, all inputs are checked for errors. Otherwise, all inputs are
    assumed to be of correct format.

Output
------
``xray_models`` : structure
    A structure including the the XRB, hot gas, and X-ray AGN models in terms
    of luminosities and count rates given the input ARF.
    The full description of the structure is as follows:

    ======================     ===========================     ================================================================================================================================================
    TAG                        TYPE                            DESCRIPTION
    ======================     ===========================     ================================================================================================================================================
    XRAY_BANDPASS              double(2, Nxray)                X-ray bandpasses :math:`[{\rm keV}]`, same as input
    XRAY_EXPOSURE [1]_         double(Nxray)                   X-ray exposure time per band :math:`[{\rm s}]`, same as input
    WAVE_REST                  double(Nwave)                   Rest-frame wavelength of the models :math:`[\mu \rm m]`
    WAVE_OBS                   double(Nwave)                   Observed-frame wavelength of the models :math:`[\mu \rm m]`
    EXP_NEG_TAU_XRAY           double(Nwave)                   Rest-frame absorption curve normalized to :math:`{\rm nH} = 10^{20}\ {\rm cm}^2`
    EXP_NEG_TAU_XRAY_MW        double(Nwave)                   Observed-frame absorption curve normalized to :math:`{\rm nH} = 10^{20}\ {\rm cm}^2`
    LNU_AGN [2]_               double(Nwave, Nmass, Nmdot)     AGN model normalized to 1 :math:`L_\odot\ {\rm Hz}^{-1}` at rest-frame 2 keV :math:`[L_\odot\ {\rm Hz}^{-1}]`
    LNU_XRB                    double(Nwave)                   XRB model normalized to 1 :math:`L_\odot\ {\rm Hz}^{-1}` over rest-frame 2-10 keV bandpass :math:`[L_\odot\ {\rm Hz}^{-1}]`
    LNU_GAS                    double(Nwave)                   Hot gas model normalized to 1 :math:`L_\odot\ {\rm Hz}^{-1}` over rest-frame 0.5-2 keV bandpass :math:`[L_\odot\ {\rm Hz}^{-1}]`
    AGN_MODEL [1]_ [2]_        double(Nwave, Nmass, Nmdot)     Instrumental count rate density produced by normalized AGN model :math:`[{\rm counts\ s^{-1}\ Hz^{-1}}]`
    XRB_MODEL [1]_             double(Nwave)                   Instrumental count rate density produced by normalized XRB model :math:`[{\rm counts\ s^{-1}\ Hz^{-1}}]`
    GAS_MODEL [1]_             double(Nwave)                   Instrumental count rate density produced by normalized hot gas model :math:`[{\rm counts\ s^{-1}\ Hz^{-1}}]`
    XRB_COUNTS_LNU [1]_        double(Nxray)                   Instrumental counts produced by normalized XRB model integrated over the supplied bandpass (Galactic absorption only) :math:`[{\rm counts}]`
    GAS_COUNTS_LNU [1]_        double(Nxray)                   Instrumental counts produced by normalized hot gas model integrated over the supplied bandpass (Galactic absorption only) :math:`[{\rm counts}]`
    GALACTIC_NH                double                          Galactic neutral hydrogen column density :math:`[10^{20}\ {\rm cm}^{-2}]`, same as input
    L2500 [3]_                 double(Nmass, Nmdot)            2500 Angstrom intrinsic luminosity grid from qsosed model, shifted to the observed frame :math:`[L_\odot\ {\rm Hz}^{-1}]`.
    AGN_MASS [3]_              double(Nmass, Nmdot)            Supermassive black hole mass grid from qsosed model :math:`[M_\odot]`
    AGN_LOGMDOT [3]_           double(Nmass, Nmdot)            Log10 of SMBH accretion rate grid from qsosed model, normalized by the Eddington rate
    REDSHIFT                   double                          Redshift of the model, same as input
    XRAY_ABS_MODEL             string                          Name of the X-ray absorption model, same as input
    XRAY_AGN_MODEL             string                          Name of the AGN model, same as input
    ======================     ===========================     ================================================================================================================================================

.. [1] ``XRAY_EXPOSURE``, ``AGN_MODEL``, ``XRB_MODEL``, ``GAS_MODEL``, ``XRB_COUNTS_LNU``, and ``GAS_COUNTS_LNU``
   are set to ``NaN`` if the ``xray_exposure`` and ``arf_*`` keywords are not specified.
.. [2] ``LNU_AGN`` and ``AGN_MODEL`` will only have the second and third dimensions if
   using the qsosed model. Additionally, the qsosed model is not normalized to 1 :math:`L_\odot\ {\rm Hz}^{-1}`.
.. [3] ``L2500``, ``AGN_MASS``, and ``AGN_LOGMDOT`` only have values if using the qsosed
   model. Otherwise, they will be ``NaN``.

Modification History
--------------------
- 2021/09/21: Created (Erik B. Monson).
- 2021/10/11: Added X-ray absorption (Erik B. Monson).
- 2022/03/16: Moved to separate file and documentation improved (Erik B. Monson).
- 2022/04/18: qsosed option now loads entire grid (Erik B. Monson)
- 2022/06/22: Major update to include new implementation (Keith Doore)
- 2022/07/07: Changed several variable names and updated documentation (Keith Doore)
- 2022/09/14: Updates to allow fitting with X-ray fluxes (Erik B. Monson)
- 2022/11/02: Galactic NH is now in units of 1e20 cm-2 (Erik B. Monson)

