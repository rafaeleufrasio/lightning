SKIRTOR_MODELS
==============

Name
----
SKIRTOR_MODELS

Purpose
-------
Generates the spectra, SEDs, and AGN parameters for the SKIRTOR model 
(Stalevski et al. 2016) given a set of filters. The model is a slice 
of :math:`\tau`-:math:`\delta`-:math:`i` space, with constant values of
:math:`R`, :math:`p`, and :math:`q`.

Calling Sequence
----------------
::

    skirtor = skirtor_models([filter_labels = , redshift = , /error_check])

Optional Inputs
---------------
``filter_labels`` : string array(Nfilters)
    The filters labels for which models should be generated.
    (Default = ``['GALEX_FUV', 'GALEX_NUV', 'SDSS_u', 'SDSS_g', 'SDSS_r',
    'SDSS_i', 'SDSS_z', '2MASS_J', '2MASS_H', '2MASS_Ks', 'IRAC_CH1',
    'IRAC_CH2', 'IRAC_CH3', 'IRAC_CH4', 'MIPS_CH1', 'PACS_green',
    'PACS_red', 'SPIRE_250', 'SPIRE_350', 'SPIRE_500']``)
``redshift`` : int, float, or double scalar
    The redshift of the model. (Default = ``0.0``)
``error_check`` : flag
    If set, all inputs are checked for errors. Otherwise, all inputs are
    assumed to be of correct format.

Output
------
``skirtor`` : structure
    A structure including the full-resolution SKIRTOR models grid (in terms
    of :math:`L_\odot\ {\rm Hz}^{-1}`) and those convolved with the filters specified
    by ``filter_labels``.
    The full description of the structure is as follows:

    ===============     ============================    =================================================================================================
    TAG                 TYPE                            DESCRIPTION
    ===============     ============================    =================================================================================================
    FILTER_LABELS       string(Nfilters)                Filters labels, same as input
    WAVE_FILTERS        double(Nfilters)                Mean wavelength of the filters :math:`[\mu \rm m]`
    MEAN_LNU            double(Nfilters, Ninc, Ntau)    Model spectra convolved with the filters :math:`[L_\odot\ {\rm Hz}^{-1}]`
    WAVE_REST           double(Nwave)                   Rest-frame wavelength of the high-res models :math:`[\mu \rm m]`
    WAVE_OBS            double(Nwave)                   Observed-frame wavelength of the high-res models :math:`[\mu \rm m]`
    LNU_TOTAL           double(Nwave, Ninc, Ntau)       Total model spectra after reprocessing of disk emission by torus :math:`[L_\odot\ {\rm Hz}^{-1}]`
    LNU_TRANSPARENT     double(Nwave, Ninc, Ntau)       Model accretion disk spectra prior to reprocessing :math:`[L_\odot\ {\rm Hz}^{-1}]`
    LNU_INTEGRATED      double(Nwave, Ninc, Ntau)       Total inclination integrated model spectra :math:`[L_\odot\ {\rm Hz}^{-1}]`
    DUST_MASS_TOTAL     double(Ninc, Ntau)              Total dust mass of each model :math:`[M_\odot]`
    L2500A_ZERO         double                          Intrinsic accretion disk L2500 at ``i=0`` :math:`[L_\odot\ {\rm Hz}^{-1}]`
    FILTERS             double(Nfilters, Nwave)         Tabulated filter transmission functions, frequency normalized
    REDSHIFT            double                          Redshift of model, same as input
    ===============     ============================    =================================================================================================

Reference
---------
`Stalevski, M., Ricci, C., Ueda, Y., et al. 2016, MNRAS, 458, 2288 <https://ui.adsabs.harvard.edu/abs/2016MNRAS.458.2288S/abstract>`_

Modification History
--------------------
- 2021/03/03: Created (Erik B. Monson)
- 2022/03/15: Documentation update (Erik B. Monson)
- 2022/05/17: Turned ``lightning_dir`` into system variable call (Keith Doore)
- 2022/05/17: Removed ``get_filters`` input of ``lightning_dir``, since folded into ``get_filters`` (Keith Doore)
- 2022/05/17: Added proper error handling (Keith Doore)
- 2022/05/17: Added ``error_check`` keyword to do error handling (Keith Doore)
- 2022/05/17: Changed function name from ``load_SKIRTOR_grid`` to ``skirtor_models`` (Keith Doore)
- 2022/05/18: Updated path to model files (Keith Doore)
- 2022/05/18: Changed model save file from ``.idl`` to ``.sav`` (Keith Doore)
- 2022/06/01: Updated to fix covering factor and changed read in file format (Erik B. Monson)
- 2022/06/09: Merged updates between versions with and without ``CF`` (Keith Doore)
- 2022/06/30: Updated documentation (Keith Doore)
- 2022/06/30: Updated default values (Keith Doore)
- 2022/08/01: Added ``/silent`` to ``mrdfits`` (Keith Doore)

