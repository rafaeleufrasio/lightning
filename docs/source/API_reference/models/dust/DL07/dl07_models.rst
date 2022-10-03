DL07_MODELS
===========

Name
----
DL07_MODELS

Purpose
-------
Generates the spectra, SEDs, and dust parameters for the Draine & Li (2007) 
dust models given a set of filters. The spectra and SEDs are tabulated as a 
function of dust model parameters.

Calling Sequence
----------------
::

    dl07 = dl07_models([filter_labels = , redshift = , /error_check])

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
``dl07`` : structure
    A structure including the full-resolution DL07 models grid (in terms
    of Lsun Hz-1) and those convolved with the filters specified
    by ``filter_labels``.
    The full description of the structure is as follows:

    =============     ===========================     ==================================================================================
    TAG               TYPE                            DESCRIPTION
    =============     ===========================     ==================================================================================
    FILTER_LABELS     string(Nfilters)                Filters labels, same as input
    WAVE_FILTERS      double(Nfilters)                Mean wavelength of the filters :math:`[\mu \rm m]`
    MEAN_LNU          double(Nfilters, NU, Nqpah)     Model spectra convolved with the filters :math:`[L_\odot\ {\rm Hz}^{-1}]`
    WAVE_REST         double(Nwave)                   Rest-frame wavelength of the high-res models :math:`[\mu \rm m]`
    WAVE_OBS          double(Nwave)                   Observed-frame wavelength of the high-res models :math:`[\mu \rm m]`
    LNU               double(Nwave, NU, Nqpah)        Observed-frame Lnu spectrum emitted by each model :math:`[L_\odot\ {\rm Hz}^{-1}]`
    U                 string(NU)                      Tabulated radiation field intensities
    MODEL             string(Nqpah)                   Milky Way models corresponding to each ``QPAH``
    QPAH              double(Nqpah)                   Tabulated PAH index values
    LBOL              double(NU, Nqpah)               Bolometric luminosity emitted by each model :math:`[L_\odot]`
    FILTERS           double(Nfilters, Nwave)         Tabulated filter transmission functions, frequency normalized
    REDSHIFT          float/double                    Redshift of model, same as input
    =============     ===========================     ==================================================================================

Reference
---------
`Draine, B. T., & Li, A. 2007, ApJ, 657, 810 <https://ui.adsabs.harvard.edu/abs/2007ApJ...657..810D/abstract>`_

Modification History
--------------------
- 2020/03/01: Created (Rafael T. Eufrasio)
- 2020/04/27: Made folder location ``lightning_dir`` vs fix location in code (Keith Doore)
- 2020/04/27: Changed ``folder`` keyword to ``lightning_dir`` for consistence between functions/procedures for ``_REF_EXTRA`` (Keith Doore)
- 2020/05/11: Corrected (1+z) term (Keith Doore)
- 2022/03/15: Added proper error handling (Keith Doore)
- 2022/03/15: Renamed variables to standard format (Keith Doore)
- 2022/03/18: Updated documentation (Keith Doore)
- 2022/04/12: Allowed for ``filter_labels`` to be scalars (Keith Doore)
- 2022/04/12: Changed ``readcol`` to skip extra beginning lines (Keith Doore)
- 2022/04/12: Replaced ``.length`` with ``n_elements`` (Keith Doore)
- 2022/05/16: Turned ``lightning_dir`` into system variable call (Keith Doore)
- 2022/05/16: Changed ``!cv`` constants call to ``!lightning_cgs`` call (Keith Doore)
- 2022/05/16: Removed ``get_filters`` input of ``lightning_dir``, since folded into ``get_filters`` (Keith Doore)
- 2022/05/17: Added ``error_check`` keyword to do error handling (Keith Doore)
- 2022/05/17: Changed function name from ``dl07_templates`` to ``dl07_models`` (Keith Doore)
- 2022/05/18: Updated path to model files (Keith Doore)
- 2022/06/29: Converted ``U`` in output to double from string (Keith Doore)
- 2022/06/29: Replaced for loops with ``reform`` for ``Lnu`` and ``mean_Lnu`` arrays (Keith Doore)

