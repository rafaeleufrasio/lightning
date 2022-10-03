BINNED_STELLAR_MODELS
=====================

Name
----
BINNED_STELLAR_MODELS

Purpose
-------
Generates the spectra, SEDs, and stellar parameters for given a set of
filters and steps boundaries. The spectra and SEDs can include or not
include nebular absorption and emission lines, and they depend on the
optionally input metallicity and initial mass function (IMF).

Calling Sequence
----------------
::

    stellar_models = binned_stellar_models([filter_labels = , redshift = , $
                            steps_bounds = , dtime_SF = , Zmetal = , IMF = , $
                            cosmology = , /no_emission_lines, /no_nebular_extinction, $
                            /error_check])

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
``steps_bounds`` : int, float, or double array(Nsteps+1)
    The age bounds to separate each star formation step :math:`[{\rm yr}]`.
    (Default = ``[0.d0, 1.d7, 1.d8, 1.d9, 5.d9, 13.6d9]``)
``dtime_SF`` : int, float, or double scalar
    The time step used when interpolating the stellar population to the age bins
    defined by ``steps_bounds`` :math:`[{\rm yr}]`. (Default = ``5.d5``)
``Zmetal`` : float or double scalar
    The metallicity in terms of Z to be used when generating the stellar models.
    (Default = ``0.02``)
``IMF`` : string scalar
    The IMF to be used when generating the stellar models. The only available
    value is ``'Kroupa01'``. (Default = ``'Kroupa01'``)
``cosmology`` : structure
    A structure containing the cosmology parameters ``H0``, ``OMEGA_M``, ``LAMBDA0``,
    ``Q0``, and ``K`` in individual tags.
    (Default = ``{H0: 70, OMEGA_M: 0.3, LAMBDA0: 0.7, Q0: -0.55, K: 0}``)
``no_emission_lines`` : flag
    If set, no emission lines are included in the output spectra.
``no_nebular_extinction`` : flag
    If set, no nebular absorption is included in the output spectra.
``error_check`` : flag
    If set, all inputs are checked for errors. Otherwise, all inputs are
    assumed to be of correct format.

Output
------
``stellar_models`` : structure
    A structure including the full-resolution non-parametric stellar models
    (in terms of :math:`L_\odot\ {\rm Hz}^{-1}`) and those convolved with the filters
    specified by ``filter_labels``.
    The full description of the structure is as follows:

    =============     ========================     ========================================================================================
    TAG               TYPE                         DESCRIPTION
    =============     ========================     ========================================================================================
    FILTER_LABELS     string(Nfilters)             Labels for filters, same as input
    WAVE_FILTERS      double(Nfilters)             Mean wavelength of the filters :math:`[\mu \rm m]`
    MEAN_LNU          double(Nfilters, Nsteps)     Model spectra convolved with the filters :math:`[L_\odot\ {\rm Hz}^{-1}]`
    WAVE_REST         double(Nwave)                Rest-frame wavelength of the high-res models :math:`[\mu \rm m]`
    WAVE_OBS          double(Nwave)                Observed-frame wavelength of the high-res models :math:`[\mu \rm m]`
    LNU               double(Nwave, Nsteps)        Observed-frame Lnu spectrum emitted within each age bin :math:`[L_\odot\ {\rm Hz}^{-1}]`
    Q0                double(Nsteps)               Ionizing flux rate emitted within each age bin :math:`[{\rm photons\ s^{-1}}]`
    MSTAR             double(Nsteps)               Surviving stellar mass within each age bin :math:`[M_\odot]`
    LBOL              double(Nsteps)               Bolometric luminosity emitted within each age bin :math:`[L_\odot]`
    BOUNDS            double(Nsteps+1)             Age bounds to separate each star formation step :math:`[{\rm yr}]`
    FILTERS           double(Nfilters, Nwave)      Tabulated filter transmission functions, frequency normalized
    REDSHIFT          float/double                 Redshift of model, same as input
    =============     ========================     ========================================================================================

Modification History
--------------------
- 2016/05/01: Created (Rafael T. Eufrasio)
- 2020/03/09: Corrected (1+z) factor (Rafael T. Eufrasio)
- 2020/04/27: Replaced if statements using ``n_elements`` on keywords to use ``keyword_set`` (Keith Doore)
- 2022/03/15: Added proper error handling (Keith Doore)
- 2022/03/15: Renamed variables to standard format (Keith Doore)
- 2022/03/18: Updated documentation (Keith Doore)
- 2022/04/12: Allowed for ``filter_labels`` to be scalars (Keith Doore)
- 2022/04/12: Allowed inputs to be integers (Keith Doore)
- 2022/04/12: Replaced ``.length`` with ``n_elements`` (Keith Doore)
- 2022/04/13: Normalized Filters to the frequency grid (Keith Doore)
- 2022/05/16: Turned ``lightning_dir`` into system variable call (Keith Doore)
- 2022/05/16: Changed ``!cv`` constants call to ``!lightning_cgs`` call (Keith Doore)
- 2022/05/16: Removed ``Filters`` output, since it is redundant (Keith Doore)
- 2022/05/16: Changed ``no_lines`` to ``emission_lines`` to match naming scheme and inverted logic accordingly (Keith Doore)
- 2022/05/16: Changed ``no_nebular`` to ``nebular_extinction`` to match naming scheme and inverted logic accordingly (Keith Doore)
- 2022/05/16: Changed ``Zmet`` to ``Zmetal`` to match naming scheme (Keith Doore)
- 2022/05/16: Added bin truncation based on redshift (Keith Doore)
- 2022/05/16: Changed model save file from ``.idl`` to ``.sav`` (Keith Doore)
- 2022/05/16: Removed ``get_filters`` input of ``lightning_dir``, since folded into ``get_filters`` (Keith Doore)
- 2022/05/17: Added ``error_check`` keyword to do error handling (Keith Doore)
- 2022/05/17: Changed function name from ``steps_stellar`` to ``binned_stellar_models`` (Keith Doore)
- 2022/05/18: Added optional ``cosmology`` input to allow for unique cosmologies in ``galage`` call (Keith Doore)
- 2022/05/18: Updated path to model files (Keith Doore)
- 2022/06/29: Changed ``emission_lines`` to back to ``no_emission_lines`` as to be default if not set (Keith Doore)
- 2022/06/29: Changed ``nebular_extinction`` to back to ``no_nebular_extinction`` as to be default if not set (Keith Doore)

