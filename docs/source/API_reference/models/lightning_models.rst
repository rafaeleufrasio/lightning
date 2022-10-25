LIGHTNING_MODELS
================

Name
----
LIGHTNING_MODELS

Purpose
-------
Calls the functions to generate each model structure (stellar, dust, AGN,
and/or X-ray) specified in the Lightning configuration structure. The
structures are then placed in a parent structure for ease in passing to
other Lightning functions.

Calling Sequence
----------------
::

    models = lightning_models(config [, filter_labels = , redshift = , lumin_dist = , $
                              xray_bandpass = , xray_exposure = , arf_E_lo = , arf_E_hi = , $
                              arf_specresp = , galactic_nh = , input_dir = , /error_check])

Input
-----
``config`` : structure
    A Lightning configuration structure. (See
    ``lightning_configure_defaults.pro`` for details and contents.)

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
``lumin_dist`` : int, float, double scalar
    The luminosity distance of the model :math:`[\rm Mpc]`. (Default = ``10``)
``xray_bandpass`` : int, float or double array(2, Nxray)
    The bandpasses for the X-ray spectrum. The first column should be the lower
    energy bound, and the second should be the upper energy bound :math:`[\rm keV]`.
``xray_exposure`` : int, float or double array(Nxray)
    The exposure time of the observations, one per band :math:`[\rm s]`.
``arf_E_lo`` : float or double array(Nchannels)
    Lower energy bounds of each channel in the ARF :math:`[\rm keV]`.
``arf_E_hi`` : float or double array(Nchannels)
    Upper energy bounds of each channel in the ARF :math:`[\rm keV]`.
``arf_specresp`` : float or double array(Nchannels)
    The spectral response of the ARF at each channel :math:`[\rm cm^2]`.
``galactic_nH`` : int, float, or double scalar
    Galactic, i.e. Milky Way, neutral Hydrogen column density along the line of
    sight :math:`[10^{20}\ \rm{cm}^{-2}]`.
``input_dir`` : string scalar
    The path to the file containing the input SED data.
``error_check`` : flag
    If set, all inputs are checked for errors. Otherwise, all inputs are
    assumed to be of correct format.

Output
------
``models`` : structure
    A structure containing the filter labels, and each model structure
    (stellar, dust, AGN, and/or X-ray) as a substructures. If a model
    component is not specified in the Lighting configuration structure,
    then the value of the substructure will be set to ``NaN``.
    The full description of the structure is as follows:

    ==============     ================     =====================================================================================
    TAG                TYPE                 DESCRIPTION
    ==============     ================     =====================================================================================
    FILTER_LABELS      string(Nfilters)     Labels for filters, same as input
    STELLAR_MODELS     structure            Stellar emission models (See ``binned_stellar_models.pro`` for details and contents.)
    ATTEN_MODELS       structure            The preloaded files for the Doore+21 attenuation.
    DUST_MODELS        structure            Dust emission models (See ``dl07_models.pro`` for details and contents.)
    XRAY_MODELS        structure            X-ray emission models (See ``xrb_xagn_models.pro`` for details and contents.)
    AGN_MODELS         structure            AGN emission models (See ``skirtor_models.pro`` for details and contents.)
    ==============     ================     =====================================================================================

Notes
-----
- When using an X-ray emission model with ``XRAY_UNIT='COUNTS'``, the optional inputs ``xray_bandpass``,
  ``xray_exposure``, ``arf_E_lo``, ``arf_E_hi``, and ``arf_specresp`` become
  required inputs.
- When using an X-ray emission model with ``XRAY_UNIT='FLUX'``, the optional inputs ``xray_bandpass``,
  becomes a required input.
- When using the ``DOORE21`` attenuation curves, the optional input ``input_dir``
  becomes a required input.

Modification History
--------------------
- 2022/07/11: Created (Keith Doore)
- 2022/08/01: Added ``/silent`` to ``mrdfits`` (Keith Doore)
- 2022/08/09: Added ``GALACTIC_NH`` as input rather than keyword inheritance from ``config``, since now in input file (Keith Doore)
- 2022/10/25: Renamed SPS to SSP (Keith Doore)

