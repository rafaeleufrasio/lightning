LIGHTNING_INPUT
===============

Name
----
LIGHTNING_INPUT

Purpose
-------
Reads in the SED data that is to be fit by Lightning from a FITS file. The
input fluxes are converted to luminosities using the user supplied redshift
or luminosity distance. These luminosities are then saved to a ``.sav`` file
for each SED to be fit by Lightning.

Calling Sequence
----------------
::

    lightning_input, input_file_fits [, cosmology = , /xray_emission, xray_unit = , xray_unc = ]

Input
-----
``input_file_fits`` : string scalar
    The name (including path) to the FITS file containing the SED fluxes
    and distances (or redshifts) in a data table. The table must be in
    the first extension. (See :ref:`fits-format-label` for full details,
    required contents, and format of the data table.)

Optional Inputs
---------------
``cosmology`` : structure
    A structure containing the cosmology parameters ``H0``, ``OMEGA_M``, ``LAMBDA0``,
    ``Q0``, and ``K`` in individual tags.
    (Default = ``{H0: 70, OMEGA_M: 0.3, LAMBDA0: 0.7, Q0: -0.55, K: 0}``)
``xray_emission`` : flag
    If set, Lightning will search for and load the additional data products necessary to fit
    an X-ray emission model.
``xray_unit`` : string scalar
    The type of X-ray data to use in fitting. Current options are ``'COUNTS'`` and ``'FLUX'``.
    (Default = ``'COUNTS'``)
``xray_unc`` : string scalar
    The errors to assume if using X-ray count data. Current options are ``'SQRT'``, ``'GEHRELS'``,
    and ``'USER'``. (Default = ``'GEHRELS'``)

Output
------
A save file (``<input_file_dir>/lightning_output/input_sav_files/lightning_input_<sed_id>.sav``)
for each SED that contains the data converted to luminosities.
The data for each SED are stored in a structure with the following format:

=======================     =================     ====================================================================================================================================
TAG                         TYPE                  DESCRIPTION
=======================     =================     ====================================================================================================================================
SED_ID                      string                Unique SED identifier
LNU_OBS                     double(Nfilters)      Luminosities of the SED for each set of filters in terms of :math:`L_\nu` :math:`[\rm{L_{\odot}\ Hz^{-1}}]`
LNU_UNC                     double(Nfilters)      Uncertainties associated with the luminosities :math:`[\rm{L_{\odot}\ Hz^{-1}}]`
FILTER_LABELS               string(Nfilters)      Names of the filters associated with each luminosity
REDSHIFT                    double                Redshift of the SED (if ``LUMIN_DIST`` was specified this is set to ``0``)
LUMIN_DIST                  double                Luminosity distance of the SED as input or converted from the redshift :math:`[\rm{Mpc}]`
XRAY_LNU_OBS [1]_           double(Nxray)         X-ray luminosities for each bandpass in terms of :math:`L_\nu` :math:`[\rm{L_{\odot}\ Hz^{-1}}]`
XRAY_LNU_UNC [1]_           double(Nxray)         Uncertainties associated with the X-ray luminosities :math:`[\rm{L_{\odot}\ Hz^{-1}}]`
XRAY_BANDPASS [1]_ [2]_     double(2, Nxray)      Bandpasses of X-ray observations: first column contains the lower energy bound, second column contains the upper. :math:`[\rm{keV}]`
XRAY_EXPOSURE [2]_          double(Nxray)         Exposure times of X-ray observations, one per band :math:`[\rm{s}]`
NET_COUNTS [2]_             double(Nxray)         Net counts in each X-ray band :math:`[\rm{counts}]`
NET_COUNTS_UNC [2]_         double(Nxray)         User-supplied uncertainty in net counts, if any :math:`[\rm{counts}]`
GALACTIC_NH [1]_ [2]_       double                Galactic (i.e., Milky Way) HI column density along the line of sight :math:`[10^{20}\ \rm{cm}^{-2}]`
ARF_E_LO [2]_               double(Nchannels)     Lower energy bounds of each channel in the ARF :math:`[\rm{keV}]`
ARF_E_HI [2]_               double(Nchannels)     Upper energy bounds of each channel in the ARF :math:`[\rm{keV}]`
ARF_SPECRESP [2]_           double(Nchannels)     Spectral response of the ARF at each channel :math:`[\rm{cm^2}]`
=======================     =================     ====================================================================================================================================

.. [1] These tags are present if the X-ray module is used with ``XRAY_UNIT = 'FLUX'``
.. [2] These tags are present if the X-ray module is used with ``XRAY_UNIT = 'COUNTS'``

Modification History
--------------------
- 2022/03/23: Created (Keith Doore)
- 2022/04/25: Updated to allow for Xray ARF file path (Keith Doore)
- 2022/04/28: Added unique ``SED_ID`` enforcement (Keith Doore)
- 2022/05/18: Added optional cosmology input to allow for unique cosmologies in ``lumdist`` (Keith Doore)
- 2022/05/18: Added check to make sure xray files existed (Keith Doore)
- 2022/05/19: Replace ``n_elements(where())`` statements in logicals with ``total()`` (Keith Doore)
- 2022/06/27: Removed missing data from SED before saving to file (Keith Doore)
- 2022/06/27: Fixed issue where xray_spec/arf file strings need to be trimmed of padded blanks (Keith Doore)
- 2022/08/01: Added ``/silent`` to ``mrdfits`` (Keith Doore)
- 2022/08/02: Undid removal of missing data from SED before saving to file to make post-processing simpler (Keith Doore)
- 2022/08/02: Added ability to use ``redshift`` if ``lumin_dist`` is ``0`` for some SEDs (Keith Doore)
- 2022/08/02: Added check to require both ``XRAY_SPEC_FILE`` and ``XRAY_ARF_FILE`` tags if including X-ray data (Keith Doore)
- 2022/08/09: Added ``GALACTIC_NH`` as input for xray emission (Keith Doore)
- 2022/08/10: Fixed bug with redshift distance being computed even with no redshifts (Keith Doore)
- 2022/09/01: Added handling for user-supplied X-ray count uncertainties (Erik B. Monson)
- 2022/09/14: Updated to allow fitting with X-ray fluxes (Erik B. Monson)
- 2022/09/15: Updated documentation (Keith Doore)
- 2022/10/24: Allowed for negative flux inputs (Keith Doore)

