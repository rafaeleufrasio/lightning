LIGHTNING_XRAY_INPUT
====================

Name
----
LIGHTNING_XRAY_INPUT

Purpose
-------
Reads in the necessary X-ray data from a spectral file (e.g., the outputs of
ACISExtract) and Auxiliary Response Function file.

Calling Sequence
----------------
::

    xray_data = lightning_xray_input(sed_id, xray_spec_file, xray_arf_file [, xray_unc =])

Inputs
------
``sed_id`` : string scalar
    A unique SED identifier
``xray_spec_file`` : string scalar
    The name (including path) to the FITS file containing the X-ray spectral
    data table, which must be in the first extension.
    The data table must have the following structure:

    ==============     ===================     ==============================================================
    TAG                TYPE                    DESCRIPTION
    ==============     ===================     ==============================================================
    ENERG_LO           float/double(Nxray)     Lower energy bound of each observation band :math:`[\rm{keV}]`
    ENERG_HI           float/double(Nxray)     Upper energy bound of each observation band :math:`[\rm{keV}]`
    NET_COUNTS         float/double(Nxray)     Net counts in each band :math:`[\rm{counts}]`
    NET_COUNTS_UNC     float/double(Nxray)     Optional, uncertainty on net counts :math:`[\rm{counts}]`
    EXPOSURE           float/double(Nxray)     Exposure time of each band :math:`[\rm{s}]`
    ==============     ===================     ==============================================================

``xray_arf_file`` : string scalar
    The name (including path) to the FITS file containing the ARF
    data table, which must be in the first extension.
    The data table must have the following structure:

    ========     =======================     ======================================================
    TAG          TYPE                        DESCRIPTION
    ========     =======================     ======================================================
    ENERG_LO     float/double(Nchannels)     Lower energy bounds of each channel :math:`[\rm{keV}]`
    ENERG_HI     float/double(Nchannels)     Upper energy bounds of each channel :math:`[\rm{keV}]`
    SPECRESP     float/double(Nchannels)     Spectral response at each channel :math:`[\rm{cm}^2]`
    ========     =======================     ======================================================

Optional Input
--------------
``xray_unc`` : string scalar
    The errors to assume if using X-ray count data. Current options are ``'SQRT'``, ``'GEHRELS'``,
    and ``'USER'``. (Default = ``'GEHRELS'``)

Output
------
``xray_data`` : structure
    This structure includes the X-ray data and associated ARF data.
    The full description of the structure is as follows:

    ==============     =================     ====================================================================================================================================
    TAGS               TYPE                  DESCRIPTION
    ==============     =================     ====================================================================================================================================
    XRAY_BANDPASS      double(2, Nxray)      Bandpasses of X-ray observations: first column contains the lower energy bound, second column contains the upper. :math:`[\rm{keV}]`
    XRAY_EXPOSURE      double(Nxray)         Exposure times of X-ray observations, one per band :math:`[\rm{s}]`
    NET_COUNTS         double(Nxray)         Net counts in each X-ray band :math:`[\rm{counts}]`
    NET_COUNTS_UNC     double(Nxray)         Uncertainty on the net counts in each X-ray band :math:`[\rm{counts}]`
    ARF_E_LO           double(Nchannels)     Lower energy bounds of each channel in the ARF :math:`[\rm{keV}]`
    ARF_E_HI           double(Nchannels)     Upper energy bounds of each channel in the ARF :math:`[\rm{keV}]`
    ARF_SPECRESP       double(Nchannels)     Spectral response of the ARF at each channel :math:`[\rm{cm}^2]`
    ==============     =================     ====================================================================================================================================

Modification History
--------------------
- 2022/06/08: Created (Erik B. Monson)
- 2022/06/08: Updated error handling (Keith Doore)
- 2022/09/01: Added handling for user-supplied X-ray count uncertainties (Erik B. Monson)
- 2022/09/01: Updated documentation (Keith Doore)

