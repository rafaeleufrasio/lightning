.. highlight:: text
.. _input-formats-label:

Input Formats
=============

.. _data-input-label:

Data Format
-----------

Below, we give a general description of the allowed data columns that can be within the input files.


SED IDs
^^^^^^^

An identifier (ID) assigned to each SED. An ID is a single string that can contain any character,
excluding spaces, tabs, and commas. Including the SED ID column is **optional** and excluding it
will result in sequential integers being used as the IDs instead. When including this column in
either of the :ref:`table-format-label`, the column name must be ``SED_ID``.


Fluxes and Filter Labels (UV-to-IR)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A set of fluxes (in terms of :math:`F_\nu`), which must be within the UV-to-IR wavelength range, are
**required** inputs for each SED. Additionally, the :math:`1\sigma` uncertainty on the fluxes and the corresponding
names (labels) of the filters are **required** for each SED. Both the fluxes and the uncertainties are
the calibrated observed values in units of :math:`\rm Jy`. The filter labels are strings indicating which
instrument and filter corresponds to each observation. (A list of filter labels understood by Lightning can
be found :ref:`here <filter_table>`.) The names and format of the fluxes and filter label columns are unique
to each of the :ref:`table-format-label` and are described in the respective sections below.

.. note::

    - If multiple SEDs are included in the input file and some of these SEDs do not have an observation in a
      given filter, setting the corresponding uncertainty to ``0`` will result in that filter being ignored during
      fitting.
    - For upper limits on fluxes, we recommend setting the flux to ``0`` and the corresponding uncertainty to 
      the :math:`1\sigma` flux upper limit. This is not the proper way to account for upper limits, which 
      requires a complex adjustment to the computation of :math:`\chi^2` (see Appendix A of `Sawicki 2012
      <https://ui.adsabs.harvard.edu/abs/2012PASP..124.1208S/abstract>`_ for details). However, it is a
      reasonable approximation.


.. _distance-input-label:

Distance Indicators
^^^^^^^^^^^^^^^^^^^

A distance indicator is a **required** input for each SED. This indicator can either be a luminosity distance
in :math:`{\rm Mpc}` or a redshift, from which a luminosity distance can be inferred. The column names in
both of the :ref:`table-format-label` are ``lumin_dist`` and ``redshift`` for the luminosity distance and
redshift, respectively. Only one of the columns is required. If both columns are included, the
luminosity distance takes precedence over the redshift column unless its value is equal to ``0``, in which
case the redshift will be used as the distance indicator instead.


.. _xray-input-label:

X-ray Data
^^^^^^^^^^

X-ray data for a set of SEDs can take one of two forms when input into Lightning. The first is in terms
of net counts (i.e., background subtracted). The second is in terms of flux (in units of
:math:`{\rm erg\ cm^{-2}\ s^{-1}}`). Each type has a different input format, and the format must match the ``XRAY_UNIT``
:ref:`configuration setting <configure-setting-label>`. Below, we describe the format for both X-ray data
forms and the corresponding input file columns. Additionally, both types of input **require** the
column ``galactic_nh``, which gives the value of the Galactic (i.e., Milky Way) HI column density
along the line of sight in :math:`10^{20}\ \rm{cm}^{-2}`.

Counts
""""""

For X-ray data in units of counts, the columns ``xray_spec_file`` and ``xray_arf_file``
are **required** inputs for each SED when using either of the :ref:`table-format-label`.
These columns must be a string containing the file name (including path) to the FITS-formatted X-ray
spectrum (e.g., the outputs of ACISExtract) and a string containing the file name (including path) to
the FITS-formatted X-ray `Auxiliary Response Function (ARF) <https://cxc.cfa.harvard.edu/ciao/dictionary/arf.html>`_,
respectively. The contents of the spectral and ARF files are given in the tables below.

**X-ray spectral file contents**:

==============     ===================     ==============================================================
TAG                TYPE                    DESCRIPTION
==============     ===================     ==============================================================
ENERG_LO           float/double(Nxray)     Lower energy bound of each observation band :math:`[\rm{keV}]`
ENERG_HI           float/double(Nxray)     Upper energy bound of each observation band :math:`[\rm{keV}]`
NET_COUNTS         float/double(Nxray)     Net counts in each band :math:`[\rm{counts}]`
NET_COUNTS_UNC     float/double(Nxray)     Optional, uncertainty on net counts :math:`[\rm{counts}]`
EXPOSURE           float/double(Nxray)     Exposure time of each band :math:`[\rm{s}]`
==============     ===================     ==============================================================

.. note::

    The ``NET_COUNTS_UNC`` tag is optional. It only needs to be provided if using user input
    count uncertainties (i.e., :ref:`configuration setting <configure-setting-label>` ``XRAY_UNC = 'USER'``)

**ARF file contents**:

========     =======================     ======================================================
TAG          TYPE                        DESCRIPTION
========     =======================     ======================================================
ENERG_LO     float/double(Nchannels)     Lower energy bounds of each channel :math:`[\rm{keV}]`
ENERG_HI     float/double(Nchannels)     Upper energy bounds of each channel :math:`[\rm{keV}]`
SPECRESP     float/double(Nchannels)     Spectral response at each channel :math:`[\rm{cm}^2]`
========     =======================     ======================================================


Flux
""""

For X-ray data in units of flux, they are input in a similar style as the UV-to-IR fluxes.
The **required** inputs for each SED are the flux(es) (in terms of :math:`F`, the integrated flux over the bandpass),
the :math:`1\sigma` uncertainty on the flux(es), and the corresponding X-ray bandpass(es).
Both the fluxes and the uncertainties are in units of :math:`{\rm erg\ cm^{-2}\ s^{-1}}`.
The X-ray bandpasses are the lower and upper energy of each observation band in :math:`\rm keV`.
The names and format of the fluxes and bandpass columns are unique
to each of the :ref:`table-format-label` and are described in the respective sections below.



.. _table-format-label:

Table Formats
-------------

Lightning expects the input data to be input as either as either an ASCII or FITS data table.
The :ref:`ascii-format-label` and :ref:`fits-format-label`, along with their allowed column
values, are described below.

.. warning::

    Failure to follow these file formats can result in Lightning producing unknown errors.
    Please make sure your data is formatted as described!


.. _ascii-format-label:

ASCII Table Format
^^^^^^^^^^^^^^^^^^

When inputting your data as an ASCII table, the file can contain an unlimited number of comment
lines at the top of the file, where ``#`` is the comment character. Additionally, the final
comment line is **required** to contain the column names. The names of the columns and their
corresponding data can be separated by either a space, tab, or comma. Additionally, the data values
in the table can be in either decimal or scientific notation form.

.. warning::

    Do not leave any blank values in a column. Values that are to be treated as missing or
    blank should have their value set to ``NaN`` instead.


Here is an example of what an ASCII table for use with Lightning could look like::

    # This is an example of the ASCII table format for input into Lightning.
    # sed_id             lumin_dist  redshift  filterA     filterA_unc  filterB    filterB_unc  ...
    NGC_5194             8.2         NaN       0.568       0.011        5.77       0.29         ...
    J123624.82+620719.2  0           0.1141    3.0747e-05  1.08e-07     1.3232e-4  1.57e-07     ...

The first three columns are the SED IDs, the luminosity distances, and redshifts, respectively.
The contents and format of these columns are described :ref:`above <data-input-label>`. The
remaining columns are the UV-to-IR fluxes and :math:`1\sigma` uncertainties.
The names of these columns indicate the associated filter label. For example, if your
first filter was ``SDSS_u``, you would replace the column names ``filterA`` and ``filterA_unc``
with ``SDSS_u`` and ``SDSS_u_unc``, respectively.

To give an example of each type of X-ray input, we will expand on the example above. For the
X-ray counts input method, the ASCII table should look like::

    # This is an example of the ASCII table format for X-ray data input into Lightning.
    # sed_id             lumin_dist  redshift  filterA     filterA_unc  filterB    filterB_unc  galactic_nh  xray_spec_file                                     xray_arf_file
    NGC_5194             8.2         NaN       0.568       0.011        5.77       0.29         1.53         <path_to_file>/NGC_5194_xray_spec.fits             <path_to_file>/NGC_5194.arf
    J123624.82+620719.2  0           0.1141    3.0747e-05  1.08e-07     1.3232e-4  1.57e-07     1.48         <path_to_file>/J123624.82+620719.2_xray_spec.fits  <path_to_file>/J123624.82+620719.2.arf

The three newly added columns give the Galactic HI column density, the X-ray spectral file, and the
X-ray ARF file as described :ref:`above <xray-input-label>`.

For the X-ray flux input method, the ``xray_spec_file`` and ``xray_arf_file`` columns will need to be
replaced with the X-ray bandpass, flux, and flux uncertainty columns. These columns are formatted
similarly to the UV-to-IR flux columns, where the ending of the column name relates the bandpass to the
corresponding flux. Updating our example to include these flux and bandpass columns, our ASCII table
should look like::

    # This is an example of the ASCII table format for X-ray data input into Lightning.
    # sed_id             lumin_dist  redshift  filterA     filterA_unc  filterB    filterB_unc  galactic_nh  xray_bandpass_l_1   xray_bandpass_u_1  xray_flux_1  xray_flux_unc_1  xray_bandpass_l_2   xray_bandpass_u_2  xray_flux_2  xray_flux_unc_2
    NGC_5194             8.2         NaN       0.568       0.011        5.77       0.29         1.53         0.5                 2.0                3.24E-02     1.54e-03         2.0                 7.0                1.83E-02     2.70e-03
    J123624.82+620719.2  0           0.1141    3.0747e-05  1.08e-07     1.3232e-4  1.57e-07     1.48         0.5                 7.0                1.47e-08     2.31e-09         NaN                 NaN                NaN          NaN

Notice that each bandpass and flux are related to each other with the ending of each column name (i.e.,
``xray_bandpass_l_1``, ``xray_bandpass_u_1``, ``xray_flux_1``, and ``xray_flux_unc_1`` are the first
X-ray bandpass indicated by the ``_1`` ending). This numbering can be increased arbitrarily for any
number of desired bandpasses, which allows for multiple bandpasses to be input for each SED.
Also, values describing a single X-ray bandpass are contained within the ``xray_bandpass_l_*``
and ``xray_bandpass_u_*`` columns, which give the lower and upper bounds of the bandpass in
:math:`{\rm keV}`, respectively. This allows for each bandpass to be unique to each SED. Finally,
if an SED has more bandpasses than another in the input catalogue, the one with less bandpasses should
have the columns of the unused numbered bandpasses set to ``NaN``, as shown in the example.

.. note::

    The order of the columns in the ASCII table does not matter. However, the column names must be
    those described above. Changing the column names (besides swapping in the appropriate filter label)
    will result in errors.



.. _fits-format-label:

FITS Table Format
^^^^^^^^^^^^^^^^^

When inputting your data as a FITS data table, the table is **required** to be in the first extension
of the FITS file. Below, we describe the format of the basic and X-ray columns separately.
However, the X-ray columns must be in the same FITS data table as the basic columns and are **required**
if using an X-ray emission model.

.. note::

    We define the array size variables here for convenience.

    - ``Nfilters`` : the number of unique filters included in the input.
    - ``Nsed``: the number of SEDs included in the input
    - ``Nxray``: the maximum number of X-ray bandpasses for an SED included in the input


**Basic Columns**:

=====================     ============================     ============================================================
Column Names              Type (Shape)                     Description
=====================     ============================     ============================================================
SED_ID                    string(Nsed)                     **Optional**, unique SED identifier
FNU_OBS                   float/double(Nfilters, Nsed)     Fluxes of each SED for each set of filters :math:`[\rm{Jy}]`
FNU_UNC                   float/double(Nfilters, Nsed)     Uncertainties associated with the fluxes :math:`[\rm{Jy}]`
FILTER_LABELS [1]_        string(Nfilters, Nsed)           Filters labels associated with each flux
REDSHIFT [2]_             int/float/double(Nsed)           Redshift of each SED
LUMIN_DIST [2]_           int/float/double(Nsed)           Luminosity distance of each SED :math:`[\rm{Mpc}]`
=====================     ============================     ============================================================

**X-ray Columns**:

=======================     ============================     ========================================================================================================================================================================
Column Names                Type (Shape)                     Description
=======================     ============================     ========================================================================================================================================================================
GALACTIC_NH                 int/float/double(Nsed)           Galactic (i.e., Milky Way) HI column density along the line of sight :math:`[10^{20}\ \rm{cm}^{-2}]`
XRAY_SPEC_FILE [3]_         string(Nsed)                     File name (including path) containing the FITS-formatted Xray spectrum
XRAY_ARF_FILE [3]_          string(Nsed)                     File name (including path) containing the Xray Auxiliary Response Function (ARF)
XRAY_BANDPASS [4]_ [5]_     float/double(2, Nxray, Nsed)     Bandpasses of X-ray observations: first index of first dimension contains the lower energy bound, second index of first dimension contains the upper. :math:`[\rm{keV}]`
XRAY_FLUX [4]_              float/double(Nxray, Nsed)        X-ray fluxes of each SED for each set of bandpasses :math:`[{\rm erg\ cm^{-2}\ s^{-1}}]`
XRAY_FLUX_UNC  [4]_         float/double(Nxray, Nsed)        Uncertainties associated with the X-ray fluxes :math:`[{\rm erg\ cm^{-2}\ s^{-1}}]`
=======================     ============================     ========================================================================================================================================================================

.. rubric:: Table Notes

.. [1] ``FILTER_LABELS`` must be a 2-D array, where the first dimension holds each unique filter label and the second dimension is
   the first dimension repeated ``Nsed`` times.
.. [2] Only one of the columns is required as described :ref:`above <distance-input-label>`.
.. [3] Only used if inputting X-ray data in units of counts (i.e., :ref:`configuration setting <configure-setting-label>` ``XRAY_UNIT = 'COUNTS'``
.. [4] Only used if inputting X-ray data in units of flux (i.e., :ref:`configuration setting <configure-setting-label>` ``XRAY_UNIT = 'FLUX'``
.. [5] The X-ray bandpasses for each SED can be unique. If an SED has more bandpasses than another in the input catalogue, the one with
   less bandpasses should have the value of the unused bandpasses set to ``NaN``.
