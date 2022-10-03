ADD_FILTERS
===========

Name
----
ADD_FILTERS

Purpose
-------
Adds a user provided filter profile to Lightning. The filter profile can
either be input as a two column text file or array. It is then formatted
for use in Lightning. Additionally, the bandpass name, instrument and
observatory, or survey will need to be specified.

Calling Sequence
----------------
::

    add_filters, filter_profile, bandpass, instrument_survey [, observatory = , $
                 unit = , /frequency, _extra=_extra]

Inputs
------
``filter_profile`` : int, float, or double array(2, Nwave) or string scalar
    If ``filter_profile`` is a string, then it contains the path and file name
    to the tabulated filter profile that is to be read. For both the tabulated
    filter profile and the directly input array, the first column must contain
    the grid of wavelengths (or frequency, see below), and the second column
    must contain the corresponding transmission profile.
``bandpass`` : string scalar
    The bandpass name to give the filter profile.
``instrument_survey`` : string scalar
    The instrument or survey associated with the filter profile.

Optional Inputs
---------------
``observatory`` : string scalar
    The observatory associated with the filter profile. Omit if
    the filter is intended to be associated with a survey. (See
    Note below.)
``unit`` : int, float, or double scalar
    The unit conversion factor needed to convert the input wavelengths
    (frequency) associated with the filter profile to microns (Hertz).
    For example, if the input wavelengths are in nanometers, ``unit``
    should be set to ``1d-3``. (Default = ``1``)
``frequency`` : flag
    If set, the input or read-in filter profile is gridded in frequency
    rather than wavelength.
``_extra`` : structure
    Additional optional inputs that are passed to ``readcol.pro`` from the
    NASA Astro library (https://idlastro.gsfc.nasa.gov/ftp/pro/misc/readcol.pro).

Output
------
A text file containing the formatted filter profile that is saved to the users
local Lightning Filters directory. Additionally, the Lightning filter label
is printed to the screen.

Notes
-----
The output file and its path within the Lightning Filters directory is named
using the scheme: ``Filters/observatory/instrument/instrument_band.txt``. In the
case of survey specific filters (e.g., 2MASS, SDSS, etc.), the pattern is
changed to ``Filters/survey/survey_band.txt``. 

Modification History
--------------------
- 2022/07/28: Created (Keith Doore)

