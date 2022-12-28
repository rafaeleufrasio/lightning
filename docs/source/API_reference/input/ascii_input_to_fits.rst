ASCII_INPUT_TO_FITS
===================

Name
----
ASCII_INPUT_TO_FITS

Purpose
-------
Reads in the SED data that is to be fit by Lightning from an ASCII file. The 
read-in data is then placed in a formatted structure and saved to the first
extension of a FITS file with the same name as the input file.

Calling Sequence
----------------
::

    ascii_input_to_fits, input_file [, input_file_fits=input_file_fits]

Input
-----
``input_file`` : string scalar
    The name (including path) to the ASCII file containing the SED fluxes 
    and distances (or redshifts) in a data table. (See :ref:`ascii-format-label`
    for full details, required contents, and format of the ASCII data table.)

Output
------
A FITS file (``<input_file_dir>/<input_file_name>.fits``) containing the
SED fluxes and distances (or redshifts) as read in from the ASCII data
table. (See :ref:`fits-format-label` for full details and format of the
FITS data table.)

Optional Output
---------------
``input_file_fits`` : string scalar
    The name (including path) of the output FITS file.

Note
----
No error handling is performed on the data in the table as this is done
with the converted FITS file data in ``lightning_input.pro``.

Modification History
--------------------
- 2022/08/02: Created (Keith Doore)
- 2022/09/19: Allowed for X-ray fluxes to be input (Keith Doore)
- 2022/12/28: Fixed bug if only reading in one SED (Keith Doore)
- 2022/12/28: Fixed bug if no x-ray data (Keith Doore)

