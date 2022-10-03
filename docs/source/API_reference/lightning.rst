LIGHTNING
=========

Name
----
LIGHTNING

Purpose
-------
Calls lower level Lightning procedures and functions to

1) generate, check, and save the configuration file;
2) load, check, and save the SED data;
3) call the Lightning fitting procedure to fit each SED using the specified configuration;
4) and post-process the fitting results.

Each call to the Lightning fitting procedure is optimized to run
in parallel depending on the allowed CPU usage.

Calling Sequence
----------------
::

    lightning, input_file [, /resume, /interactive]

Input
-----
``input_file`` : string scalar
    The name (including path) to the file containing the SED fluxes
    and distances (or redshifts) in a data table. If the file is a
    FITS file, the table must be in the first extension. (See
    :ref:`input-formats-label` for details, required contents and
    format of the data table.)

Optional Inputs
---------------
``resume`` : flag
    If set, then Lightning will resume running where it left off assuming
    some issue caused it to stop mid-run. This means Lightning would check
    what SEDs have been fit and would only fit those that were not completed.
``interactive`` : flag
    If set, then Lightning will give prompts to the command line allowing
    for a user to set the Lightning configuration interactively, rather
    than having to manually edit the configuration file.

Modification History
--------------------
- 2022/05/04: Created (Keith Doore)
- 2022/05/18: Rearranged input and configuration to allow for cosmology in input (Keith Doore)
- 2022/06/06: Moved tabulated prior check to its own procedure (Keith Doore)
- 2022/08/01: Added ``/silent`` to ``mrdfits`` (Keith Doore)
- 2022/08/02: Added ASCII input to FITS file conversion (Keith Doore)
- 2022/08/18: Added progress printing (Keith Doore)
- 2022/09/01: Swapped ``input_dir`` to absolute path, and passed ``config`` to ``lightning_input`` (Erik B. Monson)
- 2022/09/14: Updates to allow fitting with X-ray fluxes (Erik B. Monson)
- 2022/09/16: Added check for syntax errors in user edited ``lightning_configure.pro`` (Keith Doore)
- 2022/09/26: Fixed ``input_dir`` to absolute path bug (Keith Doore)

