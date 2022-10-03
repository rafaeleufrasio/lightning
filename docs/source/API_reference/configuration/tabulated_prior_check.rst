TABULATED_PRIOR_CHECK
=====================

Name
----
TABULATED_PRIOR_CHECK

Purpose
-------
Checks that all tabulated priors specified in the Lightning 
configuration structure are properly formatted. The correct format
requires that each specified tabulated prior is found within the 
tabulated prior files for each SED, has gridded values in ascending 
order, and non-negative probabilities for each grid value.

Calling Sequence
----------------
::

    tabulated_prior_check, config, sed_id

Inputs
------
``config`` : structure
    A Lightning configuration structure. (See
    ``lightning_configure_defaults.pro`` for details and contents.)
``sed_id`` : string array(Nsed)
    The ID of each SED.

Output
------
An error message stating if any errors are found within the tabulated prior structures.
Error messages are relatively specific and should help users resolve any issues.

Modification History
--------------------
- 2022/06/06: Created (Keith Doore)
- 2022/08/01: Added ``/silent`` to ``mrdfits`` (Keith Doore)
- 2022/08/02: Added check to ensure initialization range is within tabulated range (Keith Doore)

