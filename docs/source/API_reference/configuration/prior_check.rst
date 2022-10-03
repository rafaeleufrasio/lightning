PRIOR_CHECK
===========

Name
----
PRIOR_CHECK

Purpose
-------
Checks the Lightning configuration prior substructures to make
sure all inputs are correctly formatted and can be understood by
Lightning. If they are not, an error message is printed.

Calling Sequence
----------------
::

    prior_check, prior_struct, param_name, npriors, $
                 param_range, prior_options, prior_options_narg

Inputs
------
``prior_struct`` : structure
    The prior structure for a parameter as generated from the Lightning 
    configuration. The full description of the structure is as follows:

    ====================     ======================================     ==============================
    TAG                      TYPE                                       DESCRIPTION
    ====================     ======================================     ==============================
    PRIOR                    string(Npriors)                            Chosen prior distribution type
    PRIOR_ARG                int/float/double/string(Npriors, Narg)     Distribution shape arguments
    INITIALIZATION_RANGE     int/float/double(Npriors, 2)               Parameter initialization range
    ====================     ======================================     ==============================

``param_name`` : string scalar
    The name of the parameter associated with the prior structure.
``npriors`` : int, float, or double scalar
    The number of priors that should be specified for the parameter.
``param_range`` : int, float, or double array(2)
    The allowed range of the parameter contained within the prior 
    substructure, given in terms of ``[min, max]``.
``prior_options`` : string array(Nprior_options)
    The allowed names of the types of priors that the parameter can be
    set (e.g., ``'uniform'`` or ``'normal'``).
``prior_options_narg`` : int, float, or double array(Nprior_options)
    The number of arguments associated with each with each prior option.

Output
------
An error message is output if any errors are found within the prior substructures.
Error messages are relatively specific and should help users resolve any issues.

Modification History
--------------------
- 2022/04/26: Created (Keith Doore)
- 2022/05/19: Added check to ensure tabulated directory exists (Keith Doore)
- 2022/06/02: Added additional checks to confirm parameters do not result in priors with
  machine over/underflow (Keith Doore)
- 2022/06/27: Renamed inputs for consistent naming scheme (Keith Doore)
- 2022/08/01: Added check for initialization range (Keith Doore)
- 2022/08/01: Fixed bug if input multiple priors and not all had same ``Narg`` (Keith Doore)

