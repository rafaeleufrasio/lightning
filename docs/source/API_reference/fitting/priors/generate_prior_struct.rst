GENERATE_PRIOR_STRUCT
=====================

Name
----
GENERATE_PRIOR_STRUCT

Purpose
-------
Generates the priors structure to be used when computing the
prior probability of each parameter.

Calling Sequence
----------------
::

    priors = generate_prior_struct(config, sed_id [, config_nopriors=config_nopriors])

Inputs
------
``config`` : structure
    A Lightning configuration structure. (See
    ``lightning_configure_defaults.pro`` for details and contents.)
``sed_id`` : string scalar
    A unique SED identifier

Output
------
``priors`` : structure
    This structure includes the hyper-parameters needed to compute the prior
    probability for each parameter in the chosen models. It also includes
    the indices indicating what type of prior each parameter has.
    The full description of the structure is as follows:

    ================     =================     ================================================================================
    TAG                  TYPE                  DESCRIPTION
    ================     =================     ================================================================================
    PARAMETER_NAME       string(Nparam)        Name of each parameter
    IDCS.FIXED           int(...)              Indices of parameters that are fixed (size varies)
    IDCS.VARIABLE        int(...)              Indices of parameters that are not fixed (size varies)
    IDCS.UNIFORM         int(...)              Indices of parameters that have uniform priors (size varies)
    IDCS.NORMAL          int(...)              Indices of parameters that have normal priors (size varies)
    IDCS.ANALYTICAL      int(...)              Indices of parameters that have uniform or normal priors (size varies)
    IDCS.TABULATED       int(...)              Indices of parameters that have tabulated priors (size varies)
    INITIALIZE_RANGE     double(Nparam, 2)     The lower and upper bounds for the initialization range of the fitting algorithm
    FIXED                double(Nparam)        Values of fixed parameters (``NaN`` if not fixed)
    MIN_BOUND            double(Nparam)        Values of the minimum bounds for the priors
    MAX_BOUND            double(Nparam)        Values of the maximum bounds for the priors
    ANALYTICAL.WIDTH     double(Nparam)        Values of the width for the normal priors (:math:`-0.5/\sigma^2`)
    ANALYTICAL.MU        double(Nparam)        Values of the peak value for the normal priors
    ANALYTICAL.CONST     double(Nparam)        Values of the normalization constant for the uniform and normal priors
    TABULATED            structure             A structure containing the tabulated values for the tabulated priors
    ================     =================     ================================================================================

Optional Output
---------------
``config_nopriors`` : structure
    The input Lightning configuration structure edited to remove the prior 
    substructures for each parameter.

Modification History
--------------------
- 2022/06/06: Created (Keith Doore)
- 2022/06/20: Fixed issue with tabulated prior indexing of hash keys (Keith Doore)
- 2022/07/05: Allowed for other parameters with multiple priors besides ``PSI`` (Keith Doore)
- 2022/07/05: Removed parameter substructures from ``config`` after extracting values (Keith Doore)
- 2022/08/01: Added ``/silent`` to ``mrdfits`` (Keith Doore)
- 2022/08/02: Added ``initialization_range`` as ``initialize_range`` to prior structure (Keith Doore)

