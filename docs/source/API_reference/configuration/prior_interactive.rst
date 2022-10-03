PRIOR_INTERACTIVE
=================

Name
----
PRIOR_INTERACTIVE

Purpose
-------
Generates the prior type, corresponding shape arguments, and
initialization range for a prior structure in the Lightning 
configuration using interactive prompts. It additionally checks
the prompt inputs for errors and ensures proper format.

Calling Sequence
----------------
::

    prior_arg = prior_interactive(param_name, param_descript, param_range, default_prior, $
                                  default_prior_args, default_initialization, prior_options [, $
                                  /edit, prior=prior, initialization_range=initialization_range])

Inputs
------
``param_name`` : string scalar
    The name of the parameter associated with the prior.
``param_descript`` : string scalar
    The description of the parameter.
``param_range`` : int, float, or double array(2)
    The allowed range of the parameter, given in terms of ``[min, max]``.
``default_prior`` : string scalar
    The default prior type choice for the parameter.
``default_prior_args`` : int, float, or double array(5)
    The default shape arguments associated with each prior type, given
    in order of ``[fixed_value, min_bound, max_bound, normal_peak, normal_stddev]``.
``default_initialization`` : int, float, or double array(2)
    The default initialization range of the parameter indicating the minimum
    and maximum bounds for the random initialization of the fitting algorithm.
``prior_options`` : string array(Nprior_options)
    The allowed names of the types of priors that the parameter can be
    set (e.g., ``'uniform'`` or ``'normal'``).

Optional Input
--------------
``edit`` : flag
    If set, then the configuration is being edited and previous values should
    be given rather than defaults.

Output
------
``prior_arg`` : int, float, double, or string array(Narg)
    The distribution shape argument array whose size depends
    on the chosen distribution type.

Optional Outputs
----------------
``prior`` : string scalar
    The chosen prior distribution type.
``initialization_range`` : double array(2)
    The chosen initialization range of the parameter.

Modification History
--------------------
- 2022/04/26: Created (Keith Doore)
- 2022/05/05: Allowed for previous configuration to be printed with ``edit`` keyword (Keith Doore)
- 2022/05/19: Allowed for linear prior inputs to be converted to log-space (Keith Doore)
- 2022/05/31: Removed log priors (Keith Doore)
- 2022/06/27: Updated documentation (Keith Doore)
- 2022/06/27: Updated variable names (Keith Doore)
- 2022/08/01: Added initialization range (Keith Doore)

