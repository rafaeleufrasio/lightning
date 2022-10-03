PRIOR_SAMPLED_INITIALIZATION
============================

Name
----
PRIOR_SAMPLED_INITIALIZATION

Purpose
-------
Generates the initial parameter values from the prior, restricted to the
initialization range, for sampling with the MCMC or MPFIT algorithms.

Calling Sequence
----------------
::

    parameter_start = prior_sampled_initialization(priors, Nmodels [, sigma_start=sigma_start])

Inputs
------
``priors`` : structure
     A structure containing the prior hyper-parameters and initialization
     range. (See ``generate_prior_struct.pro`` for details and contents.)
``Nmodels`` : int, float, or double scalar
     The number of model parameter sets to initialize.

Output
------
``parameter_start`` : double array(Nparam, Nmodels)
    Initial parameter values for the model(s) as sampled from
    each prior distribution within initialization range.

Optional Output
---------------
``sigma_start`` : double array(Nparam)
    Initial proposal distribution sigma values for each parameter
    as determined from half of the 16th and 84th percentile width
    of each prior distribution.

Modification History
--------------------
- 2022/06/16: Created (Keith Doore)
- 2022/08/02: Updated to use initialization range to initialize within prior (Keith Doore)

