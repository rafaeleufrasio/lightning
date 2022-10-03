LIGHTNING_PRIORS
================

Name
----
LIGHTNING_PRIORS

Purpose
-------
Calculates the log probability of each prior using the input prior
shape arguments. The log probability of each parameter is then 
summed to give the total prior log probability.

Calling Sequence
----------------
::

    lnprior = lightning_priors(parameters, priors [, /error_check, lnprior_param=lnprior_param, $
                               in_bounds_idcs=in_bounds_idcs, out_bounds_idcs=out_bounds_idcs])

Inputs
------
``parameters`` : int, float, or double array(Nparam, Nmodels)
    Parameters of the model(s). The actual parameters contained in this
    array depend on the chosen model(s) during configuration.
``priors`` : structure
    A structure containing the prior hyper-parameters. (See
    ``generate_prior_struct.pro`` for details and contents.)

Optional Input
--------------
``error_check`` : flag
    If set, all inputs are checked for errors. Otherwise, all inputs are
    assumed to be of correct format.

Output
------
``lnprior`` : double array(Nmodels)
    The total prior log probability of each model(s).

Optional Outputs
----------------
``lnprior_param`` : double array(Nparam, Nmodels)
    The prior log probability for each parameter of each model.
``in_bounds_idcs`` : int array(Nmodels)
    The model indices that indicate if the given model(s) had all parameters
    within the specified bounds.
``out_bounds_idcs`` : int array(Nmodels)
    The model indices that indicate if the given model(s) had one or more
    parameters outside the specified bounds.

Notes
-----
`Uniform distribution equation <https://en.wikipedia.org/wiki/Continuous_uniform_distribution>`_ and
`Truncated normal distribution equation <https://en.wikipedia.org/wiki/Truncated_normal_distribution>`_

Modification History
--------------------
- 2022/06/02: Created (Keith Doore)

