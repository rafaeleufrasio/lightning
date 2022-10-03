LIGHTNING_MPFIT_FUNCTION
========================

Name
----
LIGHTNING_MPFIT_FUNCTION

Purpose
-------
Calculates the deviates to be minimized by MPFIT given the data
and a set of parameters.

Calling Sequence
----------------
::

    deviates = lightning_mpfit_function(parameters [, sed_data = , config_nopriors = , $
                                        models = , priors = , /error_check])

Inputs
------
``parameters`` : int, float, or double array(Nparam)
    Parameters of the model. The actual parameters contained in this
    array depend on the chosen model during configuration.

Optional Inputs
---------------
``sed_data`` : structure
    A structure containing the SED luminosities and uncertainties, filter
    labels, distances, redshifts, and optional X-ray data. (See 
    ``lightning_input.pro`` for details and contents.)
``config_nopriors`` : structure
    A Lightning configuration structure edited to remove the prior 
    substructures for each parameter. (See ``lightning_configure_defaults.pro``
    for details and contents.)
``models`` : structure
    A structure containing each model structure (stellar, dust, AGN, 
    X-ray) as a substructure. (See ``lightning_models.pro`` for details
    and contents.)
``priors`` : structure
    A structure containing the prior hyper-parameters. (See
    ``generate_prior_struct.pro`` for details and contents.)
``error_check`` : flag
    If set, all inputs are checked for errors. Otherwise, all inputs are
    assumed to be of correct format.

Output
------
``deviates`` : float or double array(Nfilters + Nxray)
    The deviates of the model for the specified parameters.

Modification History
--------------------
- 2022/08/15: Created (Keith Doore)
- 2022/09/21: Updated to allow for X-ray fluxes (Keith Doore)

