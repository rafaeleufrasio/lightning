.. _converge-good-label:

Convergence Metrics and Goodness of Fit Outputs
===============================================

Depending on your chosen fitting algorithm, the post-processed outputs will vary to include the
relevant convergence metrics and goodness of fit outputs. Below, we give a brief description of
the convergence metrics and goodness of fit outputs relevant to each fitting algorithm. (For an
in-depth discussion, see our page dedicated to explaining these :ref:`convergence-describe-label`.)
Since the MCMC algorithms use most of the same convergence metrics, we split the outputs in terms
of MCMC and MPFIT outputs.

.. note::

    We define additional array size variables here for convenience. 
    (See :ref:`here <array-size-label>` for the already defined array size variables.)

    - ``Nparallel`` : the number of parallel walkers/chains. This value is equal to
      ``NPARALLEL`` in the :ref:`configuration <configure-setting-label>`.
    - ``Nsolvers`` : the number of times MPFIT solved for the best fit SED. This value
      is equal to ``NSOLVERS`` in the :ref:`configuration <configure-setting-label>`.



To help speed up the process of checking if the fitting algorithm converged for each SED, we have
created flags for each convergence metric. The convergence metric flags are set to ``1`` if the 
corresponding convergence metric indicates potentially failed convergence and  ``0`` if it indicates
convergence has been reached. The following output, which is present for all algorithms, combines the
convergence metric flags into a single flag for a simple check of convergence.

``CONVERGENCE_FLAG`` : int
   Logical OR of all other relevant convergence flags of each algorithm. If this flag is ``0`` for
   an SED, we can confidently say that convergence has been reached. If this flag is ``1``, convergence
   may still have been reached, but it is recommended that the flagged convergence metric(s) be inspected.
   (See the :ref:`convergence-describe-label` discussion for more details on how to determine if a
   convergence metric actually indicates failed convergence.)

.. note::

    Below, we refer to convergence having "failed" when not enough walkers, chains, or solvers
    reached the same solution for us to be confident that it is the best solution, or that 
    the parameter space around the solution is well sampled.


MCMC Outputs
------------

``ACCEPTANCE_FRAC`` : double array(Nparallel)
    The acceptance fraction of each MCMC parallel walker/chain (i.e., the fraction 
    of accepted trials out of the total trials).

``ACCEPTANCE_FLAG`` : double array(Nparallel)
    A convergence flag indicating if the acceptance fraction for a given walker/chain
    is below 20% or larger than 50%. Convergence may have failed if the acceptance
    fraction is not within 20-50%.

``AUTOCORR_TIME`` : double array(Nparam)
    The integrated autocorrelation time of each model parameter, averaged over all the MCMC
    walkers/chains. Fixed parameters have a corresponding autocorrelation time of
    ``NaN``.

``BURN_IN_AUTOCORR`` : int
    The burn-in time estimated from the autocorrelation time:

    .. math::

    	{\tt BURN\_IN\_AUTOCORR} = {\rm ceiling}(2 * {\rm max}({\tt AUTOCORR\_TIME}))

``SHORT_CHAIN_FLAG`` : int
    A **non-convergence** flag indicating if ``FINAL_CHAIN_LENGTH`` in the :ref:`configuration <configure-setting-label>`
    is greater than the number of independent samples in the MCMC chain. A short chain can occur
    from any combination of a long burn-in phase, over-thinning of the chain, and/or too few MCMC
    trials. If this flag occurs in your fits, try adjusting one or more of the ``NTRIALS``, ``BURN_IN``, or
    ``THIN_FACTOR`` configuration settings to get your desired ``FINAL_CHAIN_LENGTH``.

``PVALUE`` : double
    The p-value for the fit as determined from a `posterior predictive check (PPC) 
    <https://www.jstor.org/stable/2240995>`_. This p-value gives the goodness of fit,
    which can be used to reject or accept the null hypothesis that the chosen model
    can acceptably model the given SED. Extremely low values may indicate under-fitting
    of the data (i.e. the model is not sufficiently complex to produce the data), while
    extremely high values may indicate over-fitting (i.e. the model is too complex).


Affine-Invariant MCMC Output
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Since we only check the autocorrelation time and remove stranded walkers for the affine-invariant MCMC algorithm,
these outputs are specific to the affine-invariant MCMC:

``AUTOCORR_FLAG`` : double array(Nparam)
    A convergence flag indicating if the autocorrelation time for a given parameter 
    is longer than the length of the chain divided by ``TOLERANCE`` (as given in the
    :ref:`configuration <configure-setting-label>`). Convergence may 
    have failed if the autocorrelation time multiplied by ``TOLERANCE`` is
    larger than the length of the chain.

``STRANDED_FLAG`` : double array(Nparallel)
    A flag indicating which walkers in the ensemble were considered stranded and excluded from
    the post-processed chain portion. (See the :ref:`affine-mcmc-label` description for what
    defines a stranded walker.)



Adaptive MCMC Outputs
^^^^^^^^^^^^^^^^^^^^^

Since the Brooks-Gelman and Gelman-Rubin convergence metrics are not well-defined for the affine-invariant MCMC
algorithm, these outputs are specific to the adaptive MCMC:

``GELMAN_RUBIN_R_HAT`` : double array(Nparam)
    The `Gelman-Rubin <https://ui.adsabs.harvard.edu/abs/1992StaSc...7..457G/abstract>`_ 
    convergence metric (:math:`\hat{R}`) for each parameter.

``GELMAN_RUBIN_FLAG`` : int array(Nparam)
    A convergence flag indicating if the square root of the Gelman-Rubin metric
    is greater than 1.2 for a given parameter. Convergence likely failed 
    if the square root of the metric is larger than 1.2.

``BROOKS_GELMAN_R_HAT`` : double
    The `Brooks-Gelman <https://www.tandfonline.com/doi/abs/10.1080/10618600.1998.10474787>`_
    multidimensional convergence metric (:math:`\hat{R}`).

``BROOKS_GELMAN_FLAG`` : int
    A convergence flag indicating if the square root of the Brooks-Gelman metric
    is greater than 1.2. Convergence may have failed 
    if the square root of the metric is larger than 1.2.

``NCHAIN_R_HAT`` : int
    The number of chains use to compute the Gelman-Rubin and Brooks-Gelman metrics. Chains that got stuck
    in local minima are not used to compute these metrics, as they obviously did not converge.
    Chains are determined to be stuck in local minima if their maximum log probability is 2 less than
    the parallel chain with the overall maximum log probability. The value of 2 is arbitrarily chosen
    but equates to only having a 13.5% chance of being accepted by the MCMC sampler, which is relatively low.


MPFIT Outputs
-------------

``STATUS`` : int array(Nsolvers)
    The status code as returned by MPFIT. See the MPFIT documentation for details
    on each status code.

``STATUS_FLAG`` : int array(Nsolvers)
    A convergence flag indicating if the corresponding MPFIT status code was less than or
    equal to zero, which indicates that the algorithm failed. Convergence of a solver
    definitely failed if its status code is less than or equal to zero.

``ERROR_MSG`` : string array(Nsolvers)
    An error or warning message as given by MPFIT. Typically only given if the algorithm failed.
    If no error message was given by MPFIT, this will be blank.

``ITER_FRAC`` : double array(Nsolvers)
    The fraction of the maximum iterations (``MAXITER`` as given in the
    :ref:`configuration <configure-setting-label>`) used by MPFIT to reach solution.

``ITER_FLAG`` : int array(Nsolvers)
    A convergence flag indicating if the maximum number of iteration were used by MPFIT.
    Convergence likely failed if the maximum number of iterations were used as MPFIT was
    likely still searching for the solution.

``STUCK_FRAC`` : double
    The fraction of solvers that likely got stuck in local minima. Solvers are determined
    to be stuck in local minima if their :math:`\chi^2` is 4 less than the solver with the
    overall minimum :math:`\chi^2` (i.e., best-fit solver). The value of 4 is arbitrarily chosen.

``STUCK_FLAG`` : int
    A convergence flag indicating if the majority of solvers were considered stuck in local
    minima. Convergence may have failed if the majority of solvers got stuck (i.e., did not
    reach a similar :math:`\chi^2` as the best-fit solver).

``PARAMETER_VALUES`` : double array(Nparam, Nsolvers)
    The parameter values for all solvers whose names are given in the ``PARAMETER_NAMES`` output. 
    Useful for comparing with Lightning's default convergence metrics
    (i.e., ``STUCK_FRAC`` and ``SIMILAR_FLAG``).

``SIMILAR_FLAG`` : int
    A convergence flag indicating if any non-stuck solvers had different solutions (i.e., >1% difference
    from best-fit solver's parameter values). Convergence may have failed if a reasonable portion of
    non-stuck solvers resulted in different solutions.

``NFUNC_EVALS`` : int array(Nsolvers)
    The number of ``lightning_mpfit_function.pro`` evaluations performed by MPFIT.

``PVALUE`` : double array(Nsolvers)
    The p-value for the fit as determined from a :math:`\chi^2` test using the :math:`\chi^2` and degrees
    of freedom (``DOF``) as given by MPFIT. This p-value gives a general goodness of fit. However, we caution
    against using it to reject the null hypothesis that the chosen model can acceptably
    model the given SED. Since the effective number of free parameters is lower than the actual
    number (i.e., degeneracies and covariances between parameters reduce the effectiveness), 
    the number of degrees of freedom is likely higher than what is given by MPFIT. Therefore, this
    p-value can be underestimated.

``DOF`` : int
    The number of degrees of freedom as given by MPFIT calculated as:

    .. math::

        {\tt DOF} = N_{\rm data} - N_{\rm param},

    where :math:`N_{\rm data}` is the number of data points and :math:`N_{\rm param}` is the
    number of free parameters (excludes fixed parameters). This value will be the same for each solver,
    since only the starting parameter values are changed for each solver. Also, it will likely
    overestimate the effective degrees of freedom, since paramaters degeneracies and covariances 
    between parameters reduce the effective degrees of freedom. 
