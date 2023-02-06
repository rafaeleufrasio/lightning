.. _convergence-describe-label:

Convergence Metrics
===================

.. note::

    We will reference the convergence metrics defined in the :ref:`converge-good-label`.
    So, it may be helpful to review or reference them as you read this discussion.

Convergence for our algorithms is the concept that enough walkers, chains, or solvers have
reached the same solution after their allowed iterations for us to be confident that it is
the best solution and that the parameter space around the solution is well sampled.
The metrics used to test for convergence are therefore essential to ensuring that the
algorithm successfully fit each SED. Below, we give a detailed discussion on the metrics for
each algorithm and how to determine if they indicate that convergence was or was not reached.

Before discussing the individual metrics for each algorithm, we will note that the output
``CONVERGENCE_FLAG`` is a quick and easy way to check if convergence was reached.
If the value is ``0``, the metrics for the given algorithm allow us to confidently say that
convergence has been reached and that we can trust our solution for that SED.
If the value is ``1``, convergence may still have been reached, but it is more ambiguous.
Therefore, the offending metric(s) will need to be inspected to check if it is truly indicating
lack of convergence.


Affine-Invariant MCMC
---------------------

The affine-invariant MCMC algorithm only has two metrics to determine convergence, the acceptance
fraction (``ACCEPTANCE_FRAC``) and the autocorrelation time (``AUTOCORR_TIME``). The acceptance
fraction looks at each individual walker in the ensemble and checks what fraction of the total iterations were
moves accepted for that walker. Depending on the number of free parameters, an acceptance
fraction for each walker is typically expected to be within 20-50%. Lower fractions can indicate that the
algorithm is taking too large of steps in parameter space and failing to properly sample the posterior,
while larger fractions can indicate too small of steps. If any walkers have the ``ACCEPTANCE_FLAG``
set, this does not mean the ensemble as a whole failed to converge. As discussed for the :ref:`affine-mcmc-label`,
if only a few walkers have abnormally low acceptance rates, we label them as stranded and excluded them from
the post-processed chain portion (i.e., they have no effect on convergence). Therefore, we recommend comparing
the ``STRANDED_FLAG`` with the ``ACCEPTANCE_FLAG``, and if they are both set for the same walkers only, then you
can be confident that acceptance fraction metric is **not** indicating failed convergence.

Additionally, in some of our test with more complex models, we have found that most walkers in the ensemble
can have acceptance fractions at or just below 20% when using the default ``AFFINE_A`` in the
:ref:`configure-setting-label`. This does not necessarily indicate failed convergence, especially if the
acceptance fraction is just below 20% (i.e., >18%) and consistent for the vast majority of the ensemble.
However, it does indicate that the algorithm is not sampling efficiently. Therefore, we recommend
rerunning Lightning using a slightly smaller value for ``AFFINE_A``, which should improve sampling
and increase the acceptance fraction into the expected range.

As for the autocorrelation time, it is a measure of how many steps it takes for a walker to "forget"
where it started. We recommend the `emcee Autocorrelation Analysis & Convergence
documentation <https://emcee.readthedocs.io/en/stable/tutorials/autocorr/#autocorr>`_ for more details.
To summarize, the MCMC algorithm needs to run for a number of iterations equal to some factor
(e.g., they recommend ~50) times the autocorrelation time
in order for us to trust that the autocorrelation time estimate is accurate. A factor fewer than ~50 can
cause the autocorrelation time to be underestimated, which could result in a post-processed
chain with a section of highly correlated samples.

You explicitly set what minimum value you are willing to tolerate for this factor
using the ``TOLERANCE`` value in the :ref:`configure-setting-label`, such that
if the ``AUTOCORR_TIME`` times ``TOLERANCE`` is less than ``NTRIALS``, the ``AUTOCORR_FLAG`` will be set.
We recommend using a ``TOLERANCE`` value of 50, since if the factor is above that, you can be confident
that your autocorrelation time is accurate and your walkers converged. However, if a parameter has the
``AUTOCORR_FLAG`` set, we recommend first checking the actual factor (which can be calculated by dividing
``NTRIALS`` by ``AUTOCORR_TIME``) before assuming convergence failed. If this calculated factor is >45,
the autocorrelation time estimate is likely still accurate and convergence was likely still reached.
However, below this we urge caution as ``AUTOCORR_TIME`` can become underestimated and convergence may
have failed. If you are not getting calculated factors large enough to confidently have convergence,
we recommend simply increasing ``NTRIALS`` until you get a reasonable value. (Note that this will very
likely take more additional trials than just ``AUTOCORR_TIME`` times ``TOLERANCE``, as ``AUTOCORR_TIME``
is underestimated and will increase with more trials.)


Adaptive MCMC
-------------

The adaptive MCMC algorithm has three metrics to determine convergence, the acceptance
fraction (``ACCEPTANCE_FRAC``), the Gelman-Rubin convergence metric (``GELMAN_RUBIN_R_HAT``),
and the Brooks-Gelman multidimensional convergence metric (``BROOKS_GELMAN_R_HAT``).
Exactly as the affine-invariant MCMC algorithm, the acceptance fraction
gives the fraction of accepted trials for each parallel chain. However, unlike the affine-invariant MCMC
algorithm, this metric should never be flagged as being outside the expected range of 20-50%.
This is due to the design of the adaptive algorithm. It adapts the proposal distribution to keep
the acceptance fraction between 20-50%. Therefore, if this metric is flagged for the adaptive MCMC algorithm,
you have either run for an extreme number of trials or have too large of a value for ``BETA_EXPONENT``.

The Gelman-Rubin metric is the best indicator of convergence for the adaptive MCMC algorithm.
The metric compares the within chain variance and the between chain variance
to test if multiple parallel chains have converged to the same solution for each parameter. The metric results
in a value that is greater than or equal to 1, where a value of 1 indicates that the parallel chains are identical.
As discussed in their `original paper <https://ui.adsabs.harvard.edu/abs/1992StaSc...7..457G/abstract>`_,
if the square root of this metric is less than 1.2, it can be concluded that the chains have
converged to the same solution. Therefore, if the square root of ``GELMAN_RUBIN_R_HAT`` is greater
than 1.2, the ``GELMAN_RUBIN_FLAG`` is set for the offending parameter. If the flag is set, we recommend
double checking the actual metric to see how close it is to 1.44 (i.e., :math:`1.2^2`). If it is above 1.44 but
below 1.69 (i.e., :math:`1.3^2`), then convergence may still have been reached for the parameter. However,
any value above 1.69 should be considered non-convergence for the parameter, and the SED should be refit.

The Brooks-Gelman multidimension metric is very similar to the Gelman-Rubin metric. The only difference is
that the Brooks-Gelman metric looks at all parameters collectively versus individually like the Gelman-Rubin
metric. Therefore, it is more sensitive to slight differences in the parallel chains. Just like the
``GELMAN_RUBIN_FLAG``, the ``BROOKS_GELMAN_FLAG`` is set if the ``BROOKS_GELMAN_R_HAT`` is larger than 1.44.
If this flag is set and the value of ``BROOKS_GELMAN_R_HAT`` is larger than 1.44, we still highly recommend
checking the ``GELMAN_RUBIN_R_HAT`` first before concluding that convergence has not been reached. In our
tests, we have found several cases where the ``BROOKS_GELMAN_R_HAT`` can be > 2, while all parameters can
have ``GELMAN_RUBIN_R_HAT`` values very close to 1 (i.e., < 1.05). This discrepancy is due to the increased
sensitivity of the Brooks-Gelman metric across the whole chain. Therefore, we recommend relying on the
Gelman-Rubin metric to determine if convergence has failed.


MPFIT
-----

The MPFIT algorithm has four metrics to determine convergence, the status code (``STATUS``), the iteration
fraction (``ITER_FRAC``), the stuck fraction (``STUCK_FRAC``), and the similarity.
The status code is the success status of the MPFIT algorithm. If the code is greater than 0, then the algorithm
executed successfully. This should always be the case when using the MPFIT algorithm in Lightning, since any
errors that could occur in the input or configuration should be detected by Lightning before running. Therefore,
the ``STATUS_FLAG`` will rarely occur. The only time it will is if the random initialization occurs near an
edge of the parameter bounds.

The goal for running multiple solvers for the MPFIT algorithm is the expectation that at least
a majority of them will converge to the same solution. Therefore, the
solvers that did not make it to the same solution need to be filtered out.
The iteration fraction does this by giving
the fraction of the maximum iterations used by each solver to reach their final
solution. If a solver used the maximum allotted iterations, then it was likely still searching
for the best solution before it was terminated by the algorithm. Solvers that reach the maximum
iteration have their ``ITER_FLAG`` set. Therefore, if only a small minority have their
``ITER_FLAG`` set, then convergence of the other solvers may have still occurred. In this case,
we recommend checking the next two metrics to determine if convergence has been reached.

To first check if the solvers reach the same solution, a check needs to be preformed for how many solvers
reached a similar value in :math:`\chi^2` space and how many did not. The fraction of all solvers that did
not reach a similar value of :math:`\chi^2` is given by the stuck fraction. The ``STUCK_FLAG`` is then set
if this fraction is less than 50%, meaning a minority of solvers had similar :math:`\chi^2` values. We define
a similar value of :math:`\chi^2` as values within 4 of the best fit solver. This value is completely arbitrary,
and we note that solvers at :math:`\chi^2` greater than 4 could still have reached the same solution.
Therefore, if this fraction is less than 50%, we recommend comparing the :math:`\chi^2` of each solver to
see how much worse each solver's fit was compared to the best-fit solver.

.. note::

    The :math:`\chi^2` of each solver can be recalculated from ``PVALUE`` of ``DOF`` using the
    IDL function ``CHISQR_CVF`` (i.e., ``chisqr = CHISQR_CVF(PVALUE, DOF)``).

Finally, to check if the solvers reached the same solution in parameter space, the parameter values of the
non-stuck solvers need to be compared for similarity with the best-fit solver. Parameter values that are within 1% difference
of best-fit solver’s parameter values are considered to have converged to the same solution. If parameters
have a larger difference, this can indicate that a multi-modal solution may exist and convergence to a
common solution may not be possible with MPFIT. The ``SIMILAR_FLAG`` is set if any of the non-stuck solvers
had 1% differences in solutions compared to the best-fit solver. If your ``SIMILAR_FLAG`` is set, you
will likely need to check the percent differences in parameter values using the ``PARAMETER_VALUES`` output.
The higher the percent difference the less likely that your solutions converged to the same solution.
Therefore, we recommend not settling for percent differences greater than 5%. Above this fraction, you risk
having a multi-modal solution, which MPFIT is not designed to evaluate.
